import numpy as np
import pandas as pd
import torch
from torch import nn
from torch import optim
from sklearn import metrics

k = 5 # number of components for hNMF
loss_type = ['l2', 'kl'][0]
weight_decay = 1e-5
eps = 1e-7
orthogonality_weight = 1
clf_weight = 0.1
C = 1. # regularization weight for the classification fully connected layers

class suphNMF(nn.Module):
    def __init__(
            self, X1, X2, y, n_components=k, lr=1e-3, n_iter=2000, 
            weight_decay=1e-4, clf_weight=0.1, ortho_weight=1., seed=42):
        super(suphNMF, self).__init__()
        torch.manual_seed(seed)
        np.random.seed(seed)
        self.factorized = False
        # Data
        self.X1 = torch.tensor(X1.values, dtype=torch.float32)
        self.X2 = torch.tensor(X2.values, dtype=torch.float32)
        self.y = torch.tensor(y.values, dtype=torch.float32)
        # Parameters and number of components
        self.n_components = n_components
        self.W = torch.normal(mean=3., std=1.5, size=(self.X1.shape[0], self.n_components)).clamp(min=1e-5)
        self.H1 = torch.normal(mean=3., std=1.5, size=(self.n_components, self.X1.shape[1])).clamp(min=1e-5)
        self.H2 = torch.normal(mean=3., std=1.5, size=(self.n_components, self.X2.shape[1])).clamp(min=1e-5)
        self.fc = nn.Linear(in_features=self.n_components, out_features=1)
        # Optimizer and loss functions
        self.n_iter = n_iter
        self.recon_loss_func = nn.MSELoss(reduction='mean')
        self.clf_loss_func = nn.BCEWithLogitsLoss(reduction='mean')
        self.lr = lr
        self.weight_decay = weight_decay
        self.clf_weight = clf_weight
        self.ortho_weight = ortho_weight

    def plus(self):
        self.W.data = self.W.data.clamp(min=1e-5)
        self.H1.data = self.H1.data.clamp(min=1e-5)
        self.H2.data = self.H2.data.clamp(min=1e-5)

    def fit(self):
        self.W = torch.nn.Parameter(self.W)
        self.H1 = torch.nn.Parameter(self.H1)
        self.H2 = torch.nn.Parameter(self.H2)
        self.optimizer = optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        self.clf_roc_record = []
        self.clf_loss_record = []
        self.reconloss1_record = []
        self.reconloss2_record = []
        self.ortholoss_record = []
        self.loss_record = []
        for i in range(self.n_iter):
            ## Gradient calculation
            # Zero out gradient
            self.optimizer.zero_grad()
            # Calculate factorization (reconstruction) loss (W * H - X)
            self.loss_recon1 = self.recon_loss_func(torch.mm(self.W, self.H1), self.X1)
            self.loss_recon2 = self.recon_loss_func(torch.mm(self.W, self.H2), self.X2)
            # Add L2 (or KL) regularization for orthogonality between components
            WtW = torch.mm(torch.t(self.W), self.W)
            if torch.mean(WtW) > 1e-7: # self.eps 
                WtW = WtW / torch.mean(WtW)
            self.loss_W_ortho = self.recon_loss_func(WtW/self.n_components, torch.eye(self.n_components)) * self.ortho_weight * self.n_components
            # Add classification loss
            y_pred = self.predict(self.W)
            self.loss_clf = self.clf_loss_func(y_pred, self.y)
            # No need to add L2 regularization to the clf params since added weight decay
            ## Backprop & step
            loss = self.loss_recon1 + self.loss_recon2 + self.clf_weight * self.loss_clf + self.loss_W_ortho
            loss.backward()
            self.optimizer.step()
            self.loss_record.append(loss.item())
            ## Set factorized matrix values to positive
            self.plus()
            ## Record performance
            self.clf_roc_record.append(metrics.roc_auc_score(self.y.detach().numpy(), y_pred.detach().numpy()))
            self.clf_loss_record.append(self.loss_clf.item())
            self.reconloss1_record.append(self.loss_recon1.item())
            self.reconloss2_record.append(self.loss_recon2.item())
            self.ortholoss_record.append(self.loss_W_ortho.item())
        ## Save: model object? - W, H, classifier, performances, loss curves etc.
        self.factorized = True
   
    def predict(self, W):
        """
        Make predictions using the W matrix.
        """
        return self.fc(W)

    def transform(self, X1_tf, X2_tf, y_tf=None):
        """
        Transform a given input X1/X2 into W, using the learned H1/H2.
        """
        if not self.factorized:
            print("Warning: fit has not been called yet. Using randomly initialized H.")

        self.W_tf = torch.nn.Parameter(
            torch.normal(mean=3., std=1.5, size=(X1_tf.shape[0], self.n_components)).clamp(min=1e-5)
        )
        X1_tf = torch.tensor(np.array(X1_tf), dtype=torch.float32)
        X2_tf = torch.tensor(np.array(X2_tf), dtype=torch.float32)
        if y_tf is not None:
            y_tf = torch.tensor(np.array(y_tf), dtype=torch.float32)
        self.optimizer_tf = optim.Adam([self.W_tf], lr=self.lr, weight_decay=self.weight_decay)
        self.clf_roc_record_tf = []
        self.clf_loss_record_tf = []
        self.reconloss1_record_tf = []
        self.reconloss2_record_tf = []
        self.ortholoss_record_tf = []
        self.loss_record_tf = []
        
        for i in range(self.n_iter):
            ## Gradient calculation
            # Zero out gradient
            self.optimizer_tf.zero_grad()
            # Calculate factorization (reconstruction) loss (W * H - X)
            loss_recon1_tf = self.recon_loss_func(torch.mm(self.W_tf, self.H1), X1_tf)
            loss_recon2_tf = self.recon_loss_func(torch.mm(self.W_tf, self.H2), X2_tf)
            # Add L2 (or KL) regularization for orthogonality between components
            WtW = torch.mm(torch.t(self.W_tf), self.W_tf)
            if torch.mean(WtW) > 1e-7: # self.eps 
                WtW = WtW / torch.mean(WtW)
            loss_W_ortho_tf = self.recon_loss_func(WtW/self.n_components, torch.eye(self.n_components)) * self.n_components
            # Add classification loss
            if y_tf is not None:
                y_pred = self.predict(self.W_tf)
                loss_clf_tf = self.clf_loss_func(y_pred, y_tf)
            else:
                loss_clf_tf = 0
            # No need to add L2 regularization to the clf params since added weight decay
            ## Backprop & step
            loss = loss_recon1_tf + loss_recon2_tf + self.ortho_weight * loss_W_ortho_tf
            loss.backward()
            self.optimizer_tf.step()
            self.loss_record_tf.append(loss.item())
            ## Set factorized matrix values to positive
            self.W_tf.data = self.W_tf.data.clamp(min=1e-5)
            ## Record performance
            if y_tf is not None:
                self.clf_roc_record_tf.append(metrics.roc_auc_score(y_tf, y_pred.detach().numpy()))
            self.clf_loss_record_tf.append(loss_clf_tf.item())
            self.reconloss1_record_tf.append(loss_recon1_tf.item())
            self.reconloss2_record_tf.append(loss_recon2_tf.item())
            self.ortholoss_record_tf.append(loss_W_ortho_tf.item())
        ## Save: model object? - W, H, classifier, performances, loss curves etc.
        return self.W_tf

