from itertools import product

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch import optim
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold

k = 5 # number of components for hNMF
loss_type = ['l2', 'kl'][0]
weight_decay = 1e-5
eps = 1e-7
orthogonality_weight = 1
clf_weight = 0.1
C = 1. # regularization weight for the classification fully connected layers

class suphNMF(nn.Module):
    def __init__(
            self, X1, X2, y=None, n_components=k, seed=42,
            n_iters=[1000,2000,3000], lrs=[1e-4,1e-3,1e-2], weight_decays=[1e-3, 1e-4, 1e-5], clf_weights=[0,1e-1,1e+0,1e+1,1e+2], ortho_weights=[0,1e-1,1e+0,1e+1,1e+2]):
        super(suphNMF, self).__init__()
        torch.manual_seed(seed)
        np.random.seed(seed)
        self.seed = seed
        self.factorized = False

        # Data
        self.X1 = torch.tensor(X1.values, dtype=torch.float32)
        self.X2 = torch.tensor(X2.values, dtype=torch.float32)
        if y is not None:
            self.y = torch.tensor(y.values, dtype=torch.float32)
            self.supervised = True
        else:
            self.supervised = False

        # Parameters and number of components
        self.n_components = n_components
        self.W = torch.normal(mean=3., std=1.5, size=(self.X1.shape[0], self.n_components)).clamp(min=1e-5)
        self.H1 = torch.normal(mean=3., std=1.5, size=(self.n_components, self.X1.shape[1])).clamp(min=1e-5)
        self.H2 = torch.normal(mean=3., std=1.5, size=(self.n_components, self.X2.shape[1])).clamp(min=1e-5)
        self.fc = nn.Linear(in_features=self.n_components, out_features=1)

        # Optimizer and loss functions
        self.n_iters = n_iters
        self.recon_loss_func = nn.MSELoss(reduction='mean')
        self.clf_loss_func = nn.BCEWithLogitsLoss(reduction='mean')
        self.lrs = lrs
        self.weight_decays = weight_decays
        self.clf_weights = clf_weights
        self.ortho_weights = ortho_weights

    def plus(self):
        self.W.data = self.W.data.clamp(min=1e-5)
        self.H1.data = self.H1.data.clamp(min=1e-5)
        self.H2.data = self.H2.data.clamp(min=1e-5)

    def fit_cv(self):
        ## Record CV performance
        cv_record = {}
        param_tuples = list(product(self.n_iters, self.lrs, self.weight_decays, self.clf_weights, self.ortho_weights))
        self.cv_performance = pd.DataFrame(
            index=pd.MultiIndex.from_tuples(param_tuples, names=('n_iter', 'lr', 'weight_decay', 'w_clf', 'w_ortho')), 
            columns=['performance']
        )
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.seed)
        for n_iter, lr, weight_decay, w_clf, w_ortho in param_tuples:
            fold_performances = []
            for i, (train_index, val_index) in enumerate(skf.split(self.X1, self.y)):
                X1_train, X2_train, y_train = self.X1[train_index,:], self.X2[train_index,:], self.y[train_index]
                X1_val, X2_val, y_val = self.X1[val_index,:], self.X2[val_index,:], self.y[val_index]
                out = self._fit_cv(X1_train, X2_train, y_train, X1_val, X2_val, y_val, n_iter=n_iter, lr=lr, weight_decay=weight_decay, clf_weight=w_clf, ortho_weight=w_ortho)
                cv_record[i] = out
                # for this given fold and hyperparameters, find the validation performance (ignoring number of iterations for now) 
                fold_performances.append(out['clf_val_roc'][-1])
            
            # calculate mean performance for these hyperparameters (ignoring number of iterations for now)
            self.cv_performance.loc[(n_iter, lr, weight_decay, w_clf, w_ortho),'performance'] = np.mean(fold_performances)
            self.cv_performance.to_csv('/home/srd6051/test_cv_performance.csv')

        ## Find the hyperparameters and n_iter that gave the best mean performance across the CV folds
        self.cv_performance = self.cv_performance.astype(float)
        best_params_tuple = self.cv_performance.idxmax()
        best_params = {'n_iter' : best_params_tuple[0][0], 'lr' : best_params_tuple[0][1], 'weight_decay' : best_params_tuple[0][2], 'w_clf' : best_params_tuple[0][3], 'w_ortho' : best_params_tuple[0][4]}

        return cv_record, best_params

    def _fit_cv(self, X1_train, X2_train, y_train, X1_val, X2_val, y_val, n_iter, lr, weight_decay, clf_weight, ortho_weight):
        # Hyperparams
        self.n_iter = n_iter
        self.lr = lr
        self.weight_decay = weight_decay
        self.clf_weight = clf_weight
        self.ortho_weight = ortho_weight

        # NMF matrices
        W = torch.nn.Parameter(torch.normal(mean=3., std=1.5, size=(X1_train.shape[0], self.n_components)).clamp(min=1e-5))
        H1 = torch.nn.Parameter(torch.normal(mean=3., std=1.5, size=(self.n_components, X1_train.shape[1])).clamp(min=1e-5))
        H2 = torch.nn.Parameter(torch.normal(mean=3., std=1.5, size=(self.n_components, X2_train.shape[1])).clamp(min=1e-5))

        self.optimizer = optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

        # Observed data
        self.X1_train_tr = X1_train
        self.X2_train_tr = X2_train
        self.X1_train_val = X1_val
        self.X2_train_val = X2_val

        # training record
        clf_roc_record = []
        clf_loss_record = []
        reconloss1_record = []
        reconloss2_record = []
        ortholoss_record = []
        loss_record = []

        # validation
        clf_val_roc_record = []
        clf_val_loss_record = []
        reconloss1_val_record = []
        reconloss2_val_record = []
        ortholoss_val_record = []
        loss_val_record = []

        for i in range(self.n_iter):
            ## Gradient calculation
            # Zero out gradient
            self.optimizer.zero_grad()
            # Calculate factorization (reconstruction) loss (W * H - X)
            loss_recon1 = self.recon_loss_func(torch.mm(W, H1), X1_train) / X1_train.shape[1]  # divide by number of original feature space (align scale between mutation and CNV) 
            loss_recon2 = self.recon_loss_func(torch.mm(W, H2), X2_train) / X2_train.shape[1]
            # Add L2 (or KL) regularization for orthogonality between components
            WtW = torch.mm(torch.t(W), W)
            if torch.mean(WtW) > 1e-7: # self.eps 
                WtW = WtW / torch.mean(WtW)
            loss_W_ortho = self.recon_loss_func(WtW/self.n_components, torch.eye(self.n_components)) * self.n_components
            # Add classification loss
            if self.supervised:
                y_pred = self.predict(W)
                loss_clf = self.clf_loss_func(y_pred, y_train)
            else:
                self.loss_clf = torch.tensor(0.)
            # No need to add L2 regularization to the clf params since added weight decay
            ## Backprop & step
            loss = loss_recon1 + loss_recon2 + self.clf_weight * loss_clf + self.ortho_weight * loss_W_ortho
            loss.backward()
            self.optimizer.step()
            loss_record.append(loss.item())
            ## Set factorized matrix values to positive
            self.plus()
            ## Record performance
            reconloss1_record.append(loss_recon1.item())
            reconloss2_record.append(loss_recon2.item())
            ortholoss_record.append(loss_W_ortho.item())
            if self.supervised:
                clf_roc_record.append(metrics.roc_auc_score(y_train.detach().numpy(), y_pred.detach().numpy()))
                clf_loss_record.append(loss_clf.item())

        # Record performance on validation set
        W_val = self.transform(X1_val, X2_val, y_val)
        loss_val_recon1 = self.recon_loss_func(torch.mm(W_val, H1), X1_val) / X1_val.shape[1]
        loss_val_recon2 = self.recon_loss_func(torch.mm(W_val, H2), X2_val) / X2_val.shape[1]
        WtW = torch.mm(torch.t(W_val), W_val)
        if torch.mean(WtW) > 1e-7: # self.eps 
            WtW = WtW / torch.mean(WtW)
        loss_val_W_ortho = self.recon_loss_func(WtW/self.n_components, torch.eye(self.n_components)) * self.n_components
        if self.supervised:
            y_val_pred = self.predict(W_val)
            loss_val_clf = self.clf_loss_func(y_val_pred, y_val)
            clf_val_roc = metrics.roc_auc_score(y_train.detach().numpy(), y_pred.detach().numpy())

        loss_val_overall = loss_val_recon1 + loss_val_recon2 + self.clf_weight * loss_val_clf + self.ortho_weight * loss_val_W_ortho

        clf_val_roc_record.append(clf_val_roc)
        clf_val_loss_record.append(loss_val_clf.item())
        reconloss1_val_record.append(loss_val_recon1.item())
        reconloss2_val_record.append(loss_val_recon2.item())
        ortholoss_val_record.append(loss_val_W_ortho.item())
        loss_val_record.append(loss_val_overall.item())

        return {
            "clf_roc" : clf_roc_record, "clf_loss" : clf_loss_record, "reconloss1" : reconloss1_record, "reconloss2" : reconloss2_record, "ortholoss" : ortholoss_record, "overall" : loss_record,
            "clf_val_roc" : clf_val_roc_record, "clf_val_loss" : clf_val_loss_record, "reconloss1_val" : reconloss1_val_record, "reconloss2_val" : reconloss2_val_record, "ortholoss_val" : ortholoss_val_record, "overall" : loss_val_record,
        }

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
            self.loss_recon1 = self.recon_loss_func(torch.mm(self.W, self.H1), self.X1) / self.X1.shape[1]  # divide by number of original feature space (align scale between mutation and CNV) 
            self.loss_recon2 = self.recon_loss_func(torch.mm(self.W, self.H2), self.X2) / self.X2.shape[1]
            # Add L2 (or KL) regularization for orthogonality between components
            WtW = torch.mm(torch.t(self.W), self.W)
            if torch.mean(WtW) > 1e-7: # self.eps 
                WtW = WtW / torch.mean(WtW)
            self.loss_W_ortho = self.recon_loss_func(WtW/self.n_components, torch.eye(self.n_components)) * self.n_components
            # Add classification loss
            if self.supervised:
                y_pred = self.predict(self.W)
                self.loss_clf = self.clf_loss_func(y_pred, self.y)
            else:
                self.loss_clf = torch.tensor(0.)
            # No need to add L2 regularization to the clf params since added weight decay
            ## Backprop & step
            loss = self.loss_recon1 + self.loss_recon2 + self.clf_weight * self.loss_clf + self.ortho_weight * self.loss_W_ortho
            loss.backward()
            self.optimizer.step()
            self.loss_record.append(loss.item())
            ## Set factorized matrix values to positive
            self.plus()
            ## Record performance
            self.reconloss1_record.append(self.loss_recon1.item())
            self.reconloss2_record.append(self.loss_recon2.item())
            self.ortholoss_record.append(self.loss_W_ortho.item())
            if self.supervised:
                self.clf_roc_record.append(metrics.roc_auc_score(self.y.detach().numpy(), y_pred.detach().numpy()))
                self.clf_loss_record.append(self.loss_clf.item())
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
        #if not self.factorized:
        #    print("Warning: fit has not been called yet. Using randomly initialized H.")

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
                loss_clf_tf = torch.tensor(0.)
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

