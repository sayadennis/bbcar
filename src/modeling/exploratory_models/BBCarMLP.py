import os
import sys
import numpy as np
import pandas as pd
import torch
from torch import nn
from collections import  defaultdict
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, balanced_accuracy_score
import utils

class BBCarMLP(torch.nn.Module):
    def __init__(self, X, Xval, y = None, yval=None,
                 n_iter = 10, eps = 1e-7,
                 weight_decay = 1e-5, wortho=1, C=1,
                 lr = 1e-2, verbose = False, seed=None, fn=None,
                 device=torch.device('cpu')):
        super(BBCarMLP, self).__init__()
        if seed is not None:
            torch.manual_seed(seed)
            np.random.seed(seed)
        self.n_iter = n_iter
        self.weight_decay = weight_decay
        self.lr = lr
        self.verbose = verbose
        self.eps = eps
        self.fn = fn
        self.C = C
        self.device = device
        self.celoss = -1
        self.report = defaultdict(list)
        self.__initfact__(X, Xval, y, yval)

    def __initfact__(self, X, Xval, y, yval):
        self.n = torch.tensor(X.shape, dtype=int)[0].to(self.device)
        self.m = torch.tensor(X.shape, dtype=int)[1].to(self.device)
        self.ntr = torch.tensor(len(y), dtype=int).to(self.device)
        self.nval = torch.tensor(len(yval), dtype=int).to(self.device)
        self.nte = self.n - self.ntr - self.nval
        self.X = torch.from_numpy(X).float().to(self.device)
        self.Xval = torch.from_numpy(Xval).float().to(self.device)

        self.y = torch.from_numpy(y).long().to(self.device)
        w = 1 / pd.Series(y).value_counts(normalize=True).sort_index().to_numpy()
        self.loss_cls = nn.CrossEntropyLoss(weight=torch.from_numpy(w).float()).to(self.device) # weight=torch.from_numpy(w).float()
        
        self.fc = nn.Linear(self.m, len(np.unique(y))).to(self.device)
        
        self.yval = torch.from_numpy(yval).long().to(self.device)
        self.opt = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

    def to(self,device):
        self.device = device
        self.X = self.X.to(device)
        self.y = self.y.to(device)
        self.yval = self.yval.to(device)        
        return super(BBCarMLP, self).to(device)

    def __autograd__(self,epoch):
        """
           autograd update, with gradient projection
        """
        self.opt.zero_grad()
        if self.y is not None:
            self.celoss = self.loss_cls(self.fc(self.X), self.y)
            l = self.celoss * self.n # 01/24/20
            # print('cross entropy: %.4f' % (self.loss_cls(self.fc(self.Wtr), self.y)))
            for p in self.fc.parameters():
                l = l + p.pow(2).sum() * self.C 
                # print('complexity: %.4f' % (p.pow(2).sum()))

        l.backward()
        self.opt.step()
        return l.item()

    def __update__(self,epoch):
        l = self.__autograd__(epoch)

        self.report['epoch'].append(epoch)
        self.report['loss'].append(l)
        if self.verbose and epoch % 100 == 0:
            print("%d\tloss: %.4f"%(epoch,l))

    def fit(self):
        it = range(self.n_iter)
        # for autograd solver
        best_val_acc = 0
        for e in it:
            self.__update__(e)
            # here using pinverse seems to mess up GPU/CPU, it's really the pinverse that's taking a lot of CPU.
            if e >= self.n_iter - 1:
                y_val_pred = self.__predict__(self.Xval)
                acc = balanced_accuracy_score(self.yval.detach().cpu().numpy(),
                                    y_val_pred.detach().cpu().numpy())
                if acc > best_val_acc:
                    best_val_acc = acc
                    if self.fn is not None:
                        torch.save({'epoch': e+1,
                                    'state_dict': self.state_dict(),
                                    'optimizer': self.opt.state_dict(),
                                    'best_val_acc': best_val_acc,
                                    'report': self.report,
                                    'celoss': self.celoss,
                        }, self.fn)
                        # print('best_val_acc: %.4f, test_acc: %.4f' % (best_val_acc, test_acc))
        self.decomposed = True
        return self

    def show_report(self):
        return pd.DataFrame(self.report)

    def predict(self, X):
        X = torch.from_numpy(X).float()
        if self.device.type == 'cuda':
            X = X.to(self.device)
        y = self.__predict__(X)
        if self.device.type == 'cuda':
            return y.detach().cpu().numpy()
        else:
            return y.detach().numpy()

    def __predict__(self, X):
        y = torch.argmax(self.fc(X), dim=1)
        return y

    # def __init__(self, input_size, hidden_size):
    #     super(BBCarMLP, self).__init__()
    #     self.input_size = input_size
    #     self.hidden_size  = hidden_size
    #     self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
    #     self.relu = torch.nn.ReLU()
    #     self.fc2 = torch.nn.Linear(self.hidden_size, 1)
    #     self.sigmoid = torch.nn.Sigmoid()
    #     ## alternative A 
    #     # self.layers = nn.Sequential(
    #     #     nn.Flatten(),
    #     #     nn.Linear(32 * 32 * 3, 64),
    #     #     nn.ReLU(),
    #     #     nn.Linear(64, 32),
    #     #     nn.ReLU(),
    #     #     nn.Linear(32, 10)
    #     #     )
    #     ## alternative B
    #     # self.lin1 = nn.Linear(784, 512, bias=True) 
    #     # self.lin2 = nn.Linear(512, 256, bias=True)
    #     # self.lin3 = nn.Linear(256, 10, bias=True)
    #     ## alternative C
    #     # self.dropout = nn.Dropout(0.2)


    
    # def forward(self, x):
    #     hidden = self.fc1(x)
    #     relu = self.relu(hidden)
    #     output = self.fc2(relu)
    #     output = self.sigmoid(output)
    #     ## alternative A 
    #     # return self.layers(x)
    #     ## alternative B
    #     # x = xb.view(-1,784) 
    #     # x = F.relu(self.lin1(x))
    #     # x = F.relu(self.lin2(x))
    #     # return self.lin3(x)
    #     return output

