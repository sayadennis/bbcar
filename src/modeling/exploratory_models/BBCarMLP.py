import os
import sys
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn.functional import relu
from collections import  defaultdict
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, balanced_accuracy_score
import utils

class BBCarMLP(torch.nn.Module):
    def __init__(self, X, Xval, y = None, yval=None,
                 n_iter = 10, eps = 1e-7, hidden_layers=None,
                 weight_decay = 1e-5, C=1, contrastive=False, wcont=1e-2, # weight of contrastive loss
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
        self.hidden_layers = hidden_layers
        self.device = device
        self.celoss = -1
        self.contrastive = contrastive
        self.wcont = wcont
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

        if self.hidden_layers==3:
            self.h1_size = 400
            self.h2_size = 200
            self.h3_size = 20 # try 40 next time 
            self.fc1 = nn.Linear(self.m, self.h1_size).to(self.device)
            self.fc2 = nn.Linear(self.h1_size, self.h2_size).to(self.device)
            self.fc3 = nn.Linear(self.h2_size, self.h3_size).to(self.device)
            self.fc4 = nn.Linear(self.h3_size, len(np.unique(y))).to(self.device)
        elif self.hidden_layers==2:
            self.h1_size = 400
            self.h2_size = 40
            self.fc1 = nn.Linear(self.m, self.h1_size).to(self.device)
            self.fc2 = nn.Linear(self.h1_size, self.h2_size).to(self.device)
            self.fc3 = nn.Linear(self.h2_size, len(np.unique(y))).to(self.device)
        else:
            KeyError('Unknown hidden layer size: %s' % self.hidden_layers)
        
        self.yval = torch.from_numpy(yval).long().to(self.device)
        self.opt = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)
    
    def forward(self, X):
        if self.hidden_layers==3:
            out = relu(self.fc1(X))
            out = relu(self.fc2(out))
            out = relu(self.fc3(out))
            self.finallayer = self.fc3(out)
            out = torch.softmax(self.fc4(out), dim=1)
        elif self.hidden_layers==2:
            out = relu(self.fc1(X))
            out = relu(self.fc2(out))
            self.finallayer = self.fc3(out)
            out = torch.softmax(self.fc3(out), dim=1)
        return out

    def to(self,device):
        self.device = device
        self.X = self.X.to(device)
        self.y = self.y.to(device)
        self.yval = self.yval.to(device)        
        return super(BBCarMLP, self).to(device)
    
    def get_pairs(self,ix): # get non-overlapping pairs to calculate distances
        pairs = []
        for i in ix:
            for j in ix:
                if ((i,j) not in pairs) and ((j,i) not in pairs):
                    pairs.append((i,j))
        return pairs

    def __autograd__(self,epoch):
        """
           autograd update, with gradient projection
        """
        self.opt.zero_grad()
        if self.y is not None:
            self.celoss = self.loss_cls(self.forward(self.X), self.y)
            l = self.celoss * self.n # 01/24/20
            # add contrastive loss
            if self.contrastive=='euclidean':
                contloss = 0
                # use pairs of self.finallayer and y 
                for (i,j) in self.get_pairs(range(self.X.shape[0])):
                    dist = abs(self.y[i]-self.y[j]) # calculate distance (0 or 1)
                    if dist==0: # if same label
                        contloss += torch.cdist(torch.reshape(self.finallayer[i,],(1,-1)), torch.reshape(self.finallayer[j,],(1,-1)))[0][0]
                    elif dist==1:
                        contloss -= torch.cdist(torch.reshape(self.finallayer[i,],(1,-1)), torch.reshape(self.finallayer[j,],(1,-1)))[0][0]
                l += self.wcont * contloss
            # print('cross entropy: %.4f' % (self.loss_cls(self.fc(self.Wtr), self.y)))
            if self.hidden_layers==3:
                for fc in [self.fc1, self.fc2, self.fc3, self.fc4]:
                    for p in fc.parameters():
                        l = l + p.pow(2).sum() * self.C 
            if self.hidden_layers==2:
                for fc in [self.fc1, self.fc2, self.fc3]:
                    for p in fc.parameters():
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
        y = torch.argmax(self.forward(X), dim=1)
        return y
