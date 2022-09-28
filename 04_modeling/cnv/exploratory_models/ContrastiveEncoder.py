import numpy as np
import pandas as pd
import torch
from torch import nn
import utils
from collections import defaultdict
from sklearn.metrics import accuracy_score, balanced_accuracy_score

## this implements ScanMap transductive learning with demographics information and orthogonality constraints
## assumes X is passed as [Xtrain; Xtest], and only y_train is passed
class ContrastiveEncoder(nn.Module):
    def __init__(self, Xtr, Xval, ytr, yval, n_iter=10, eps=1e-7,
                 floss='l2', weight_decay=1e-5, wortho=1, C=1,
                 lr=1e-2, verbose=False, seed=None, fn=None,
                 device=torch.device('cpu')):
        super(ContrastiveEncoder, self).__init__()
        # define parameters that are plugged in here 
        if seed is not None:
            torch.manual_seed(seed)
            np.random.seed(seed)
        self.n_iter = n_iter
        # self.k = k
        self.floss = floss
        self.weight_decay = weight_decay
        self.lr = lr
        self.verbose = verbose
        self.eps = eps
        # self.wcls = wcls # weight of classification loss
        self.wortho = wortho
        self.fn = fn
        self.C = C
        self.fitted = False
        # self.dcf = cf.shape[1]
        self.device = device
        self.celoss = -1
        self.report = defaultdict(list)
        self.__initfact__(Xtr, Xval, ytr, yval)

    def __initfact__(self, Xtr, Xval, ytr, yval):
        self.ntr = torch.tensor(len(ytr), dtype=int)[0].to(self.device)
        self.nval = torch.tensor(len(yval), dtype=int).to(self.device)
        self.m = torch.tensor(Xtr.shape, dtype=int)[1].to(self.device)
        self.Xtr = torch.from_numpy(Xtr).float().to(self.device)
        
        if self.floss == 'l2':
            self.loss_fac = utils.l2
        elif self.floss == 'kl':
            self.loss_fac = utils.kl_div

        self.ytr = torch.from_numpy(ytr).long().to(self.device)
        self.yval = torch.from_numpy(yval).long().to(self.device)
        w = 1 / pd.Series(ytr).value_counts(normalize=True).sort_index().to_numpy()
        self.loss_cls = nn.CrossEntropyLoss(weight=torch.from_numpy(w).float()).to(self.device) # weight=torch.from_numpy(w).float()

        ## Insert architecture here ## (i think)
        self.fc = nn.Linear(self.m, len(torch.unique(self.y))).to(self.device) # 4 = len(torch.unique(self.y))**2
        self.opt = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

    def to(self,device):
        self.device = device
        self.Xtr = self.Xtr.to(device)
        self.Xval = self.Xval.to(device)
        self.y = self.y.to(device)
        self.yval = self.yval.to(device)        
        return super(ContrastiveEncoder, self).to(device)

    def plus(self,X):
        X[X < 0] = 0 # self.eps
        return X

    def __autograd__(self,epoch):
        """
           autograd update, with gradient projection
        """
        self.opt.zero_grad()
        ## add classification loss 
        self.celoss = self.loss_cls(self.fc(Xtr), self.y)
        l = l + self.celoss * self.n # 01/24/20 # * self.wcls 
        # print('cross entropy: %.4f' % (self.loss_cls(self.fc(self.Wtr), self.y)))
        ## add L2 regularization 
        for p in self.fc.parameters():
            # has two parameters (cls #, ft #) and (cls #), they are the weight and biases, /self.k to add normalized weight regularization
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
                y_val_pred = self.__predict__(Xval)
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
        self.fitted = True
        return self

    def show_report(self):
        return pd.DataFrame(self.report)

    def fit_transform(self): # need edit: transform not applicable for this case?
        if not self.fitted:
            self.fit()
            # detach all params including the linear fc layer
            for p in self.parameters():
                p.requires_grad = False
        return
        # if self.device.type == 'cuda':
        #     return [self.Wtr.detach().cpu().numpy(), self.Wval.detach().cpu().numpy(), self.Wte.detach().cpu().numpy(), self.H.detach().cpu().numpy()]
        # else:
        #     return [self.Wtr.detach().numpy(), self.Wval.detach().numpy(), self.Wte.detach().numpy(), self.H.detach().numpy()]

    def predict(self, X):
        X = torch.from_numpy(X).float()
        # cf = torch.from_numpy(cf).float()                        
        if self.device.type == 'cuda':
            X = X.to(self.device)
            # cf = cf.to(self.device)
        y = self.__predict__(X)
        if self.device.type == 'cuda':
            return y.detach().cpu().numpy()
        else:
            return y.detach().numpy()

    def __predict__(self, X):
        y = torch.argmax(self.fc(X), dim=1)
        return y
