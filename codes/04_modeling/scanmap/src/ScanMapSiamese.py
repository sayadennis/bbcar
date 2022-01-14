import torch
from torch import nn
from collections import  defaultdict
import numpy as np
import utils
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score

## this implements ScanMapSiamese transductive learning with demographics information and orthogonality constraints
## assumes X is passed as [Xtrain; Xtest], and only y_train is passed
class ScanMapSiamese(nn.Module):
    def __init__(self, X, cf = None, cfval = None, y = None, yval=None,
                 k = 10, n_iter = 10, eps = 1e-7,
                 floss = 'l2', weight_decay = 1e-5, wortho=1, wcls=0.1, C=1,
                 lr = 1e-2, verbose = False, seed=None, fn=None,
                 device=torch.device('cpu')):
        super(ScanMapSiamese, self).__init__()
        if seed is not None:
            torch.manual_seed(seed)
            np.random.seed(seed)
        self.n_iter = n_iter
        self.k = k
        self.floss = floss
        self.weight_decay = weight_decay
        self.lr = lr
        self.verbose = verbose
        self.eps = eps
        self.wcls = wcls # weight of classification loss
        self.wortho = wortho
        self.fn = fn
        self.C = C
        self.decomposed = False
        self.dcf = cf.shape[1]
        self.device = device
        self.celoss = -1
        self.report = defaultdict(list)
        self.__initfact__(X, cf, cfval, y, yval)

    def __initfact__(self, X, cf, cfval, y, yval):
        # self.n,self.m = X.shape
        # self.ntr = len(y)
        # self.nval = len(yval)
        self.n = torch.tensor(X.shape, dtype=int)[0].to(self.device)
        self.m = torch.tensor(X.shape, dtype=int)[1].to(self.device)
        self.ntr = torch.tensor(len(y), dtype=int).to(self.device)
        self.nval = torch.tensor(len(yval), dtype=int).to(self.device)
        self.nte = self.n - self.ntr - self.nval
        self.X = torch.from_numpy(X).float().to(self.device)
        self.cf = torch.from_numpy(cf).float().to(self.device)
        self.cfval = torch.from_numpy(cfval).float().to(self.device)
        
        self.scale = torch.tensor(torch.mean(self.X) / self.k).to(self.device)
        Wtr = torch.abs(torch.rand([self.ntr,self.k]).to(self.device) * self.scale).to(self.device)
        Wval = torch.abs(torch.rand([self.nval,self.k]).to(self.device) * self.scale).to(self.device)
        Wte = torch.abs(torch.rand([self.nte,self.k]).to(self.device) * self.scale).to(self.device)
        H = torch.abs(torch.rand([self.k,self.m])).to(self.device)
        self.Wtr = torch.nn.Parameter(Wtr)
        self.Wval = torch.nn.Parameter(Wval)
        self.Wte = torch.nn.Parameter(Wte)
        self.H = torch.nn.Parameter(H)
        self.identity = torch.eye(self.k, device=self.device)
        
        if self.floss == 'l2':
            self.loss_fac = utils.l2
        elif self.floss == 'kl':
            self.loss_fac = utils.kl_div

        self.y = torch.from_numpy(y).long().to(self.device)
        self.y_pairs = torch.cat([torch.repeat_interleave(self.y, self.ntr, dim=0).reshape(1,-1), self.y.repeat(1, self.ntr)], dim=0)
        self.fourclass = []
        for i in range(self.y_pairs.shape[1]):
            if torch.all(torch.eq(self.y_pairs[:,i],torch.tensor([0,0]).to(self.device))):
                self.fourclass.append(0)
            elif torch.all(torch.eq(self.y_pairs[:,i],torch.tensor([1,0]).to(self.device))):
                self.fourclass.append(1)
            elif torch.all(torch.eq(self.y_pairs[:,i],torch.tensor([0,1]).to(self.device))):
                self.fourclass.append(2)
            elif torch.all(torch.eq(self.y_pairs[:,i],torch.tensor([1,1]).to(self.device))):
                self.fourclass.append(3)

        w = 1 / pd.Series(self.fourclass).value_counts(normalize=True).sort_index().to_numpy() # initially y instead of fourclass
        self.loss_cls = nn.CrossEntropyLoss(weight=torch.from_numpy(w).float()).to(self.device) # weight=torch.from_numpy(w).float()
        ## EDIT FOR SIAMESE: change the dimensions of nn.Linear 
        self.fc = nn.Linear(self.k*2 + self.dcf*2, 4).to(self.device) # 4 = len(torch.unique(self.y))**2
            
        self.yval = torch.from_numpy(yval).long().to(self.device)
        self.opt = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

    def to(self,device):
        self.device = device
        self.X = self.X.to(device)
        self.cf = self.cf.to(device)
        self.cfval = self.cfval.to(device)
        self.y = self.y.to(device)
        self.yval = self.yval.to(device)        
        self.y_pairs = self.y_pairs.to(device)
        self.fourclass = self.fourclass.to(device)
        return super(ScanMapSiamese, self).to(device)

    def plus(self,X):
        X[X < 0] = 0 # self.eps
        return X

    def __autograd__(self,epoch):
        """
           autograd update, with gradient projection
        """
        self.opt.zero_grad()
        A = torch.cat([self.Wtr, self.Wval, self.Wte])
        # factorization loss
        l = self.loss_fac(A @ self.H, self.X) * self.k * self.n #   01/24/20
        # add l2 regularization for orthogonality
        AtA = torch.mm(torch.t(A), A)
        if torch.mean(AtA) > self.eps:
            AtA = AtA/torch.mean(AtA)
        l += self.loss_fac(AtA/self.k, self.identity) * self.wortho * self.k # scale invariant orthogonal
        # add classification loss 
        if self.y is not None:
            # EDIT FOR SIAMESE: create pairs of Wtr+cf, create new four-class y for these pairs 
            # Wtr_pairs = torch.cat([torch.repeat_interleave(self.Wtr, self.ntr, dim=0), self.Wtr.repeat(self.ntr, 1)], dim=1)
            # cf_pairs = torch.cat([torch.repeat_interleave(self.cf, self.ntr, dim=0), self.cf.repeat(self.ntr, 1)], dim=1)
            # y_pairs_onehot = nn.functional.one_hot(torch.tensor(fourclass)).to(self.device)
            ## Test code: sample from the pairs to reduce computation time
            # sample_ix = torch.randint(low=0, high=y_pairs.shape[0], size=(self.k*20,))
            self.celoss = self.loss_cls(
                self.fc(
                    torch.cat([
                        torch.cat([torch.repeat_interleave(self.Wtr, self.ntr, dim=0), self.Wtr.repeat(self.ntr, 1)], dim=1), 
                        torch.cat([torch.repeat_interleave(self.cf, self.ntr, dim=0), self.cf.repeat(self.ntr, 1)], dim=1)
                        ], dim=1) # [sample_ix,:]
                ), 
                torch.tensor(self.fourclass).to(self.device) # [sample_ix]
            ) # y_pairs_onehot
            l = l + self.celoss * self.wcls * self.n # 01/24/20
            # print('cross entropy: %.4f' % (self.loss_cls(self.fc(self.Wtr), self.y)))
            for p in self.fc.parameters():
                # has two parameters (cls #, ft #) and (cls #), they are the weight and biases, /self.k to add normalized weight regularization
                l = l + p.pow(2).sum() * self.C 
                # print('complexity: %.4f' % (p.pow(2).sum()))

        l.backward()
        self.opt.step()
        ## grad projection
        self.Wtr.data = self.plus(self.Wtr.data)
        self.Wval.data = self.plus(self.Wval.data)
        self.Wte.data = self.plus(self.Wte.data)
        self.H.data = self.plus(self.H.data)
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
                y_val_pred = self.__predict__(self.Wval, self.cfval)
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

    def fit_transform(self):
        if not self.decomposed:
            self.fit()
            # detach all params including the linear fc layer
            for p in self.parameters():
                p.requires_grad = False
        if self.device.type == 'cuda':
            return [self.Wtr.detach().cpu().numpy(), self.Wval.detach().cpu().numpy(), self.Wte.detach().cpu().numpy(), self.H.detach().cpu().numpy()]
        else:
            return [self.Wtr.detach().numpy(), self.Wval.detach().numpy(), self.Wte.detach().numpy(), self.H.detach().numpy()]

    def predict(self, w, cf):
        w = torch.from_numpy(w).float()
        cf = torch.from_numpy(cf).float()                        
        if self.device.type == 'cuda':
            w = w.to(self.device)
            cf = cf.to(self.device)
        y = self.__predict__(w, cf)
        if self.device.type == 'cuda':
            return y.detach().cpu().numpy()
        else:
            return y.detach().numpy()

    def __predict__(self, w, cf):
        # EDIT FOR SIAMESE: create pairs of w+cf
        index_map = torch.cat([torch.repeat_interleave(torch.arange(w.shape[0]), w.shape[0], dim=0).reshape(1,-1), torch.arange(w.shape[0]).repeat(1, w.shape[0])], dim=0)
        w_pairs = torch.cat([torch.repeat_interleave(w, w.shape[0], dim=0), w.repeat(w.shape[0], 1)], dim=1)
        cf_pairs = torch.cat([torch.repeat_interleave(cf, cf.shape[0], dim=0), cf.repeat(cf.shape[0], 1)], dim=1)
        # y_pairs_onehot = torch.argmax(self.fc(torch.cat([w_pairs, cf_pairs],1)), dim=1)
        y_pairs_fourclass = torch.argmax(self.fc(torch.cat([w_pairs, cf_pairs],1)), dim=1)
        # EDIT FOR SIAMESE: add lines to convert predictions on pairs to preds and probs on individuals
        y_pairs = torch.tensor([], dtype=int)
        for i in range(len(y_pairs_fourclass)):
            if torch.eq(y_pairs_fourclass[i],torch.tensor([0]).to(self.device)):
                y_pairs = torch.cat([y_pairs, torch.tensor([0,0])], dim=0)
            elif torch.eq(y_pairs_fourclass[i],torch.tensor([1]).to(self.device)):
                y_pairs = torch.cat([y_pairs, torch.tensor([1,0])], dim=0)
            elif torch.eq(y_pairs_fourclass[i],torch.tensor([2]).to(self.device)):
                y_pairs = torch.cat([y_pairs, torch.tensor([0,1])], dim=0)
            elif torch.eq(y_pairs_fourclass[i],torch.tensor([3]).to(self.device)):
                y_pairs = torch.cat([y_pairs, torch.tensor([1,1])], dim=0)
        y_pairs = y_pairs.reshape(-1,2)
        y_pairs = torch.t(y_pairs)
        y = torch.tensor([], dtype=int)
        for i in range(w.shape[0]):
            mapping = torch.where(index_map==i)
            votes = []
            for j,k in zip(mapping[0], mapping[1]):
                votes.append(y_pairs[j,k])
            if torch.sum(torch.tensor(votes))/len(votes) >= 0.5:
                y = torch.cat([y, torch.tensor([1])], dim=0)
            else:
                y = torch.cat([y, torch.tensor([0])], dim=0)
        return y
