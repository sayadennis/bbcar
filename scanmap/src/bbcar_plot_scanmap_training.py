import os
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt

dn = '/projects/b1122/saya/scanmap_data/train_record/genethres_clin_mut'

nc = 50
wcls = 50
C = 1.
# lr = 1.

# fn = 'scanmap_k%d_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
fn = 'scanmap_k%d_wcls%s_C%s.p' % (nc, wcls, C)

training = torch.load(os.path.join(dn, fn), map_location=torch.device('cpu'))

plt.plot(training['report']['loss'])
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Learning curve: k=%d, wcls=%s, C=%s' % (nc, wcls, C)) # , lr=%s // , lr
plt.savefig(os.path.join(dn, 'learning_curve_scanmap_k%d_wcls%s_C%s.png' % (nc, wcls, C))) # add lr 
