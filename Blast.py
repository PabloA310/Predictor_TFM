#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Blast curves

@author: Pablo Aliaga Gaspar
"""


#%%
# Initial stuff

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score


#%%
DF = pd.read_csv("eValues_ltp.txt", sep = "\t")

DF['p-value'] = 1 - np.exp(-DF['e-value'])
DF['probability'] = 1 - DF['p-value']
DF['x'] = 1 - DF['e-value']


#%%
precision_model, recall_model, _ = precision_recall_curve(DF['condBin'],
    DF['probability'])
ap_model = average_precision_score(DF['condBin'], DF['probability'])
plt.plot(recall_model, precision_model, label = 
         'Blast, AP = {:.2f}'.format(ap_model))  
plt.title("Precision/Recall")
plt.ylabel("Precision")
plt.xlabel("Recall")
plt.legend(loc = "lower left")
plt.savefig("precrec_Blast_ltp.png")
plt.clf()

fpr_model, tpr_model, _ = roc_curve(DF['condBin'], DF['probability'])
model_auc = roc_auc_score(DF['condBin'], DF['probability'])
plt.plot(fpr_model, tpr_model, label = 'Blast, AUC = {:.2f}'.format(model_auc))
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate (Recall)")
plt.title("ROC curve")
plt.legend(loc = "best")
plt.savefig('roc_Blast_ltp.png')


  