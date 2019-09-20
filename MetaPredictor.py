#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Machine learning using a meta-predictor (K-mers/Protein features)

@author: Pablo Aliaga Gaspar
"""

#%%
# Initial stuff

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve
from sklearn.metrics import recall_score
from tabulate import tabulate
from optparse import OptionParser
from optparse import OptionGroup

def main():
    """Optparse function for running the script in BASH"""
    usage = "usage: %prog [options] arg1, arg2, arg3"
    parser = OptionParser(usage)
    parser.add_option("-f", "--file", 
                      action = "store",
                      dest = "filename", 
                      help = "Name of 'data.txt' file with TAB as " 
                      "separator and with format: [First column: Condition], " 
                      "[Second column: Binary condition], [Thir column: Label],"
                      " [Rest: Features]. Use example: filename = 'DF.txt', "
                      "only introduce DF without ''")
    parser.add_option("-s", "--seq",
                      action = "store",
                      dest = "seqname",
                      help = "Give 'seq_data.txt' file with TAB as separator "
                      "and with format: [First column: Binary condition], "
                      "[Second column: Sequence].")
    parser.add_option("-d", "--dataframe",
                      action = "store_true",
                      dest = "table_name",
                      default = False,
                      help = "Create a metric data.frame [default: %default]")
    parser.add_option("-p", "--precrec",
                      action = "store_true",
                      dest = "precrec_name",
                      default = False,
                      help = "Create a precision/recall curve "
                      "[default: %default]")
    parser.add_option("-r", "--roc",
                      action = "store_true",
                      dest = "roc_name",
                      default = False,
                      help = "Create a ROC curve [default: %default]")
    
    par_group = OptionGroup(parser, "Parameters")
    par_group.add_option("-n", "--neighbor",
                         action = "store",
                         type = "int",
                         dest = "k_n",
                         default = 3,
                         help = "Int number to select n_neighbours "
                         "for KNN model [default: %default]")
    par_group.add_option("-c", "--csvc",
                         action = "store",
                         type = "float",
                         dest = "c_svc",
                         default = 1,
                         help = "Float number to select C for SVC "
                         "model [default: %default]")
    par_group.add_option("-a", "--alpha",
                         action = "store",
                         type = "float",
                         dest = "alpha",
                         default = 0.0001,
                         help = "Float number to select alpha for Multilayer "
                         "perceptrons model [default: %default]")
    parser.add_option_group(par_group)
    
    (options, args) = parser.parse_args()
    return parser, options  

def get_models(c_svc, k_n, alpha):
    """Generate a library of base learners."""
    svc = SVC(C = c_svc, gamma = 'auto', kernel = 'linear',
              probability = True)
    knn = KNeighborsClassifier(n_neighbors = k_n)
    mlp = MLPClassifier(alpha = alpha, max_iter = 1000, random_state = 42)
    models = {'svm': svc,
              'knn': knn,
              'mlp-nn': mlp
              }
    return models

def getKmers(sequence, k):
    """Get K-mers from a sequence, this function is used as auxiliar in
    --create_text function."""
    return [sequence[x:x + k].lower() for x in range(
            len(sequence) - k + 1)]

def create_text(seqFile, k):
    """Prepare data for --CountVectorizer() from table of sequences"""
    prot = pd.read_csv(seqFile, sep = "\t")
    prot['words'] = prot.apply(lambda x: getKmers(x['sequence'], k),
        axis = 1)
    prot = prot.drop('sequence', axis = 1)
    prot_texts = list(prot['words'])
    for item in range(len(prot_texts)):
        prot_texts[item] = ' '.join(prot_texts[item])
    y = prot.iloc[:,0].values
    return prot_texts, y

def prepare_data(DF):
    """Prepare data and condition from data.frame."""
    DF_data = DF.iloc[:, 3 : len(DF.columns)]
    DF_condition = DF.iloc[:,1]
    return(DF_data, DF_condition)
    
def scale_data(total_data, train_data):
    """Scale data with mean and standard deviation."""
    mean_on_train = train_data.mean(axis = 0)
    std_on_train = train_data.std(axis = 0)
    data_scaled = (total_data - mean_on_train) / std_on_train
    return(data_scaled)
    
def train_base_learners(base_learners, inp1, inp2, out1, out2):
    """Train all base learners in the library."""
    for i, (name, m) in enumerate(base_learners.items()):
        if name != 'svm':
            m.fit(inp2, out2)
        else:
            m.fit(inp1, out1)

def predict_base_learners(pred_base_learners, inp1, inp2):
    """Generate a prediction matrix."""
    P = np.zeros((inp1.shape[0], len(pred_base_learners)))

    for i, (name, m) in enumerate(pred_base_learners.items()):
        if name != 'svm':
            p = m.predict_proba(inp2)
            P[:, i] = p[:, 1]
        else:
            p = m.decision_function(inp1)
            P[:, i] = p
    return P

def ensemble_predict(base_learners, meta_learner, inp1, inp2):
    """Generate predictions from the ensemble."""
    P_pred = predict_base_learners(base_learners, inp1, inp2)
    return P_pred, meta_learner.predict_proba(P_pred)[:, 1]
    
def make_precRec(predict_data, test_condition):
    """Make a preccision/recall curve for the ensemble and obtain average
    precision for it."""
    precision, recall, _ = precision_recall_curve(test_condition, predict_data)
    ap = average_precision_score(test_condition, predict_data)
    plt.plot(recall, precision, label = 'Ensemble, AP = {:.2f}'.format(ap))
    plt.title("Precision/Recall")
    plt.ylabel("Precision")
    plt.xlabel("Recall")
    plt.legend(loc = "lower left")
    plt.savefig('meta_precrec_' + aux.name + '.png')
    plt.clf()   
    
def make_roc(predict_data, test_condition):
    """Make a receiver operating characteristics (ROC) curve 
    and obtain the area under the curve (AUC) for the ensemble."""
    fpr, tpr, _ = roc_curve(test_condition, predict_data)
    auc = roc_auc_score(test_condition, predict_data)
    plt.plot(fpr, tpr, label = 'Ensemble, AUC = {:.2f}'.format(auc))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate (Recall)")
    plt.title("ROC curve")
    plt.legend(loc = "best")
    plt.savefig('meta_roc_' + aux.name + '.png')

def make_metdf(train_data, train_condition, test_data,
               test_condition, predict_data):
    """Make a data.frame of metrics for the ensemble."""
    tr_score = meta_learner.score(train_data, train_condition)
    test_score = meta_learner.score(test_data, test_condition)
    f_score = f1_score(test_condition, meta_learner.predict(test_data))
    auc = roc_auc_score(test_condition, predict_data)
    rec = recall_score(test_condition, meta_learner.predict(test_data))
    dictionary = {"Training set score" : tr_score, 
              "Test set score" : test_score,
              "F-1 score" : f_score,
              "Recall" : rec,
              "Area under the curve" : auc}
    metdf = pd.DataFrame.from_dict(dictionary, orient = "index", 
                                   columns = ['Ensemble'])
    return(metdf)
        

#%%
"""
Import and prepare data
"""
# Run optparse function
parser, options = main()

# Features data
# Load data
aux = pd.DataFrame()
aux.name = options.filename

# Import and prepare data
DF = pd.read_csv(aux.name + ".txt", sep = "\t")
DF_data, DF_condition = prepare_data(DF)

# Split the dataset for testing
X_train2, X_test2, y_train2, y_test2 = train_test_split(
DF_data, DF_condition, random_state = 42)

# Scale data
X_train2_scaled = scale_data(X_train2, X_train2)
X_test2_scaled = scale_data(X_test2, X_train2)

# Split again for training base models and meta-learner
X_train2_base, X_pred2_base, y_train2_base, y_pred2_base = train_test_split(
    X_train2_scaled, y_train2, test_size=0.5, random_state=42)

#%%
# K-mers data
# Import and prepare data
prot_texts, y = create_text(options.seqname, k = 5) 

# Count K-mers and split data
cv = CountVectorizer()
X = cv.fit_transform(prot_texts)       
X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state = 42)

# Split again for training base models and meta-learner
X_train_base, X_pred_base, y_train_base, y_pred_base = train_test_split(
    X_train, y_train, test_size=0.5, random_state=42)

#%%
"""
Prepare and train models
"""
# Base learners
base_learners = get_models(options.c_svc, options.k_n, options.alpha)

#%%
# Define meta learner
meta_learner = RandomForestClassifier(n_estimators=100, 
                                      max_depth=3, random_state=42)

#%%
# Train base learners
train_base_learners(base_learners, X_train_base, X_train2_base,
                    y_train_base, y_train2_base)
P_base = predict_base_learners(base_learners, X_pred_base, X_pred2_base)

#%%.
# Train meta learner
meta_learner.fit(P_base, y_pred_base)

P_pred, p = ensemble_predict(base_learners, meta_learner, X_test, X_test2)

#%%
"""
Plots and metrics for model evaluation
"""
# Data.frame with metrics for each model
if options.table_name is not False:
    metricDF = make_metdf(P_base, y_pred_base, P_pred, y_test, p)
    with open('meta_metric_' + aux.name + '.tsv', "w") as f:
        f.write(tabulate(metricDF, headers= metricDF.columns.values,
                         tablefmt="tsv"))

#%%
# Precision/recall curve and average precision
if options.precrec_name is not False:
    make_precRec(p, y_test)

#%%
# Receiver operating characteristics (ROC) curve and area under the curve (AUC)
if options.roc_name is not False:
    make_roc(p, y_test)
