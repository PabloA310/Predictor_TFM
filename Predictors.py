# -*- coding: utf-8 -*-
"""
Title: Machine learning using protein features

@author: Pablo Aliaga Gaspar
"""

  
#%%
# Initial stuff

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import f1_score
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
                      "only introduce 'DF' without '.txt'")
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
    
    model_group = OptionGroup(parser, "Models")                     
    model_group.add_option("-k", "--knn",
                           action = "store_true",
                           dest = "knn",
                           default = False,
                           help = "Apply KNN model [default: %default]")
    model_group.add_option("-l", "--logreg",
                           action = "store_true",
                           dest = "logreg",
                           default = False,
                           help = "Apply Logistic regression model "
                           "[default: %default]")
    model_group.add_option("-t", "--tree",
                           action = "store_true",
                           dest = "tree",
                           default = False,
                           help = "Apply Decision tree model [default: %default]")
    model_group.add_option("-m", "--mlp",
                           action = "store_true",
                           dest = "mlp",
                           default = False,
                           help = "Apply Multilayer perceptrons model "
                           "[default: %default]")
    model_group.add_option("-g", "--gnb",
                           action = "store_true",
                           dest = "gnb",
                           default = False,
                           help = "Apply Gaussian Naive Bayes model "
                           "[default: %default]")
    parser.add_option_group(model_group)
    
    par_group = OptionGroup(parser, "Parameters")
    par_group.add_option("-n", "--neighbor",
                         action = "store",
                         type = "int",
                         dest = "best_n",
                         default = None,
                         help = "Int number to select n_neighbours "
                         "for KNN model [default: Best N for each data]")
    par_group.add_option("-c", "--clog",
                         action = "store",
                         type = "float",
                         dest = "c_log",
                         default = 1,
                         help = "Float number to select C for Logistic "
                         "regression model [default: %default]")
    par_group.add_option("--mdepth",
                         action = "store",
                         type = "int",
                         dest = "max_depth",
                         default = 8,
                         help = "Int number to select max depth for Decision "
                         "tree model [default: %default]")
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
    
def test_best_n(train_data, train_condition, test_data, test_condition):
    """Test the best n_neighbors (1-10) for the model KNN"""
    test_accuracy = []
    neighbors_settings = range(1, 10)
    for n_neighbors in neighbors_settings:
        knn = KNeighborsClassifier(n_neighbors = n_neighbors)
        knn.fit(train_data, train_condition)
        test_accuracy.append(knn.score(test_data, test_condition))
    best_n = test_accuracy.index(max(test_accuracy)) + 1
    return best_n
    
def crossval_model(model, total_data, total_condition):
    """Obtain the cross validation score for a model."""
    kfold = KFold(n_splits = 5, shuffle = True, random_state = 42)
    cv_score = cross_val_score(model, total_data, total_condition, cv = kfold)
    cv_score_mean = cv_score.mean()
    return(cv_score, cv_score_mean)

def make_precRec(test_data, test_condition, model_list, 
                 model_list_names):
    """Make a preccision/recall curve for each model and obtain average
    precision for them."""
    model_list_zip = zip(model_list, model_list_names)
    for model, model_name in model_list_zip:
        precision_model, recall_model, _ = precision_recall_curve(
        test_condition, model.predict_proba(test_data)[:, 1])
        ap_model = average_precision_score(
                test_condition, model.predict_proba(test_data)[:, 1])
        plt.plot(recall_model, precision_model, label =
                 '{}, AP = {:.2f}'.format(model_name, ap_model))  
    plt.title("Precision/Recall")
    plt.ylabel("Precision")
    plt.xlabel("Recall")
    plt.legend(loc = "lower left")
    plt.savefig('precrec_' + aux.name + '.png')
    plt.clf()   
    
def make_roc(test_data, test_condition, model_list, model_list_names):
    """Make a receiver operating characteristics (ROC) curve 
    and obtain the area under the curve (AUC) for each model."""
    model_list_zip = zip(model_list, model_list_names)
    for model, model_name in model_list_zip:
        fpr_model, tpr_model, _ = roc_curve(
                test_condition, model.predict_proba(test_data)[:, 1])
        model_auc = roc_auc_score(test_condition, model.predict_proba(
                test_data)[:, 1])
        plt.plot(fpr_model, tpr_model, label = '{}, AUC = {:.2f}'.format(
                model_name, model_auc))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate (Recall)")
    plt.title("ROC curve")
    plt.legend(loc = "best")
    plt.savefig('roc_' + aux.name + '.png')

def make_metdf(train_data, train_condition, test_data, test_condition,
               model_list, model_list_names,
               crossval_mean_list):
    """Make a data.frame of metrics for each model."""
    tr_score_list = []
    test_score_list = []
    f_score_list = []
    auc_list = []
    rec_list = []
    for model in model_list:
        tr_score_list.append(model.score(train_data, train_condition))
        test_score_list.append(model.score(test_data, test_condition))
        f_score_list.append(f1_score(test_condition, model.predict(test_data)))
        rec_list.append(recall_score(test_condition, model.predict(test_data)))
        auc_list.append(roc_auc_score(test_condition, 
                                      model.predict_proba(test_data)[:, 1]))
    dictionary = {"Training set score" : tr_score_list, 
              "Test set score" : test_score_list,
              "Cross validation" : crossval_mean_list,
              "F-1 score" : f_score_list,
              "Recall" : rec_list,
              "Area under the curve" : auc_list}
    metdf = pd.DataFrame.from_dict(dictionary, orient = "index", 
                           columns = model_list_names)
    return(metdf)

#%%
"""
Import and prepare data
"""
# Run optparse function
parser, options = main()

# Load data
aux = pd.DataFrame()
aux.name = options.filename

DF = pd.read_csv(aux.name + '.txt', sep = "\t")
DF_data, DF_condition = prepare_data(DF)

# Split the dataset for testing
X_train, X_test, y_train, y_test = train_test_split(
DF_data, DF_condition, random_state = 0)

# Scale data
X_train_scaled = scale_data(X_train, X_train)
X_test_scaled = scale_data(X_test, X_train)

#%%
"""
Different prediction models
"""
# List to be called later in functions
crossval_mean_list = []
model_list = []
model_list_names = []

# K nearest neighbors(KNN)
if options.knn is True:
    if options.best_n is None:
        n = test_best_n(X_train_scaled, y_train, X_test_scaled, y_test)
    else:
        n = options.best_n
    knn = KNeighborsClassifier(n_neighbors = n)
    knn.fit(X_train_scaled, y_train)

    # Cross validation
    crossval_knn, crossval_knn_mean = crossval_model(knn, DF_data,
                                                    DF_condition)
    crossval_mean_list.append(crossval_knn_mean)
    model_list.append(knn)
    model_list_names.append("KNN")
    
#%%
# Logistic regression
if options.logreg is True:
    logreg = LogisticRegression(C = options.c_log, solver = 'liblinear')
    logreg.fit(X_train_scaled, y_train)

    # Cross validation
    crossval_logreg, crossval_logreg_mean = crossval_model(logreg, DF_data,
                                                           DF_condition)
    crossval_mean_list.append(crossval_logreg_mean)
    model_list.append(logreg)
    model_list_names.append("Logistic Regression")

#%%
# Decision Tree Classifier
if options.tree is True:
    tree = DecisionTreeClassifier(max_depth = options.max_depth,
                                  random_state = 0)
    tree.fit(X_train_scaled, y_train)

    # Cross validation
    crossval_tree, crossval_tree_mean = crossval_model(tree, DF_data,
                                                       DF_condition)
    crossval_mean_list.append(crossval_tree_mean)
    model_list.append(tree)
    model_list_names.append("Decision Tree")
    
#%%
# Neural network - Multilayer perceptrons (MLP)
if options.mlp is True:
    mlp = MLPClassifier(alpha = options.alpha, max_iter = 1000, 
                        random_state = 0)
    mlp.fit(X_train_scaled, y_train)

    # Cross validation
    crossval_mlp, crossval_mlp_mean = crossval_model(mlp, DF_data, 
                                                     DF_condition)
    crossval_mean_list.append(crossval_mlp_mean)
    model_list.append(mlp)
    model_list_names.append("MLP")
    
#%%
# Gaussian Naive Bayes Classifier
if options.gnb is True:
    gnb = GaussianNB()
    gnb.fit(X_train_scaled, y_train)

    # Cross validation
    crossval_gnb, crossval_gnb_mean = crossval_model(gnb, DF_data,
                                                     DF_condition)
    crossval_mean_list.append(crossval_gnb_mean)
    model_list.append(gnb)
    model_list_names.append("GNB")
    
#%%
"""
Plots and metrics for model evaluation
"""
# Data.frame with metrics for each model
if options.table_name is not False:
    metricDF = make_metdf(X_train_scaled, y_train, X_test_scaled, y_test, 
                      model_list, model_list_names, crossval_mean_list)
    with open('metric_' + aux.name + '.tsv', "w") as f:
        f.write(tabulate(metricDF, headers= metricDF.columns.values,
                         tablefmt="tsv"))

#%%
# Precision/recall curve and average precision
if options.precrec_name is not False:
    make_precRec(X_test_scaled, y_test, model_list, model_list_names)

#%%
# Receiver operating characteristics (ROC) curve and area under the curve (AUC)
if options.roc_name is not False:
    make_roc(X_test_scaled, y_test, model_list, model_list_names)

