#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import itertools
import warnings
warnings.filterwarnings('ignore')
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, r2_score
from sklearn.metrics import f1_score
from sklearn import metrics

TbData = pd.read_csv("../intermediates/Tb/result_aln_pf_Tb_tm_labeled.tsv", sep="\t")
TbData["FracOfPf"] = (TbData["qend"] -TbData["qstart"] + 1)/TbData["qlen"]
TbData["p_evalue"] = -np.log(TbData["evalue"] + 1e-300)
TbData["MismatchPerc"] = TbData["mismatch"]/TbData["alnlen"] *100
TbData["GapOpenPerc"] = TbData["gapopen"]/TbData["alnlen"] *100

ColumnNames = ['fident', 'alnlen', "MismatchPerc", 'GapOpenPerc', 'qlen', 'bits',
       'alntmscore', "FracOfPf","p_evalue", 'query', 'target','qstart',
       'qend', 'tstart', 'tend', 'tlen','PredPF', 'PF', 'PFstart', 'PFend', 'evalue', 'Status']
TbData = TbData[ColumnNames]
LabeledData = TbData[TbData["Status"]!=-1]
LabeledData = LabeledData.sample(frac=1, random_state=1).reset_index(drop=True)

X = LabeledData.values
y = LabeledData['Status'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size = 0.9, test_size=0.1, shuffle=False)
X_cv, X_test, y_cv, y_test = train_test_split(X_test, y_test, train_size = 0.5, test_size=0.5, shuffle=False)

Params = {"max_depth":[ None, 5, 7, 9],
"criterion":["gini", "entropy", "log_loss"],
"n_estimators" :[20, 100, 500, 1000]}

reportFile = open("../intermediates/ReportOfPerformanceCV.txt",'w')
reportFile.write("Model\tAccuracy\tF1-score\tAUC\n")
for dep, cri, n_es in itertools.product(Params["max_depth"], Params["criterion"],Params["n_estimators"]):
    rf = RandomForestClassifier(n_estimators=n_es, max_depth = dep, criterion= cri, n_jobs=-1, random_state=1)
    rf.fit(X_train[:,:8], y_train)
    CV_pred = rf.predict(X=X_cv[:,:8])
    CV_pred_prob = rf.predict_proba(X_cv[:,:8])[:, 1]
    acc = rf.score(X_cv[:,:8], y_cv)
    f1  = f1_score(y_cv, CV_pred)
    #fpr, tpr, _ = metrics.roc_curve(y_cv, CV_pred_prob)
    auc = metrics.roc_auc_score(y_cv, CV_pred_prob)

    reportFile.write("{}_estimators_{}_max_depth_{}_split_criterion\t".format(n_es, dep, cri,acc ))
    reportFile.write("{}\t{}\t{}\n".format(acc, f1, auc))
reportFile.close()
