#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:30:57 2019

@author: julienh
"""
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.preprocessing import StandardScaler  

from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

from descriptors import create_is_train

def classification_all(descriptors,labels):
    
    N_sample, N_features = descriptors.shape
    is_train = create_is_train(N_sample)
    
    scaler = StandardScaler()
    scaler.fit(descriptors[is_train])    
    descriptors = scaler.transform(descriptors)

    
    clf = RandomForestClassifier(n_estimators=100, max_depth=2,
                                 random_state=0)
    clf.fit(descriptors[is_train], labels[is_train])
    predictions = clf.predict(descriptors[np.invert(is_train)])
    
    print(clf.feature_importances_)
    
    print(confusion_matrix(labels[np.invert(is_train)],predictions))  
    print(classification_report(labels[np.invert(is_train)],predictions))  
    print(accuracy_score(labels[np.invert(is_train)], predictions))  