#!/usr/bin/env python
#_*_coding:utf-8_*_


import numpy as np
import pandas as pd

feature_index = pd.read_csv("featureindex.csv")

def select_features(allfeatures):
    print('Feature selection...')
    new_features = []
    orignal_data = pd.DataFrame(allfeatures)
    for i in list(feature_index):
        new_features.append(orignal_data[int(i)-1])
    features = np.array(new_features).T

    return features







