#!/usr/bin/env python
# _*_coding:utf-8_*_

from flask import Flask,render_template,request
import feature_extraction
import read_fasta_sequences
import numpy as np
import pandas as pd
import joblib
import feature_selection
import sys,re
from lightgbm import LGBMClassifier
from check_sequences import *

# initialize flask application
app = Flask(__name__)

@app.route('/')
def main():
   return render_template('prediction.html')

@app.route('/submit',methods=("GET", "POST"))
def predict():
    if request.method == "POST":
        sequencedata = request.form["sequence"]

        fastas, sequence_name = read_fasta_sequences.read_protein_sequences(sequencedata)
        check=check_sequences(sequencedata, fastas)
        if check == False:

            return render_template('predicting.html')
        else:
            encodings = feature_extraction.get_features(fastas)
            #print(encodings)
            #np.savetxt("feature_extraction.csv",encodings, fmt='%5s', delimiter=',')
            newdataset = feature_selection.select_features(encodings)
            #print(newdataset)
            #np.savetxt("feature_selection.csv",newdataset, fmt='%5s', delimiter=',')
            X = newdataset
            scale = joblib.load(r"./models/scaler.pkl")
            X = scale.transform(X)
            model = joblib.load(r"./models/LGBM_model.pkl")
            y_pred_prob = model.predict_proba(X)
            #np.savetxt('LGBM_hybrid_train.csv', y_pred_prob, fmt='%5s', delimiter=',')
            df_out = pd.DataFrame(np.zeros((y_pred_prob.shape[0], 3)),
                                  columns=["Sequence name", "Predicted epitope", "Probability"])
            print(df_out)
            y_pred = model.predict(X)
            for i in range(y_pred.shape[0]):
                df_out.iloc[i, 0] = str(sequence_name[i])
                if y_pred[i] == 1:
                    df_out.iloc[i, 1] = "yes"
                    df_out.iloc[i, 2] = "%.2f%%" % (y_pred_prob[i, 1] * 100)
                if y_pred[i] == 0:
                    df_out.iloc[i, 1] = "no"
                    df_out.iloc[i, 2] = "%.2f%%" % (y_pred_prob[i, 0] * 100)

            df_out.to_excel("results.xls")
            results_all = pd.read_excel('results.xls', header=0, index_col=0).to_numpy()
            return render_template('result.html', Fasta=sequencedata, result=results_all, Predicting="Output results ")
            return render_template('home.html')

if __name__=="__main__":
    app.run(debug=True, port=7005)


