"""
Author: Sailee Sansgiri (sailee.s.sansgiri@jyu.fi)
"""

import pdb
import csv
import traceback
import sklearn
import sklearn.ensemble
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

X_COLUMNS = ['Maximum_ROM','Passive_ROM','Passive_PT','median_fast_DSRT_sol','median_fast_DSRT_MG', 'EMG_fast_normalized_sol', 'EMG_fast_normalized_mg'] #input hyperreflexia+mecha
Y_COLUMNS = ['pt_10'] #output hyperreflexia
Y_GROUP   = 'Group' #CP or TD
Y_COL1 = Y_COLUMNS[0]


def feat_importances(df, title):
    """
    Params
    ------
    df   : pd.Dataframe of a particular group
    title: str, For plotting purposes 
    """

    try:
        df_X = df[X_COLUMNS]
        df_Y = df[Y_COLUMNS]
        df_Y1=df[Y_COL1]


        X=df_X.to_numpy()   
        Y1=df_Y1.to_numpy()
    
        clf1=sklearn.ensemble.RandomForestRegressor(n_estimators=100, oob_score=True)
        clf1.fit(X,Y1)
        score1=clf1.score(X,Y1)
       

        columns = ['Group'] + ['Output'] + list(df_X.columns) + ['OOB Score']+['R2 Value']
        data    = [[title, Y_COL1] + list(clf1.feature_importances_) + [clf1.oob_score_] + [score1]]
                    
        return np.array(columns), np.array(data)
    
    except:
        traceback.print_exc()
        pdb.set_trace()

if __name__ == "__main__":

    try:
        df = pd.read_excel('ML_hyper-resistance_pt10.xlsx')
        df.dropna(inplace=True)

        df_CP = df[df[Y_GROUP] == 'CP']
        df_TD = df[df[Y_GROUP] == 'TD']

        columns, feat_imp_cp = feat_importances(df_CP, title='CP')
        columns, feat_imp_td = feat_importances(df_TD, title='TD')

        columns, feat_imp = feat_importances(df, title='CP+TD')
        
        df_feat_imp = pd.DataFrame(feat_imp, columns=columns)

        df_feat_imp2 = pd.DataFrame(np.vstack((feat_imp_cp, feat_imp_td)), columns=columns)
        df_feat_imp.to_csv('feature_importance_pt10.csv', index=False)
        df_feat_imp2.to_csv('feature_importance_pt10_CPTD.csv', index=False)
        pdb.set_trace()
    
    except:
        traceback.print_exc()
        pdb.set_trace()




pdb.set_trace()
