#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 18:03:24 2020

@author: shuvomsadhuka
"""

from biopandas.pdb import PandasPdb
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
from scipy.stats import mannwhitneyu
from scipy import stats
import matplotlib.pyplot as plt
import sys


def make_dataframe(df, data, sequence):
    """Concatenates two files (a protein index and residue-mutation score file) into a single Pandas dataframe object.
    INPUTS
    =======
    df: Pandas DataFrame
      A DataFrame giving the residue numbers, locations, and amino acid information for each residue in a specific protein
    data: Pandas DataFrame
      A DataFrame storing all the mutations with residue numbers
    sequence: array
        an array storing the 3-to-1 letter encodings of amino acids
    RETURNS
    ========
    df_tot: a concatenated DataFrame of all mutations and corresponding locations, residue numbers, etc.
    """

    # extract residue number and name of each mutation from column
    data['residue_number'] = data[0].str.strip().str[3:].str[:-1]
    data = data[data['residue_number'].astype(int).isin(list(set(df['residue_number'])))]
    data['residue_name_convert'] = data[0].str.strip().str[2:3]
    
    #print(data)
    
    # get x, y, z coords of each residue as mean coordinate of all atoms making up residue
    for i in ['x', 'y', 'z']:
        coord = i + '_coord'
        pos_df = df.groupby('residue_number', as_index=False)[coord].mean()
        df = df.join(pos_df.set_index('residue_number'), on='residue_number', how='right', lsuffix='_left', rsuffix='_right')
    
    pd.set_option('display.max_columns',None)
    #print(df)
    
    # drop all columns except residue name, number, and coordinates
    df = df[['residue_name', 'residue_number', 'x_coord_right', 'y_coord_right', 'z_coord_right']]
    df = df.drop_duplicates(subset='residue_number', keep="first")
    
    # 3-to-1 amino acid string conversion
    #df['residue_name_convert'] = list(sequence['residue_name'])
    #print(data)
    # type casting
    data = data.rename(columns={1: "score", 2: "p-val"})
    data = data.astype({'residue_number': 'int64'})
    
    # merge datasets
    df_tot = pd.merge(df, data, on="residue_number")
    
    df_tot['risk_prob'] = df_tot['residue_number']
    
    df_tot.loc[df_tot['score'] < 0, 'risk_prob'] = df_tot['p-val']/2
    df_tot.loc[df_tot['score'] > 0, 'risk_prob'] = 1-df_tot['p-val']/2
    
    #print("HELLO")
    #print(df_tot)
    #print(df_tot.count)
    return df_tot


def get_prob_vec(df, risk):
    """Returns a vector of the average distance of each mutation of a specific type (risk or protective) to all other mutations
    INPUTS
    =======
    df: Pandas DataFrame
      A concatenated DataFrame of all mutations and corresponding locations, residue numbers, etc.
    risk: boolean
      A boolean indicating whether mutation group is risk (True) or protective (False)
    RETURNS
    ========
    vec: a vector of the average distance of each mutation of a specific type (risk or protective) to all other mutations
    """
    # subset dataframe based on risk or non-risk
    if risk:
        return df['risk_prob'].to_numpy()
    else:
        return(np.ones(len(df['risk_prob'])-df['risk_prob'].to_numpy))
        


def prob_mean_dist(prots, prob_vec):
    """Returns the mean distance over all mutations (of a specific type), for permutation test
    INPUTS
    =======
    prots: array
      An nx3 array of coordinates of mutations of a specific type
    RETURNS
    ========
    a float of the mean distance in 3D space of all mutations to each other
    """
    # get pairwise distances of all rows
    dists_prots = euclidean_distances(prots, prots)
    #print(dists_prots)
    prob_mat = np.outer(prob_vec, prob_vec)
 
    dists_mat = np.multiply(dists_prots, prob_mat)
    
    m = dists_mat.shape[0]
    strided = np.lib.stride_tricks.as_strided
    s0,s1 = dists_mat.strides
    
    Y = strided(dists_mat.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)
    X = strided(prob_mat.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)

    return(np.mean(Y)/np.mean(X))


def run_prob_permutation(muts_df, df, prob_vec, num_runs, trials=0):
    #print(muts_df)
    x = muts_df.iloc[:,[2,3,4]].to_numpy()
    emp_risk = prob_mean_dist(x, prob_vec)
    emp_prot = prob_mean_dist(x, np.ones(len(prob_vec)) - prob_vec)
    
    for_perms = prob_vec
    
    dist_risks = []
    dist_prots = []
    
    for i in range(num_runs):
        np.random.shuffle(for_perms)
        
        dist_risks.append(prob_mean_dist(x, for_perms))
        dist_prots.append(prob_mean_dist(x, np.ones(len(for_perms)) - for_perms))
            
    ret1 = stats.percentileofscore(dist_risks, emp_risk)
    ret2 = stats.percentileofscore(dist_prots, emp_prot)
    
    if (ret1 < 5 or ret2 < 5) and trials==0:
        ret = run_prob_permutation(muts_df, df, prob_vec, 50000, trials=1)
        ret1 = ret[0]
        ret2 = ret[1]
    return(ret1, ret2)
    
