#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:35:53 2020

@author: shuvomsadhuka
"""
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
import os
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb


def distances_pairwise(df):
    muts = df[['x', 'y', 'z']].to_numpy()
    return(euclidean_distances(muts, muts))


def covariance_matrix(dists, t):
    #cov = 1/(dists+1)
    cov = np.exp(-np.square(dists)/(2*t**2))
    return(cov)

def select_mutations(df, num_muts):
    return(df.sample(num_muts))


def set_pi_and_tau(var_to_explain, pi_sum_ratio, comb_maf, eps=1):
    # see math - solve system of linear eqs with (x+y)/(x+y+1) = var_to_explain
    # and x/(x+y) = pi_sum_ratio, then rescale pi
    sum_pi_tau = var_to_explain/(1 - var_to_explain)
    unscaled_pi = sum_pi_tau * pi_sum_ratio
    scaled_pi = np.sqrt(unscaled_pi / (2 * (1 - comb_maf) * comb_maf))
    tau_sq = sum_pi_tau - unscaled_pi
    return(scaled_pi, tau_sq)


def run_logistic(X, y):
    logitr = LogisticRegression(penalty='none', fit_intercept=False)
    logitr.fit(X.T, y)
    return(logitr)


def run_linear(X, y):
    linear = LinearRegression( fit_intercept=False)
    linear.fit(X.T, y)
    return(linear)


def make_protein_df(protein_id):
    protein_dir = protein_id[0:2] + '/' + protein_id[2:4] + '/' + protein_id[4:6]
    file_repo = '../Scale_Up/SWISS-MODEL_Repository/' + protein_dir + '/swissmodel/'
    
    if os.path.isdir(file_repo):
        pdb_file = file_repo + str(os.listdir(file_repo)[0])
        
        ppdb = PandasPdb().read_pdb(pdb_file)
        protein_structure = pd.DataFrame(ppdb.df['ATOM'])
        
        for i in ['x', 'y', 'z']:
            coord = i + '_coord'
            pos_df = protein_structure.groupby(['residue_number', 'chain_id'], as_index=False)[coord].mean()
            protein_structure = pd.merge(protein_structure, pos_df, on=['residue_number', 'chain_id'])
        
        protein_structure = protein_structure.rename(columns={'x_coord_y': 'x', 'y_coord_y': 'y', 'z_coord_y': 'z'})
        protein_structure = protein_structure.drop_duplicates(subset=['x', 'y', 'z'], keep="first") 
        return(protein_structure)
    
    else:
        raise SyntaxError('Protein directory specified does not exist!')
    