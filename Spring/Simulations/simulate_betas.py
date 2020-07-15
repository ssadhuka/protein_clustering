#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:58:53 2020

@author: shuvomsadhuka
"""
import os
from biopandas.pdb import PandasPdb
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
#from scipy.stats import multivariate_normal
import scipy
import numpy as np

def genotype_matrix(cov, num_subjects):
    maf = np.random.uniform(0, 0.5, size=cov.shape[0])
    
    genotypes = np.array(list(map(lambda freq: np.random.binomial(2, freq, size=num_subjects), maf)))
    #print(genotypes.shape)
    return(genotypes)
   


def make_betas(cov):
    betas = np.random.multivariate_normal(np.zeros(cov.shape[0]), cov)
    #print(betas)
    return(betas)


def phenotype_vector(genotypes, betas):
    prob_disease = scipy.special.expit(np.dot(genotypes.T, betas))


def covariance_matrix(dists):
    cov = 1/(dists+1)
    return(cov)


def select_mutations(df, num_muts):
    return(df.sample(num_muts))


def distances_pairwise(df):
    muts = df[['x', 'y', 'z']].to_numpy()
    return(euclidean_distances(muts, muts))


def fit_regression(genotypes, cov):
    second_order = genotypes.T.dot(genotypes)
    print(second_order.shape)
    #print()
    pass
    

def likelihood_ratio_test():
    pass
    

def main():
    protein_id = 'Q01113'
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
        
        #protein_structure = select_mutations(protein_structure, 10)
        dists = distances_pairwise(protein_structure)
        cov = covariance_matrix(dists)
        betas = make_betas(cov)
        
        genotypes = genotype_matrix(cov, num_subjects=10000)
        phenotypes = phenotype_vector(genotypes, betas)
        fit_regression(genotypes, cov)
        
        



if __name__ == '__main__':
    main()