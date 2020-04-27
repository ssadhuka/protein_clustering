#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
    
    print(data)
    
    # get x, y, z coords of each residue as mean coordinate of all atoms making up residue
    for i in ['x', 'y', 'z']:
        coord = i + '_coord'
        pos_df = df.groupby('residue_number', as_index=False)[coord].mean()
        df = df.join(pos_df.set_index('residue_number'), on='residue_number', how='right', lsuffix='_left', rsuffix='_right')
    
    # drop all columns except residue name, number, and coordinates
    df = df[['residue_name', 'residue_number', 'x_coord_right', 'y_coord_right', 'z_coord_right']]
    df = df.drop_duplicates(subset='residue_number', keep="first")
    
    # 3-to-1 amino acid string conversion
    #df['residue_name_convert'] = list(sequence['residue_name'])
    
    
    
    # type casting
    data = data.rename(columns={1: "score", 2: "p-val"})
    data = data.astype({'residue_number': 'int64'})
    
    # merge datasets
    df_tot = pd.merge(df, data, on="residue_number")
    #print(df_tot.count)
    return df_tot


def get_dist_vec(df, risk):
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
        df_dist = df[df['score'] < 0]
    
    else: 
        df_dist = df[df['score'] > 0]
    
    df_dist = df_dist.drop_duplicates(subset=['x_coord_right'], keep='first')
    
    # isolate coordinates of mutations
    muts = np.array(df_dist[df_dist.columns[2:5]])
    nrow = muts.shape[0]
    
    # get pairwise distances of all rows, remove if row is duplicated (i.e. don't want to calculate distance of mutation A to itself)
    X = euclidean_distances(muts, muts)
    X = np.ma.masked_equal(X,0)
    X = X.compressed().reshape(nrow,nrow-1)
    
    # distances vector is row means of distance matrix
    vec = np.mean(X, axis=1)
    return vec


def mean_dist(prots):
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
    prots_vec = np.mean(dists_prots, axis=1)
    
    Y = dists_prots
    Y = np.ma.masked_equal(Y,0)
    nrow = Y.shape[0]
    
    Y = Y.compressed().reshape(nrow,nrow-1)
    return(np.mean(Y))


def run_permutation(muts_df, df, emp_risk, emp_prot, num_runs):
    """Run a permutation test, permuting mutations of specific type over all possible residues
    INPUTS
    =======
    muts_df: Pandas DataFrame
        An nx3 array of coordinates of mutated residues
    df: Pandas DataFrame
      An nx3 array of coordinates of all residues on protein
    emp_risk:
        a float representing the mean distance between all risk mutations in our original dataset
    emp_prot:
        a float representing the mean distance between all protective mutations in our original dataset
    num_runs:
        number of times to permute the dataset
    RETURNS
    ========
    two percentile scores of empirically observed distance between mutations relatively to distribution from permutation test
    lower distances (possibly) indicate higher clustering
    """
    dist_prots = []
    dist_risks = []
    
    
    len_prots = len(muts_df[muts_df['score'] > 0])
    len_tots = len(muts_df[muts_df['score'] > 0]) + len(muts_df[muts_df['score'] < 0])

    
    x = df[["x_coord", "y_coord", "z_coord"]].to_numpy()
    for i in range(num_runs):
        np.random.shuffle(x)
        
        prots = x[0:len_prots,:]
        prots = np.unique(prots, axis=0)
        #print(prots.shape)
        dist_prots.append(mean_dist(prots))
       
        risks = x[len_prots:len_tots,:]
        risks = np.unique(risks, axis=0)
        #print(risks.shape)
        dist_risks.append(mean_dist(risks))
    
    dist_risks.append(emp_risk)
    dist_prots.append(emp_prot)
        
    # plot risk mutations and empirical distance observed
    # plt.hist(dist_risks, bins=50)
    # plt.axvline(x=emp_risk, ymin=0, ymax=1, color = 'r')
    # plt.xlabel('Average Distance between Risk Mutations')
    # plt.ylabel('Count')
    # plt.show()
    # plt.savefig('risk_cluster.eps')
    
    return(stats.percentileofscore(dist_risks, emp_risk), stats.percentileofscore(dist_prots, emp_prot))

def main(argv):
    ppdb = PandasPdb().read_pdb(argv[1])
    df = pd.DataFrame(ppdb.df['ATOM'])
    sequence = ppdb.amino3to1()
    
    data = pd.read_csv(argv[2], header = None, sep="\t")
    num_runs = argv[3]
    
    df = make_dataframe(df, data, sequence)
    
    risk = get_dist_vec(df, True)
    prot = get_dist_vec(df, False)
    
    mw = mannwhitneyu(risk, prot)
    #print(mw)
    
    perm = run_permutation(df, np.mean(risk), np.mean(prot), num_runs)
    #print(perm)
    
if __name__ == '__main__':
    main(sys.argv)
    
    