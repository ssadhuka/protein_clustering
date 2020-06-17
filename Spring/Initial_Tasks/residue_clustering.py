#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:30:08 2020

@author: shuvomsadhuka
"""

from biopandas.pdb import PandasPdb
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import sys
from collections import Counter
pd.set_option('display.max_rows', 150)

class ResidueClustering(object):
    
    def __init__(self, protein_df, mutation_df, dicti):
        dfs = self.make_protein(protein_df, mutation_df, dicti)
        self.all_pos_df = dfs[0]
        self.muts_df = dfs[1]
        self.set_properties()
        
    def make_protein(self, protein_df, mutation_df, dicti):
        """Concatenates two files (a protein index and residue-mutation score file) into a single Pandas mutation_dfframe object.
        INPUTS
        =======
        protein_df: Pandas mutation_dfFrame
          A DataFrame giving the residue numbers, locations, and amino acid information for each residue in a specific protein
        mutation_df: Pandas mutation_dfFrame
          A mutation_dfFrame storing all the mutations with residue numbers
        dicti: dictionary
            a dictionary storing the 3-to-1 letter encodings of amino acids
        RETURNS
        ========
        df_tot: a concatenated DataFrame of all mutations and corresponding locations, residue numbers, etc.
        """
        
        # extract residue number and name of each mutation from column
        mutation_df['residue_number'] = mutation_df[0].str.strip().str[3:].str[:-1]
        mutation_df = mutation_df[mutation_df['residue_number'].astype(int).isin(list(set(protein_df['residue_number'])))]
        mutation_df['residue_name_convert'] = mutation_df[0].str.strip().str[2:3]
        
        # get x, y, z coords of each residue as mean coordinate of all atoms making up residue
        for i in ['x', 'y', 'z']:
            coord = i + '_coord'
            pos_df = protein_df.groupby(['residue_number', 'chain_id'], as_index=False)[coord].mean()
            #print(pos_df)
            #df = df.join(pos_df.set_index('residue_number'), on=['residue_number', 'chain_id'], how='right', lsuffix='_left', rsuffix='_right')
            protein_df = pd.merge(protein_df, pos_df, on=['residue_number', 'chain_id'])
        #print(df_new)
    	#df['chain_id'] = df['chain_id_right']
        
        protein_df = protein_df.rename(columns={'x_coord_y': 'x_coord_right', 'y_coord_y': 'y_coord_right', 'z_coord_y': 'z_coord_right'})
        
        # drop all columns except residue name, number, and coordinates
        protein_df = protein_df[['residue_name', 'residue_number', 'x_coord_right', 'y_coord_right', 'z_coord_right', 'chain_id']]
        protein_df = protein_df.drop_duplicates(subset=['x_coord_right', 'y_coord_right', 'z_coord_right'], keep="first") 
        # 3-to-1 amino acid string conversion
        #df['residue_name_convert'] = list(sequence['residue_name'])
        #print(df)
    
        # type casting
        mutation_df = mutation_df.rename(columns={1: "score", 2: "p-val"})
        
        #print(df[df['residue_number'] == 7])
        mutation_df = mutation_df.astype({'residue_number': 'int64'})
    
        # merge datasets
        df_tot = pd.merge(protein_df, mutation_df, on="residue_number")
        df_tot['residue_name'] = df_tot.replace({'residue_name': dicti})
        
        df_tot = df_tot.drop_duplicates(subset=['x_coord_right', 'y_coord_right', 'z_coord_right'], keep='first')
        df_tot['dummy'] = (df_tot['residue_name'] == df_tot['residue_name_convert']).astype(int)
        df_tot = df_tot[df_tot['dummy'] == 1]
        
        df_tot['risk_prob'] = df_tot['residue_number']
        df_tot.loc[df_tot['score'] < 0, 'risk_prob'] = df_tot['p-val']/2
        df_tot.loc[df_tot['score'] > 0, 'risk_prob'] = 1-df_tot['p-val']/2
        return(protein_df, df_tot)
    
    
    def set_properties(self):
        l1 = list(self.muts_df['residue_name'])
        l2 = list(self.muts_df['residue_name_convert'])
        self.num_match = len([x for x, j in zip(l1, l2) if x == j])
        
        self.var_risk = len(list(set(self.muts_df[self.muts_df['score'] > 0]['score'])))
        self.var_prot = len(list(set(self.muts_df[self.muts_df['score'] < 0]['score'])))
        self.prob_vec = np.array(self.muts_df['risk_prob'])
        
        
    def set_dist_vec(self, local_distance, probability, risk, t=0):
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
        if risk:
            prob_vec = self.prob_vec
        else:
            prob_vec = 1 - self.prob_vec
        
        df_dist = self.muts_df.drop_duplicates(subset=['x_coord_right', 'y_coord_right', 'z_coord_right'], keep='first')
        
        muts = np.array(df_dist[df_dist.columns[2:5]])
        #print(muts)
        X = euclidean_distances(muts, muts)
        
        if probability:
            prob_mat = np.outer(prob_vec, prob_vec)
            X = np.multiply(X, prob_mat)
        
        if local_distance:
            X = -np.exp(X)/(2*t**2)
        
        X = X[~np.eye(X.shape[0],dtype=bool)].reshape(X.shape[0],-1)
        
        if risk and probability: self.emp_risk = round(np.sum(X)/np.sum(prob_mat))
        if risk and not probability: self.emp_risk = round(np.mean(X), 3)
        if not risk and probability: self.emp_prot = round(np.sum(X)/np.sum(prob_mat), 3)
        if not risk and not probability: self.emp_prot = round(np.mean(X), 3)
            
    
    def mean_dist(self, prots, local_distance, t=0):
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
        X = euclidean_distances(prots, prots)
        #if probability:
        #    prob_mat = np.outer(prob_vec, prob_vec)
        #    X = np.multiply(X, prob_mat)
        
        if local_distance:
            X = -np.exp(X)/(2*t**2)
        
        X = X[~np.eye(X.shape[0],dtype=bool)].reshape(X.shape[0],-1)
        return(round(np.mean(X), 3))
    
    
    def prob_mean_dist(self, prots, local_distance, prob_vec, t=0):
        # get pairwise distances of all rows
        #print(prots)
        X = euclidean_distances(prots, prots)
        prob_mat = np.outer(prob_vec, prob_vec)
        X = np.multiply(X, prob_mat)
        
        if local_distance:
            X = -np.exp(X)/(2*t**2)
        
        X = X[~np.eye(X.shape[0],dtype=bool)].reshape(X.shape[0],-1)
        return(round(np.sum(X)/np.sum(prob_mat), 3))
    
    
    def make_stratum(self, k, probability):
        """Permutes residues within a k-mer stratum.  For example, if k=2, permutes all mutated di-mers in
        a protein.
        INPUTS
        =======
        df: array
          An nx5 DataFrame giving the coordinates, directionality, and p-value on each mutation in a protein.
        k: int
          An integer for the desired stratum to permute.  E.g. k=2 permute di-mer mutations.
        RETURNS
        ========
        a DataFrame of permuted residues within k-mer stratum
        """
        
        # subset dataframe by stratum
        df_k = self.all_pos_df[self.all_pos_df['stratum'] == k]
        muts_df_k = self.muts_df[self.muts_df['stratum'] == k]
        num_muts = len(list(muts_df_k['score']))
        
        
        if not df_k.empty:
            # convert coordinates of mutations to matrix
            x = df_k.iloc[:,[2,3,4]].to_numpy()
            
            # residues corresponding to the same mutation should fit onto one row in the matrix - i.e.
            # group residues according to the mutation controlling them
            strat_fixed_mat = x[0::k,:]
            for i in range(1, k):
                strat_fixed_mat = np.concatenate([strat_fixed_mat, x[i::k,:]], axis=1)
            
            # permute rows of the stratified matrix - equivalent to permuting mutations
            np.random.shuffle(strat_fixed_mat)
            
            # rebuild matrix by splitting matrix k times corresponding to k-mer
            strat_mat_split = np.split(strat_fixed_mat, k, axis=1)
            post_split = np.concatenate(strat_mat_split, axis=0)
        
            # convert to dataframe
            df_new = pd.DataFrame(post_split)
            
            # the assignment groups coordinates by residue number
            assignments = [n for n in range(int(x.shape[0]/k))] * k
            df_new['assignments'] = assignments
            
            # sort by assignment - want residues from same mutation to be consecutive rows in matrix
            df_new = df_new.sort_values(by=['assignments'])
            df_new = df_new.head(num_muts)
        
        # assign scores
        df_new['score'] = list(muts_df_k['score'])
        df_new['risk_prob'] = list(muts_df_k['risk_prob'])
        return(df_new)
    
    
    def stratified_permutation(self, probability):
        """Returns stratified permutation across all k-mers in a protein.
        INPUTS
        =======
        df: array
          An nx5 DataFrame giving the coordinates, directionality, and p-value on each mutation in a protein.
        RETURNS
        ========
        a DataFrame of permuted residues within all k-mer strata.
        """
        
        # mark strata by counting k-mer number - use dictionary comprehension
        strata = Counter(list(self.all_pos_df['residue_number']))
        self.muts_df['stratum'] = self.muts_df['residue_number'].map(strata)
        self.all_pos_df['stratum'] = self.all_pos_df['residue_number'].map(strata)

        # build out stratified permutation
        df_perm = self.make_stratum(int(list(set(self.all_pos_df['stratum']))[0]), probability)
       
        #print(list(set(df['stratum'])))
        if len(list(set(self.muts_df['stratum']))) > 1:
            for k in list(set(self.muts_df['stratum']))[1:]:
                df_perm = df_perm.append(self.make_stratum(k, probability))
        return(df_perm)
    
    
    def get_prob_vec(self, df, risk):
        if risk:
            return(df['risk_prob'].to_numpy())
        else:
            return(1.0 - df['risk_prob'].to_numpy())
        
        
    def make_box(self):
        min_x, max_x = self.muts_df.min(axis=0)['x_coord_right'], self.muts_df.max(axis=0)['x_coord_right']
        min_y, max_y = self.muts_df.min(axis=0)['y_coord_right'], self.muts_df.max(axis=0)['y_coord_right']
        min_z, max_z = self.muts_df.min(axis=0)['z_coord_right'], self.muts_df.max(axis=0)['z_coord_right']
        
        self.all_pos_df = self.all_pos_df[(self.all_pos_df['x_coord_right'] >= min_x) & (self.all_pos_df['x_coord_right'] <= max_x)]
        self.all_pos_df = self.all_pos_df[(self.all_pos_df['y_coord_right'] >= min_y) & (self.all_pos_df['y_coord_right'] <= max_y)]
        self.all_pos_df = self.all_pos_df[(self.all_pos_df['z_coord_right'] >= min_z) & (self.all_pos_df['z_coord_right'] <= max_z)]
        
        
    def run_permutation(self, num_runs, local_distance, probability, trials=0, t=0):
        dist_prots = []
        dist_risks = []
        
        len_prots = len(self.muts_df[self.muts_df['score'] > 0])
        len_tots = len(self.muts_df[self.muts_df['score'] > 0]) + len(self.muts_df[self.muts_df['score'] < 0])
        #x = muts_df[[2:5]].to_numpy()
        #x = muts_df.iloc[:,[2,3,4]].to_numpy()
        #prots = x[0:len_prots,:]
        #prots = np.unique(prots, axis=0)
        #print(mean_dist(prots))
        #risks = x[len_prots:len_tots,:]
        #risks = np.unique(risks, axis=0)
        #emp_risk = mean_dist(risks)
        
        for i in range(num_runs):
            
            #np.random.shuffle(x)
            x = self.stratified_permutation(probability=probability)
            #x = muts_df
            if probability:
                prob_vec_risk = self.get_prob_vec(x, risk=True)
                prob_vec_prot = self.get_prob_vec(x, risk=False)
                dist_prots.append(self.prob_mean_dist(x, local_distance, prob_vec_prot, t=t))
                dist_risks.append(self.prob_mean_dist(x, local_distance, prob_vec_risk, t=t))
                 
            else:    
                x_prots = x[x['score'] > 0].iloc[:, [0,1,2]].to_numpy()
                x_risks = x[x['score'] < 0].iloc[:, [0,1,2]].to_numpy()
                    
                dist_prots.append(self.mean_dist(x_prots, local_distance, probability))
                dist_risks.append(self.mean_dist(x_risks, local_distance, probability))
                #prots = x[0:len_prots,:]
                #prots = np.unique(prots, axis=0)
                #dist_prots.append(mean_dist(prots))
               
                #risks = x[len_prots:len_tots,:]
                #risks = np.unique(risks, axis=0)
                #dist_risks.append(mean_dist(risks))
        
        #print(np.array(self.muts_df['risk_prob']))
        dist_risks.append(self.emp_risk)
        dist_prots.append(self.emp_prot)
        #print(prob_vec_prot)
        #print(dist_prots)

        ret1 = stats.percentileofscore(dist_risks, self.emp_risk, kind='mean')
        ret2 = stats.percentileofscore(dist_prots, self.emp_prot, kind='mean')
        #if (ret1 < 5 or ret2 < 5) and trials==0:
        #        print(ret1, ret2)
        #        ret = self.run_permutation(5000, local_distance, probability, trials=1, t=t)
        #        ret1 = ret[0]
        #        ret2 = ret[1]
        #elif (ret1 < 0.5 or ret2 < 0.5) and trials==1:
        #        ret = self.run_permutation(50000, local_distance, probability, trials=2, t=t)
        #        ret1 = ret[0]
        #        ret2 = ret[1]
        
        return(ret1, ret2)
        
    
    def execute(self, local_distance=False, probability=False, fix_residues=False, t=0):
        if fix_residues:
            self.all_pos_df = self.muts_df
        
        #all_pos_df = self.make_box(all_pos_df, muts_df)
        
        if not (self.var_prot < 2 or self.var_risk < 2):
            self.set_dist_vec(local_distance=local_distance, probability=probability, risk=True, t=t)
            self.set_dist_vec(local_distance=local_distance, probability=probability, risk=False, t=t)
            print(self.emp_risk)
            
            rets = self.run_permutation(100, local_distance=local_distance, probability=probability, t=t)
            ret1 = rets[0]
            ret2 = rets[1]
        else:
            ret1 = 'EMPTY'
            ret2 = 'EMPTY'
        
        return(ret1, ret2)
        
        
        
        
        
        
        