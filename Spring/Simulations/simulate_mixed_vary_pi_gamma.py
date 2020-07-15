#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:45:33 2020

@author: shuvomsadhuka
"""
import numpy as np
from scipy.special import logit, expit
from scipy.linalg import sqrtm
from scipy.stats import chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
import os
import pandas as pd
from biopandas.pdb import PandasPdb



class SimulatePhenotypes(object):
    def __init__(self, num_subjects, num_mutations, pct_causal, alpha, pi):
        self.num_subjects = num_subjects
        self.num_mutations = num_mutations
        self.alpha = np.array(alpha).reshape(len(alpha), 1)
        self.make_pis(pi, pct_causal)
        self.make_mafs()
        self.make_genotypes()
        self.make_covariates()
        self.simulate_phenotypes()
    
    def make_mafs(self):
        self.maf = np.random.uniform(0, 0.5, size=self.num_mutations)
        
    def make_pis(self, pi, pct_causal):
        indices = np.random.choice(np.arange(len(pi)), replace=False, size=int(len(pi) * (1 - pct_causal)))
        for index in indices: pi[index] = 0
        self.pi = np.array(pi).reshape(len(pi), 1)
        
    def make_covariates(self):
        X1 = np.random.normal(0, 1, size=self.num_subjects)
        X2 = np.random.binomial(1, 0.5, size=self.num_subjects)
        X0 = np.ones(shape=self.num_subjects)
        self.covariates = np.concatenate([X0, X1, X2]).reshape((3, self.num_subjects))
        
    def make_genotypes(self):
        self.genotypes = np.array(list(map(lambda freq: np.random.binomial(2, freq, size=self.num_subjects), self.maf)))
        
    def simulate_phenotypes(self):
        noise = np.random.normal(size=self.num_subjects).reshape((self.num_subjects, 1))
        self.phenotypes = (self.covariates.T @ self.alpha + self.genotypes.T @ self.pi) + noise
        self.phenotypes = expit(self.phenotypes)
        #print(self.pi)
    
    
    
class ScoreTest(object):
    def __init__(self):
        pass
    

class SKAT(object):
    def __init___(self):
        pass
    

class SKAT_O(object):
    def __init__(self):
        pass
    
    
class Burden(object):
    def __init__(self):
        pass
    

def run_logistic(X, y):
    logitr = LogisticRegression(penalty='none', fit_intercept=False)
    logitr.fit(X.T, y)
    return(logitr)
    
    # do a sanity check that predict_proba agrees with just straight calculated mu hats
    '''mu_hat = logitr.predict_proba(X.T)[:,1]
    sigma_hat = mu_hat * (1 - mu_hat)
    D_hat = np.diag(sigma_hat)
    
    U_pi = genotypes @ (phenotypes.flatten() - mu_hat.flatten())
    
    inner = (D_hat @ covariates.T) @ np.linalg.inv(covariates @ D_hat @ covariates.T) @ (covariates @ D_hat)
    inner = D_hat - inner
    Sigma = genotypes @ inner @ genotypes.T
    U_pi_standardized = np.linalg.inv(sqrtm(Sigma)) @ U_pi
    
    
    mc = 0
    print(U_pi_standardized)
    for i in range(10000):
        ex = np.random.multivariate_normal(mean=np.zeros(len(U_pi_standardized)), cov=np.identity(len(U_pi_standardized)))
        
        if np.linalg.norm(ex) > np.linalg.norm(U_pi_standardized):
            mc += 1
            
    print(mc/10000)
    '''
            
def test_pi(X, y, genotypes, covariates, lr_alpha):
    mu_hat = lr_alpha.predict_proba(X.T)[:,1]
    sigma_hat = mu_hat * (1 - mu_hat)
    D_hat = np.diag(sigma_hat)
    
    U_pi = genotypes @ (y.flatten() - mu_hat.flatten())
    
    inner = (D_hat @ covariates.T) @ np.linalg.pinv(covariates @ D_hat @ covariates.T) @ (covariates @ D_hat)
    inner = D_hat - inner
    Sigma = genotypes @ inner @ genotypes.T
    U_pi_standardized = np.linalg.pinv(sqrtm(Sigma)) @ U_pi
    U_pi_chi = np.sum(np.square(U_pi_standardized))
    
    return(1-chi2.cdf(U_pi_chi, len(U_pi_standardized)))
    

def test_tau(X, y, genotypes, covariates, distance_kernel, lr_pi):
    mu_hat = lr_pi.predict_proba(X.T)[:,1]
    
    kernel = genotypes.T @ distance_kernel @ genotypes
    
    

def distances_pairwise(df):
    muts = df[['x', 'y', 'z']].to_numpy()
    return(euclidean_distances(muts, muts))


def covariance_matrix(dists):
    cov = 1/(dists+1)
    return(cov)




pi = [0.1 for i in range(10)]   


ps = []
for i in range(1000):
    SP = SimulatePhenotypes(num_subjects=500, num_mutations=10, pct_causal = 1.0, alpha=[0.1, 0.5, 0.5], pi=pi)
    print(SP.pi)
    phenotypes = np.where(SP.phenotypes > 0.5, 1, 0)
    genotypes = SP.genotypes
    covariates = SP.covariates
    
    Xs_0 = covariates
    Xs_1 = np.concatenate((covariates, genotypes))
    lr_alpha = run_logistic(Xs_0, phenotypes.flatten())
    lr_pi = run_logistic(Xs_1, phenotypes.flatten())
    ps.append(test_pi(Xs_0, phenotypes.flatten(), genotypes, covariates, lr_alpha))
# need to test null model with pi = 0 and gamma = 0 and then null model with only gamma = 0
# first test is U_pi, second test is S_tau^2

theoretical = [(i+1)/1000 for i in range(1000)]
theoretical = -np.log10(theoretical)
ps = -np.log10(np.sort(ps))

plt.scatter(theoretical, ps)
p1= [0,3.5]
p2=[0, 3.5]
plt.plot(p1,p2, c='r')
plt.xlabel('Theoretical')
plt.ylabel('Observed')



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
