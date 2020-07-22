#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:31:37 2020

@author: shuvomsadhuka

This file simulates phenotypes for a specified number of individuals (n) given some parameters, including 
their genotypes, a protein structure-based distance kernel, and relative weights on fixed and random
effects.  We assume a mixed effect model, similar to "A Unified Mixed-Effects Model for Rare-Variant Association 
in Sequencing Studies" (2013) with modifications for the distance kernel.

We developed a statistical procedure to test whether mutations in a certain gene tend to cluster
in their corresponding protein structures.  We use the simulated data from this file to evaluate
the power and error rates of our procedure.  Code for the actual test is given in tau_pi_test.py.

INPUTS
=======
num_subjects: integer
  An integer specifying the number of individuals in the sample.
num_mutations: integer
  We consider mutations on one gene at a time, so num_mutations specifies the number of loci that are mutated
  in the gene of choice.
pct_causal: float
  Percent of mutations in the gene that are actually causal.
alpha: vector
  Vector of coefficients from random covariates.
pi: scalar
  Fixed effect coefficient, equivalent to burden beta coefficient.
gamma: vector
  Random effect vector, sampled from N(0, tau^2).
Z_mat: vector
  Weight vector (could be matrix, in which case pi becomes a vector) of weights of loci for fixed effects.
dist_kernel: matrix
  Symmetric matrix f(D) with each element f(D)_ij = f(d_ij), or a function of the distance of residue i from
  residue j.
tau_sq: scalar
  Variance of gamma, pre-selected to keep variance explained by pi vs gamma constant.


"""

from scipy.special import expit
from scipy.spatial import distance_matrix
import pandas as pd
import numpy as np


def distance_to_domain(df, causal_domain, i):
    df['distance'] = (df['x'].subtract(causal_domain[i][0]))**2 + (df['y'].subtract(causal_domain[i][1]))**2 
    df['distance'] += (df['z'].subtract(causal_domain[i][2]))**2 
    df['distance'] = np.sqrt(df['distance'])
    return(df)


def select_causal_domains(df, radius, num_domains):
    its = 0
    min_dist = 0
    if num_domains == 1:
        causal_domain = df.sample(n=num_domains).to_numpy()
        min_dist = 2*radius + 5
        
    while (its < 10) and (min_dist < 2*radius + 1) and (num_domains > 1):
        causal_domain = df.sample(n=num_domains).to_numpy()
        dists = distance_matrix(causal_domain, causal_domain)
        min_dist = np.min(dists[np.nonzero(dists)])
        its += 1

    if min_dist < 2*radius + 1:
        raise ValueError('Protein structure cannot support %s domains' % num_domains)
    for n in range(num_domains):
        structure = distance_to_domain(df, causal_domain, n)
        structure['domain'].mask(structure['distance'] < radius, n+1, inplace=True)
    return(structure)

 
def select_domains(df, radius, total_muts, prob_causal, num_domains):
    structure = df[['x', 'y', 'z']]
    structure['domain'] = [0 for i in range(len(structure))]
    
    structure = select_causal_domains(structure, radius, num_domains)

    num_domain_muts = (structure.domain.values != 0).sum()
    non_domain_muts = structure[structure['domain'] == 0].sample(total_muts-num_domain_muts)
    non_domain_muts['causal'] = [0 for i in range(total_muts - num_domain_muts)]
    non_domain_muts['maf'] = np.random.uniform(0, 0.5, size=total_muts-num_domain_muts)  
    
    domain_muts_tot = structure[structure['domain'] == 1]
    num_domain_muts_tot = (structure.domain.values == 1).sum()
    
    causals = []
    domain_muts_tot['maf'] = np.random.uniform(0, 0.5, size=num_domain_muts_tot)
    domain_muts_tot = domain_muts_tot.sample(frac=1)
    for i in range(int(prob_causal*num_domain_muts_tot)): causals.append(1)
    for j in range(num_domain_muts_tot - int(prob_causal*num_domain_muts_tot)): causals.append(0)
    domain_muts_tot['causal'] = causals
    
    for k in range(1, num_domains):
        num_domain_muts = (structure.domain.values == k+1).sum()
        domain_muts = structure[structure['domain'] == k+1]
        domain_muts['maf'] = np.random.uniform(0, 0.5, size=num_domain_muts)
        domain_muts = domain_muts.sample(frac=1)
        
        causals = []
        for i in range(int(prob_causal*num_domain_muts)): causals.append(k+1)
        for j in range(num_domain_muts - int(prob_causal*num_domain_muts)): causals.append(0)
        domain_muts['causal'] = causals
        domain_muts_tot = pd.concat([domain_muts_tot, domain_muts])
    
    muts = pd.concat([domain_muts_tot, non_domain_muts]).sort_index()    
    return(muts)


def set_betas(df, var_to_explain, domain_vars=None, eps=1):
    df_causal = df[df['causal'] != 0]
    
    scaled_sum_betas = var_to_explain/(1 - var_to_explain)
    #domain_sum_betas = scaled_sum_betas * domain_vars
    scaled_beta = np.sqrt(scaled_sum_betas/np.sum((1 - df_causal['maf']) * (df_causal['maf'])))
    df['beta'] = np.where(df['causal'] != 0, scaled_beta, 0)
    df['beta'].mask(df['causal'] == 2, -df['beta'], inplace=True)
    return(df)


class SimD(object):
    def __init__(self, num_subjects, df, alpha, phenotype_type):
        self.num_subjects = num_subjects
        self.df = df
        self.alpha = alpha
        self.phenotype_type = phenotype_type
        self.set_variables()
        self.make_genotypes()
        self.make_covariates()
        self.simulate_phenotypes()
    
    def set_variables(self):
        self.mafs = np.array(self.df['maf'])
        self.betas = np.array(self.df['beta'])
    
    def make_covariates(self):
        X1 = np.random.normal(0, 1, size=self.num_subjects)
        X2 = np.random.binomial(1, 0.5, size=self.num_subjects)
        X0 = np.ones(shape=self.num_subjects)
        self.covariates = np.concatenate([X0, X1, X2]).reshape((3, self.num_subjects))
        
    def make_genotypes(self):
        self.genotypes = np.array(list(map(lambda freq: np.random.binomial(2, freq, size=self.num_subjects), self.mafs)))
      
    def simulate_phenotypes(self):
        noise = np.random.normal(size=self.num_subjects).reshape((self.num_subjects, 1)).T
        
        Gbeta = self.genotypes.T @ self.betas
        self.phenotypes = self.covariates.T @ self.alpha + Gbeta + noise
        print(self.phenotypes.shape)
        if self.phenotype_type == 'logistic':
            self.phenotypes = expit(self.phenotypes)
        #print(self.pi)
        
    @staticmethod
    def cholesky_decomp(cor_matrix):
        try:
            L = np.linalg.cholesky(cor_matrix)
        except np.linalg.linalg.LinAlgError:
            l, Q = np.linalg.eigh(cor_matrix)
            l[l < 0] = 0
            L = np.dot(Q, np.diag(np.sqrt(l)))
            #renormalize diagonal to be 1
            D2 = np.diag(1/np.sqrt(np.diag(np.dot(L, L.transpose()))))
            L = np.dot(D2, L)
        return(L)
    
    
    