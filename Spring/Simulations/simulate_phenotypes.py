#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:45:33 2020

@author: shuvomsadhuka

This file simulates phenotypes for a specified number of individuals (n) given some parameters, including 
their genotypes, a protein structure-based distance kernel, and relative weights on fixed and random
effects.  We assume a mixed effect model, similar to "A Unified Mixed-Effects Model for Rare-Variant Association 
in Sequencing Studies" (2013) with modifications for the kernel.

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


import numpy as np
from scipy.special import expit
#from scipy.linalg import sqrtm
#from scipy.stats import chi2
#from sklearn.linear_model import LogisticRegression
#from sklearn.metrics.pairwise import euclidean_distances
#import matplotlib.pyplot as plt
#import os
#import pandas as pd
#from biopandas.pdb import PandasPdb


class SimP(object):
    def __init__(self, num_subjects, num_mutations, pct_causal, alpha, pi, gamma, Z_mat, 
                 dist_kernel, tau_sq, phenotype_type):
        self.phenotype_type = phenotype_type
        self.tau_sq = tau_sq
        self.num_subjects = num_subjects
        self.num_mutations = num_mutations
        self.alpha = np.array(alpha).reshape(len(alpha), 1)
        self.Z_mat = Z_mat
        self.dist_kernel = dist_kernel
        self.make_pis_and_gammas(pi, gamma, pct_causal)
        self.make_mafs()
        self.make_genotypes()
        self.make_covariates()
        self.simulate_phenotypes()
    
    def make_mafs(self):
        self.maf = np.random.uniform(0, 0.5, size=self.num_mutations)
        self.combined_maf = np.mean(self.maf)
        
    def make_pis_and_gammas(self, pi, gamma, pct_causal):
        indices = np.random.choice(np.arange(len(pi)), replace=False, size=int(len(pi) * (1 - pct_causal)))
        for index in indices: 
            #pi[index] = 0
            #gamma[index] = 0
            pass
        self.pi = np.array(pi).reshape(len(pi), 1)
        self.gamma = np.array(gamma).reshape(len(gamma), 1)
        
    def make_covariates(self):
        X1 = np.random.normal(0, 1, size=self.num_subjects)
        X2 = np.random.binomial(1, 0.5, size=self.num_subjects)
        X0 = np.ones(shape=self.num_subjects)
        self.covariates = np.concatenate([X0, X1, X2]).reshape((3, self.num_subjects))
        
    def make_genotypes(self):
        self.genotypes = np.array(list(map(lambda freq: np.random.binomial(2, freq, size=self.num_subjects), self.maf)))
        self.genotypes_fixed = (self.genotypes.T @ self.Z_mat).T
      
    def simulate_phenotypes(self):
        noise = np.random.normal(size=self.num_subjects).reshape((self.num_subjects, 1))
        
        GD = self.genotypes.T @ self.cholesky_decomp(self.dist_kernel) @ self.gamma
        self.phenotypes = (self.covariates.T @ self.alpha + self.genotypes_fixed.T @ self.pi) + GD + noise
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
    


    







