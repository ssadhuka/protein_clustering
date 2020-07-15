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



class SimP(object):
    def __init__(self, num_subjects, num_mutations, pct_causal, alpha, pi, gamma, Z_mat, dist_kernel):
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
        
    def make_pis_and_gammas(self, pi, gamma, pct_causal):
        indices = np.random.choice(np.arange(len(pi)), replace=False, size=int(len(pi) * (1 - pct_causal)))
        for index in indices: 
            #pi[index] = 0
            gamma[index] = 0
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
        GD = self.genotypes.T @ sqrtm(self.dist_kernel) @ self.gamma
        self.phenotypes = (self.covariates.T @ self.alpha + self.genotypes_fixed.T @ self.pi) + GD + noise
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
    


    







