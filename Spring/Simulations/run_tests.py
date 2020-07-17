#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:53:37 2020

@author: shuvomsadhuka
"""


from utils import *
from tau_pi_test import *
from simulate_phenotypes import SimP
import numpy as np
import math
from scipy.special import logit



def main(var_ratio, pi_sum_ratio, num_muts, num_indivs):
    protein_id = 'Q01113'
    protein_structure = make_protein_df(protein_id)
    protein_structure = select_mutations(protein_structure, num_muts)
    
    dists = distances_pairwise(protein_structure)
    #distance_kernel = covariance_matrix(dists, 6)
    distance_kernel = np.eye(100)
    
    p_pis = []
    p_gammas = []
    pi_fixed, tau_sq = set_pi_and_tau(var_ratio, pi_sum_ratio, 0.25, eps=1)
    alpha = [0.1, 0.5, 0.5]
    prevalence = np.exp(logit(alpha[0]))
    
    
    for i in range(500):
        try:
            #gamma = np.random.normal(0, np.sqrt(tau_sq), size=100)
            #pi_fixed = [0]
            #gamma = [0 for i in range(100)]
            
            #for i in range(100):
            simp = SimP(num_subjects=num_indivs, num_mutations=num_muts, pct_causal = 0.5, alpha=alpha, pi=[pi_fixed], 
                        gamma=gamma, 
                        Z_mat = np.ones((100,1)),
                        dist_kernel = distance_kernel)
            
            
            # hold disease prevalance constant
            threshold = np.sort(simp.phenotypes)[math.floor(len(simp.phenotypes) * prevalence)]
            phenotypes = np.where(simp.phenotypes > threshold, 1, 0)
            #phenotypes = np.where(simp.phenotypes > 0.5, 1, 0)
            genotypes = simp.genotypes
            covariates = simp.covariates
            
            Xs_0 = covariates
            #Xs_1 = np.concatenate((covariates, genotypes))
            Xs_2 = np.concatenate((covariates,(genotypes.T @ np.ones((100,1))).T))
            #Xs_3 = (genotypes.T @ np.random.multivariate_normal(np.zeros(200), np.eye(200)).reshape(100, 2)).T
            #Xs_4 = genotypes
            #print(Xs_1)
            
            lr_alpha = run_logistic(Xs_0, phenotypes.flatten())
            lr_pi = run_logistic(Xs_2, phenotypes.flatten())
            p_pi = test_pi(Xs_0, phenotypes.flatten(), genotypes, covariates, lr_alpha)
            p_tau = test_tau(Xs_2, phenotypes.flatten(), genotypes,  covariates, distance_kernel, lr_pi,
                     temp_dir = 'temp_files/pi_0_gamma_0.txt',
                     out_dir = 'type1_sims/pi_0_gamma_0.txt')
            print(p_pi, p_tau)
            p_pis.append(p_pi)
            p_gammas.append(p_tau)
            
        except Exception as e:
            print(e)
            continue
        
    df = pd.DataFrame(list(zip(p_pis, p_gammas)), columns=['p_pi', 'p_gamma'])
    #out_dir = 'power/both_nonzero/var_ratio_' + str(var_ratio) + '_pi_ratio_' + str(pi_sum_ratio) + '/'
    #out_dir += str(num_muts) + 'muts_' + str(num_indivs) + 'indiv.csv'
    out_dir = 'type1_error/type1sims_1000indiv_100muts.csv'
    df.to_csv(out_dir)
    
if __name__ == '__main__':
    #var_ratios = [(i+1)*0.01 for i in range(5)]
    var_ratios = [0.00]
    pi_sum_ratios = [0.75]
    indivs = [1000]
    
    for vr in var_ratios:
        for psr in pi_sum_ratios:
            for indiv in indivs:
                main(var_ratio=vr, pi_sum_ratio=psr, num_muts=100, num_indivs=indiv)
                
            