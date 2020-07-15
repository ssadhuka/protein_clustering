#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:53:37 2020

@author: shuvomsadhuka
"""


from utils import *
from tau_pi_test import *
from simulate_phenotypes import SimP



def main():
    protein_id = 'Q01113'
    protein_structure = make_protein_df(protein_id)
    protein_structure = select_mutations(protein_structure, 100)
    
    dists = distances_pairwise(protein_structure)
    distance_kernel = covariance_matrix(dists, 6)
    #distance_kernel = np.eye(10)
    
    p_pis = []
    p_gammas = []
    for i in range(1000):
        try:
            pi = [0]
            gamma = [-0.1, 0.1, -0.05, -0.05, 0, 0.01, 0.17, -0.02, 0.1, 0] * 10
            #gamma = [0 for i in range(100)]
            
            #for i in range(100):
            simp = SimP(num_subjects=500, num_mutations=100, pct_causal = 0.5, alpha=[0.1, 0.5, 0.5], pi=pi, gamma=gamma, 
                        Z_mat = np.ones((100,1)),
                        dist_kernel = distance_kernel)
            
            phenotypes = np.where(simp.phenotypes > 0.5, 1, 0)
            genotypes = simp.genotypes
            covariates = simp.covariates
            
            Xs_0 = covariates
            Xs_1 = np.concatenate((covariates, genotypes))
            Xs_2 = np.concatenate((covariates,(genotypes.T @ np.ones((100,1))).T))
            Xs_3 = (genotypes.T @ np.random.multivariate_normal(np.zeros(200), np.eye(200)).reshape(100, 2)).T
            Xs_4 = genotypes
            #print(Xs_1)
            
            lr_alpha = run_logistic(Xs_0, phenotypes.flatten())
            lr_pi = run_logistic(Xs_2, phenotypes.flatten())
            pi = test_pi(Xs_0, phenotypes.flatten(), genotypes, covariates, lr_alpha)
            tau = test_tau(Xs_2, phenotypes.flatten(), genotypes,  covariates, distance_kernel, lr_pi,
                     temp_dir = 'temp_files/pi_0_gamma_0.txt',
                     out_dir = 'type1_sims/pi_0_gamma_0.txt')
            p_pis.append(pi)
            p_gammas.append(tau)
            
        except Exception as e:
            print(e)
            continue
        
    df = pd.DataFrame(list(zip(p_pis, p_gammas)), columns=['p_pi', 'p_gamma'])
    df.to_csv('power/nonzero_gamma_only/100muts_500indiv_gamma_mixed.csv')
    
if __name__ == '__main__':
    main()