#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 18:27:20 2020

@author: shuvomsadhuka
"""


from utils import *
from tau_pi_test import *
from simulate_domains import *
from simulate_phenotypes import SimP
import numpy as np
import math
from scipy.special import logit
import warnings
warnings.filterwarnings("ignore")
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 100)



def main(var_ratio, pi_sum_ratio, num_muts, num_indivs):
    protein_id = 'Q01113'
    protein_structure = make_protein_df(protein_id)
    alpha = [0.1, 0.5, 0.5]
    prevalence = np.exp(logit(alpha[0]))
    domains = 2
    
    power = [0.5, 1, 2, 3]
    t = [3,6,9,12]
                
    #distance_kernel = np.eye(100)
    p_pis = []
    p_taus_dist = []        
    p_taus_id = []    
    p_taus = {}
    for alph in power:
        for t in ts:
            key = str(alph) + '_' + str(t)
            p_taus[key] = []
    p_taus['id'] = []
    p_taus['fixed'] = []
    
    for i in range(300):
        #try:
        muts = select_domains(protein_structure, radius=10, total_muts=100, prob_causal=0.5, 
                             num_domains=domains)
        muts = set_betas(muts, var_ratio)
        dists = distances_pairwise(muts)
        id_kernel = np.eye(100)

        simd = SimD(num_indivs, muts, alpha, phenotype_type='logistic')

        # hold disease prevalance constant
        threshold = np.sort(simd.phenotypes)[::-1][math.floor(len(simd.phenotypes) * prevalence)]
        phenotypes = np.where(simd.phenotypes > threshold, 1, 0)

        genotypes = simd.genotypes
        covariates = simd.covariates
        
        Xs_0 = covariates
        #Xs_1 = np.concatenate((covariates, genotypes))
        Xs_2 = np.concatenate((covariates,(genotypes.T @ np.ones((100,1))).T))
        #Xs_3 = (genotypes.T @ np.random.multivariate_normal(np.zeros(200), np.eye(200)).reshape(100, 2)).T
        #Xs_4 = genotypes
        #print(Xs_1)
        
        lr_alpha = run_logistic(Xs_0, phenotypes.flatten())
        lr_pi = run_logistic(Xs_2, phenotypes.flatten())
        p_pi = test_pi(Xs_0, phenotypes.flatten(), genotypes, covariates, lr_alpha,
                       regression = 'logistic')
        
        for alph in power:
            for t in ts:
                key = str(alph) + '_' + str(t)
                distance_kernel = covariance_matrix(dists, alph, t)
                p_tau_dist = test_tau(Xs_2, phenotypes.flatten(), genotypes,  covariates, distance_kernel, lr_pi,
                         temp_dir = 'temp_files/pi_0_gamma_0.txt',
                         out_dir = 'type1_sims/pi_0_gamma_0.txt',
                         regression = 'logistic')
                p_taus[key].append(float(p_tau_dist))
                
        
        p_tau_id = test_tau(Xs_2, phenotypes.flatten(), genotypes,  covariates, id_kernel, lr_pi,
                 temp_dir = 'temp_files/pi_0_gamma_0.txt',
                 out_dir = 'type1_sims/pi_0_gamma_0.txt',
                 regression = 'logistic')
        p_taus['id'].append(float(p_tau_id))
        p_taus['fixed'].append(float(p_pi))
        
        print(p_tau_dist, p_tau_id)
            
        #except Exception as e:
        #    print(e)
        #    continue
        
    #df = pd.DataFrame(list(zip(p_pis, p_taus_dist, p_taus_id)), columns=['p_pi', 'p_tau_dist', 'p_tau_id'])
    df = pd.DataFrame.from_dict(p_taus, orient='index').transpose()
    out_dir = 'power/v1_experiment/var_ratio_' + str(var_ratio) + '/opp_'+ str(domains) + 'domains_'
    out_dir += str(num_muts) + 'muts_' + str(num_indivs) + 'indiv.csv'
    #out_dir = 'type1_error/type1sims_500indiv_100muts_dist_kernel.csv'
    #df.to_csv(out_dir)
    
if __name__ == '__main__':
    var_ratios = [0.015, 0.02, 0.025, 0.03, 0.035]
    alphas = [0.5, 1, 2, 3, 4]
    ts = [3,6,9,12]
    indivs = [500]
    
    for vr in [0.035]:
        main(var_ratio=vr, pi_sum_ratio=0.0, num_muts=100, num_indivs=500)
                
            
            
        