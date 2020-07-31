#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 23:35:32 2020

@author: shuvomsadhuka
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def rejection_region(df, alpha):
    threshold = stats.chi2.isf(q=alpha, df=4)
    reject_fisher_dist = (df['fisher_dist'] > threshold).sum()/len(df['fisher_dist'])
    reject_fisher_id = (df['fisher_id'] > threshold).sum()/len(df['fisher_dist'])
    
    reject_pi = (df['p_pi'] < alpha).sum()/len(df['p_pi'])
    reject_tau_dist = (df['p_tau_dist'] < alpha).sum()/len(df['fisher_dist'])
    reject_tippett_dist = (df['tippett_dist'] < alpha).sum()/len(df['fisher_dist'])
    
    reject_tau_id = (df['p_tau_id'] < alpha).sum()/len(df['fisher_id'])
    reject_tippett_id = (df['tippett_id'] < alpha).sum()/len(df['fisher_id'])
    return([reject_pi, reject_tau_dist, reject_tau_id, reject_fisher_dist,
            reject_fisher_id, reject_tippett_dist, reject_tippett_id])



def make_fisher(df):
    df['p_tau_dist'][df['p_tau_dist'] == 0] = 0.0000000001
    df['p_tau_id'][df['p_tau_id'] == 0] = 0.0000000001
    try:
        df['p_pi'] = pd.to_numeric(df['p_pi'])
    except:
        df['p_pi'] = df['p_pi'].apply(lambda st: st[st.find("(")+1:st.find("+")])
        df['p_pi'] = pd.to_numeric(df['p_pi'])
    #print(df)
    
    df['fisher_dist'] = -2 * np.log(df['p_tau_dist']) - 2 * np.log(df['p_pi'])
    df['fisher_id'] = -2 * np.log(df['p_tau_id']) - 2 * np.log(df['p_pi'])
    return(df)


def make_tippett(df):
    print(1-np.square(df[['p_pi', 'p_tau_dist']].min(axis=1)))
    df['tippett_dist'] = 1 - np.square(1 - df[['p_pi', 'p_tau_dist']].min(axis=1))
    df['tippett_id'] = 1 - np.square(1 - df[['p_pi', 'p_tau_id']].min(axis=1))
    return(df)


def main():
    plt.figure(figsize=(8, 6))
    expld_var = [0.02, 0.025]
    pi_ratios = [0]
    tests = ['fixed', 'random_dist', 'random_id','fisher_dist', 'fisher_id',
             'tippett_dist','tippett_id']
    colors = {'fixed': 'blue', 'random_dist': 'orange', 'random_id': 'darkorange',
              'fisher_dist': 'green', 'tippett_dist': 'red', 
              'fisher_id': 'darkgreen', 'tippett_id': 'darkred'}
    vers = ['v1_experiment']
    alpha = 0.05
    
    for v in vers:
        for i, var in enumerate(expld_var):
            #read_dir = '../../type1_error/type1sims_500indiv_100muts.csv'
            read_dir = 'power/' + str(v) + '/var_ratio_' + str(var)
            read_dir += '/1domains_100muts_' + str(500) + 'indiv.csv'
            df_in = pd.read_csv(read_dir, index_col=0)
            
            rejections = rejection_region(df_in, alpha)
            if i == 0:
                power = rejections
            else:
                power = np.vstack((power, rejections))
        for i, test in enumerate(tests):
            lbl = str(test)
            if i == 0:
                plt.plot(expld_var, power[:,i], marker='o', label = lbl, color=colors[test])
            if i > 0 and i % 2 == 0:
                plt.plot(expld_var, power[:,i], marker='o', label = lbl, color=colors[test])
            elif i > 0 and i % 2 == 1:
                plt.plot(expld_var, power[:,i], marker='x', label = lbl, color=colors[test])
            
    
    plt.xlabel('Total Variance Explained')
    plt.ylabel('Power')
    title = 'Power Simulations, alpha = ' + str(alpha) + ' n_subjects = ' + str(500)
    plt.title(title)
    plt.legend()
    
    save_to = 'power_v1_initial_500_0-05_redo.png'
    #plt.savefig(save_to, dpi=400)
        
        
        
        
if __name__ == '__main__':
    main()