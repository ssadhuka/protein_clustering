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
    reject_fisher = (df['fisher'] > threshold).sum()/len(df['fisher'])
    
    reject_pi = (df['p_pi'] < alpha).sum()/len(df['fisher'])
    reject_gamma = (df['p_gamma'] < alpha).sum()/len(df['fisher'])
    reject_tippett = (df['tippett'] < alpha).sum()/len(df['fisher'])
    return([reject_pi, reject_gamma, reject_fisher, reject_tippett])


def make_fisher(df):
    df['p_gamma'][df['p_gamma'] == 0] = 0.0000000001
    try:
        df['p_pi'] = pd.to_numeric(df['p_pi'])
    except:
        df['p_pi'] = df['p_pi'].apply(lambda st: st[st.find("(")+1:st.find("+")])
        df['p_pi'] = pd.to_numeric(df['p_pi'])
    print(df)
    
    df['fisher'] = -2 * np.log(df['p_gamma']) - 2 * np.log(df['p_pi'])
    return(df)


def make_tippett(df):
    df['tippett'] = 1 - (1 - np.square(df[['p_pi', 'p_gamma']].max(axis=1)))
    return(df)


def main():
    plt.figure(figsize=(8, 6))
    expld_var = 0.02
    pi_ratios = [0.0, 0.25, 0.5, 0.75, 1.0]
    tests = ['fixed', 'random', 'fisher', 'tippett']
    colors = {'fixed': 'blue', 'random': 'orange', 'fisher': 'green', 'tippett': 'red'}
    vers = ['v0.5_debug_vars']
    
    for v in vers:
        for i, pr in enumerate(pi_ratios):
            #read_dir = '../../type1_error/type1sims_500indiv_100muts.csv'
            read_dir = 'power/' + str(v) + '/var_ratio_0.02_pi_ratio_' + str(pr)
            read_dir += '/100muts_' + str(500) + 'indiv.csv'
            df_in = pd.read_csv(read_dir, index_col=0)
            df_in = make_fisher(df_in)
            df_in = make_tippett(df_in)
            
            rejections = rejection_region(df_in, 0.00000001)
            if i == 0:
                power = rejections
            else:
                power = np.vstack((power, rejections))
           # print(pr, rejections)
        for i, test in enumerate(tests):
            legend_lbl = test + ' ' + str(v)
            if v == 'v0.5':
                plt.plot(pi_ratios, power[:,i], marker='o', label=legend_lbl, color=colors[test])
            if v == 'v0.5_identity_kernel':
                plt.plot(pi_ratios, power[:,i], marker='x', label=legend_lbl, color=colors[test])
            if v == 'v0.5_debug_vars':
                plt.plot(pi_ratios, power[:,i], marker='x', label=legend_lbl, color=colors[test])
    
    plt.xlabel('Percent Explained by Fixed')
    plt.ylabel('Power')
    title = 'Power Simulations, Total Explained Variance = ' + str(expld_var)
    plt.title(title)
    plt.legend()
    
    save_to = str(expld_var) + 'power_sims_debug_alpha1e-8.png'
    plt.savefig(save_to, dpi=400)
        
        
        
        
if __name__ == '__main__':
    main()