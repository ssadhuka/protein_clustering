#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:05:07 2020

@author: shuvomsadhuka
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 100)


def make_fisher_tippett_ablation(df, alpha):
    try:
        df['fixed'] = pd.to_numeric(df['fixed'])
    except:
        df['fixed'] = df['fixed'].apply(lambda st: st[st.find("(")+1:st.find("+")])
        df['fixed'] = pd.to_numeric(df['fixed'])
        
    fisher_threshold = stats.chi2.isf(q=alpha, df=4)
    fishers, tippetts, fixed, randoms = {}, {}, {}, {}
    for col in df.columns:
        df[col][df[col] == 0] = 0.0000000001
        if col != 'fixed':
            new_col_f = col + '_fisher'
            new_col_t = col + '_tippett'
            df[new_col_f]  = -2 * np.log(df[col]) - 2 * np.log(df['fixed'])
            df[new_col_t] = 1 - np.square(1 - df[['fixed', col]].min(axis=1))
            
            fishers[col] = (df[new_col_f] > fisher_threshold).sum()/len(df[new_col_f])
            tippetts[col] = (df[new_col_t] < alpha).sum()/len(df[new_col_t])
            randoms[col] = (df[col] < alpha).sum()/len(df[col])
    fixed['fixed'] = (df['fixed'] < alpha).sum()/len(df['fixed'])
    print(fishers, tippetts, randoms, fixed)
    return(fishers, tippetts, randoms, fixed)
    
            
def main():
    plt.figure(figsize=(8, 6))
    var = 0.025
    pi_ratios = [0]
    tests = ['fixed', 'random_dist', 'random_id','fisher_dist', 'fisher_id',
             'tippett_dist','tippett_id']
    colors = {'fixed': 'blue', 'random_dist': 'orange', 'random_id': 'darkorange',
              'fisher_dist': 'green', 'tippett_dist': 'red', 
              'fisher_id': 'darkgreen', 'tippett_id': 'darkred'}
    vers = ['v1_experiment']
    alpha = 0.01
    domains = 1
    power = 2
    pow_key = str(power)
    
    ts = [3, 6, 9, 12]
    p_fishers, p_tippetts, p_randoms, p_fixed = [], [], [], []
    for v in vers:
        #read_dir = '../../type1_error/type1sims_500indiv_100muts.csv'
        read_dir = 'power/' + str(v) + '/var_ratio_' + str(var)
        read_dir += '/1domains_100muts_' + str(500) + 'indiv.csv'
        df_in = pd.read_csv(read_dir, index_col=0)
        fishers, tippetts, randoms, fixed = make_fisher_tippett_ablation(df_in, alpha)
    
        for key in fishers.keys():
            if key.startswith(pow_key) or key.endswith('id'):
                if key.startswith(pow_key):
                    print(key)
                    p_fishers.append(fishers[key])
                    p_tippetts.append(tippetts[key])
                    p_randoms.append(randoms[key])
                    
                    
        plt.plot(ts, p_fishers, c='green', label = 'fisher_dist', marker='x')
        plt.plot(ts, p_tippetts, c='red', label = 'tippett_dist', marker='x')
        plt.plot(ts, p_randoms, c='cyan', label = 'random_dist', marker='x')
        plt.plot(ts, [float(fishers['id']) for i in range(4)], c='darkgreen', marker='o',
                 label = 'fisher_id')
        plt.plot(ts, [float(tippetts['id']) for i in range(4)], c='darkred', marker='o',
                 label = 'tippett_id')
        plt.plot(ts, [float(randoms['id']) for i in range(4)], c='darkcyan', marker='o',
                 label = 'random_id')
        plt.plot(ts, [float(fixed['fixed']) for i in range(4)], c='black', label = 'fixed')
            
    
    plt.xlabel('t constant')
    plt.ylabel('Power')
    title = 'Power Simulations, alpha = ' + str(alpha) + ', domains = ' + str(domains) + ', exponent = ' + str(power)
    title += ', Explained Variance = '+ str(var)
    plt.title(title)
    plt.legend()
    
    save_to = 'plots/power_v1_experiment_alpha' + str(alpha) + '_domains' + str(domains) + '_exponent' + str(power)
    save_to += '_var' + str(var) + '.png'
    plt.savefig(save_to, dpi=400)
    
if __name__ == '__main__':
    main()