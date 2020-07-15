#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:48:18 2020

@author: shuvomsadhuka
"""
import numpy as np
import pandas as pd
from scipy import stats

def type1(df, test, alpha):
    if test == 'fisher':
        threshold = stats.chi2.isf(q=alpha, df=4)
        type1 = (df['fisher'] > threshold).sum()/len(df['fisher'])
        return(type1)
    
    type1 = (df[test] < alpha).sum()/len(df['fisher'])
    return(type1)



def make_fisher(df):
    df['p_gamma'][df['p_gamma'] == 0] = 0.0000000001
    print(df['p_gamma'][61])
    df['fisher'] = -2 * np.log(df['p_gamma']) - 2 * np.log(df['p_pi'])
    return(df)


def make_tippett(df):
    df['tippett'] = 1 - (1 - np.square(df[['p_pi', 'p_gamma']].max(axis=1)))
    return(df)


def main():
    df = pd.read_csv('power/nonzero_gamma_only/100muts_500indiv_gamma_mixed_1.csv')
    df = make_fisher(df)
    df = make_tippett(df)
    
    type1_fisher = type1(df, 'fisher', 0.001)
    type1_tippett = type1(df, 'tippett', 0.001)
    type1_pi = type1(df, 'p_pi', 0.001)
    type1_gamma = type1(df, 'p_gamma', 0.001)
    print(type1_pi)
    print(type1_gamma)
    print(type1_fisher)
    print(type1_tippett)
    
    
if __name__ == '__main__':
    main()