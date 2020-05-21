#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 13:45:30 2020

@author: shuvomsadhuka
"""

import pandas as pd
import numpy as np 
import pylab 
import scipy.stats as stats
import matplotlib.pyplot as plt


measurements = pd.read_csv('results_5_19_20_3.csv')#.drop_duplicates(keep=False)

#measurements = measurements.drop_duplicates(subset=['gene_name'])

print(len(set(measurements['gene_name'])))

#subset = pd.read_csv('all_2.csv').drop_duplicates(keep=False)
#print(len(set(list(measurements['gene_name']))))



#pd.set_option('display.max_columns', None)
#measurements['permutation prot'] = measurements['permutation risk']

measurements['permutation risk'] = measurements['permutation prot']

measurements = measurements[~measurements['permutation risk'].isin(['EMPTY'])]
measurements = measurements.dropna()
measurements['permutation risk'] = pd.to_numeric(measurements['permutation risk'])
measurements['permutation -log10'] = -np.log10(measurements['permutation risk']/100)
#measurements = measurements.drop_duplicates(subset=['gene_name'])


theoretical = [(i+1)/len(measurements['permutation -log10']) for i in range(len(measurements['permutation -log10']))]
theoretical = -np.log10(theoretical)
#stats.probplot(measurements['Mann-Whitney p-val']*2, dist=stats.uniform, plot=pylab)
#pylab.show()



measurements = measurements.sort_values(by=['permutation -log10'], ascending=False)
measurements['theoretical'] = theoretical

#measurements = measurements[measurements['gene_name'].isin(list(subset['gene_name']))]
#print(measurements)

#measurements = measurements[measurements['Unnamed: 0'] == 0.0]



#measurements.to_csv('qq_bad1.csv')
p1= [0, 5]
p2=[0, 5]
plt.plot(p1,p2, c='r')
plt.scatter(measurements['theoretical'], measurements['permutation -log10'])
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Observed Quantiles')
plt.show()

print(list(measurements.head(20)['gene_name']))

