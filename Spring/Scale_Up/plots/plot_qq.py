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


measurements = pd.read_csv('all_2.csv').drop_duplicates(keep=False)
print(measurements)


#measurements['permutation risk'] = measurements['permutation prot']

measurements = measurements[~measurements['permutation risk'].isin(['EMPTY'])]
measurements = measurements.dropna(subset=['permutation risk'])
measurements['permutation risk'] = pd.to_numeric(measurements['permutation risk'])
measurements['permutation -log10'] = -np.log10(measurements['permutation risk']/100)


theoretical = [(i+1)/len(measurements['permutation -log10']) for i in range(len(measurements['permutation -log10']))]
theoretical = -np.log10(theoretical)
#stats.probplot(measurements['Mann-Whitney p-val']*2, dist=stats.uniform, plot=pylab)
#pylab.show()



measurements = measurements.sort_values(by=['permutation -log10'], ascending=False)
measurements['theoretical'] = theoretical

#measurements.to_csv('qq._lambda.csv')
p1= [0, 4]
p2=[0, 4]
plt.plot(p1,p2, c='r')
plt.scatter(measurements['theoretical'], measurements['permutation -log10'])
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Observed Quantiles')
plt.show()

