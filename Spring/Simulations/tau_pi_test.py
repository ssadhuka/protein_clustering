#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:32:07 2020

@author: shuvomsadhuka
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from scipy.linalg import sqrtm
from scipy.stats import chi2
import os
#import rpy2
#from rpy2.robjects.packages import importr
#stats = importr('stats')
#base = importr('base')
#compquadform = importr('CompQuadForm')

    
            
def test_pi(X, y, genotypes, covariates, lr_alpha, regression):
    if regression == 'logistic':
        mu_hat = lr_alpha.predict_proba(X.T)[:,1]
        sigma_hat = mu_hat * (1 - mu_hat)
        D_hat = np.diag(sigma_hat)
    if regression == 'linear':
        mu_hat = lr_alpha.predict(X.T)
        sigma_hat = np.sum(np.square(y.flatten() - mu_hat.flatten()))/len(y.flatten())
        D_hat = np.diag([sigma_hat for i in range(len(y.flatten()))])
        
    U_pi = genotypes @ (y.flatten() - mu_hat.flatten())
    
    inner = (D_hat @ covariates.T) @ np.linalg.pinv(covariates @ D_hat @ covariates.T) @ (covariates @ D_hat)
    inner = D_hat - inner
    Sigma = genotypes @ inner @ genotypes.T
    U_pi_standardized = np.linalg.pinv(sqrtm(Sigma)) @ U_pi
    U_pi_chi = np.sum(np.square(U_pi_standardized))
    
    return(1-chi2.cdf(U_pi_chi, len(U_pi_standardized)))
    

def test_tau(X, y, genotypes, covariates, distance_kernel, lr_pi, temp_dir, out_dir, regression):
    if regression == 'logistic':
        mu_hat = lr_pi.predict_proba(X.T)[:,1]
        sigma_hat = mu_hat * (1 - mu_hat)
        D_hat = np.diag(sigma_hat)
    if regression == 'linear':
        mu_hat = lr_pi.predict(X.T)
        sigma_hat = np.sum(np.square(y.flatten() - mu_hat.flatten()))/len(y.flatten())
        D_hat = np.diag([sigma_hat for i in range(len(y.flatten()))])
    
    kernel = genotypes.T @ distance_kernel @ genotypes
    S_tau_2 = (y.flatten() - mu_hat.flatten()).T @ kernel @ (y.flatten() - mu_hat.flatten())
    
    def check_symmetric(a, tol=1e-3):
        return np.all(np.abs(a-a.T) < tol)
    
    inner = np.linalg.pinv(X @ D_hat @ X.T)
    P_hat = D_hat - (D_hat @ X.T @ inner @ X @ D_hat)
    
    P_hat_sqrt = sqrtm(P_hat)
    # pay attention to this line - do we need dist kernel in between?
    mixture = P_hat_sqrt @ (genotypes.T @ distance_kernel @ genotypes) @ P_hat_sqrt
    #print(check_symmetric(mixture))
    
    # mixture should definitely be PSD because it can be factored as a matrix M * M.T
    #mixture_with_tolerance = mixture + 

    '''
    try:
        L = np.linalg.cholesky(mixture)
    except np.linalg.linalg.LinAlgError:
        l, Q = np.linalg.eigh(mixture)
        l[l < 0] = 0
        L = np.dot(Q, np.diag(np.sqrt(l)))
        #print(np.dot(L, L.transpose()))
        D2 = np.diag(1/np.sqrt(np.diag(np.dot(L, L.transpose()))))
        L = np.dot(D2, L)
    '''
    w,v = np.linalg.eigh(mixture)
    w = np.sort(w)[::-1]
    eig_ratio = w[1]/w[0]
    i = 1
    eigs_davies = [w[0]]
    
    while (eig_ratio > 1e-5) and (i < len(w)-2) and (w[i] > 0):
        eig_ratio = w[i+1]/w[i]
        eigs_davies.append(w[i])
        i += 1
    return(davies_approximation(S_tau_2, eigs_davies, temp_dir, out_dir))

    

def fisher_procedure():
    pass


def tippett_procedure():
    pass


def davies_approximation(S_tau_2, eigs, temp_dir, out_dir):
    # for now, we only consider eigenvalues with multiplicity 1 as in Sun paper
    temp_relative = 'for_davies_alg/' + temp_dir
    out_relative = 'for_davies_alg/' + out_dir
    with open(temp_relative, 'w') as f:
        f.write("%s\n" % S_tau_2)
        for eig in eigs:
            f.write("%s\n" % eig)
    f.close()
    
    #os.environ['R_HOME'] = ';ibrary/Frameworks/R.framework/Resources'
    Rscript_call = '/usr/local/bin/Rscript for_davies_alg/davies_approx.R ' + temp_dir + ' ' +  out_dir
    os.system(Rscript_call)
    
    f = open(out_relative, 'r')
    p_gamma = f.readline()
    return(p_gamma)
    


def liu_approximation():
    pass