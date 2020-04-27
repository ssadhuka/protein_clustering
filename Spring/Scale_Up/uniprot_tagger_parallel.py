#!/broad/software/dotkit/bin/dkenv --dotkits=Python-3.6 --command=python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 22:51:32 2020

@author: shuvomsadhuka
"""

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
import sys
import os
sys.path.append('../Initial_Tasks')
import fix_permutations as tests
import argparse


def run_test(prot_dir, protein, data):
        results_df = pd.DataFrame(columns=['gene_name', 'uniprot_ID', 'Mann-Whitney p-val', 'permutation risk', 'permutation prot'])
        file_repo = 'SWISS-MODEL_Repository/' + prot_dir + '/swissmodel/'
        if os.path.isdir(file_repo):
                try:
                        print(file_repo)
                        pdb_file = file_repo + str(os.listdir(file_repo)[0])
                                
                        ppdb = PandasPdb().read_pdb(pdb_file)
                        df = pd.DataFrame(ppdb.df['ATOM'])
                        sequence = ppdb.amino3to1()
                        
                        protein_spec_df = data[data['uniprot_repo'] == prot_dir]
                        gene_name = protein_spec_df['gene'].values[0]
                        uniprot_ID = protein_spec_df['uniprot'].values[0]
                        protein_spec_df = protein_spec_df[['mutation', 'effect_size', 'p-value', 'transition']]
                        
        
                        df_write = "protein_structs/" +  gene_name + '.csv'
                        df.to_csv(df_write, header=None, index=None, sep='\t')
        
                        write_to_dir = 'protein_mutation_locs_txts/' + gene_name + '.T2D.txt'
                        protein_spec_df.to_csv(write_to_dir, header=None, index=None, sep='\t')
                        protein_df = pd.read_csv(write_to_dir, header = None, sep="\t")
                        muts_df = tests.make_dataframe(df, protein_df, sequence) 
                        if not (muts_df[muts_df['score'] > 0].empty or muts_df[muts_df['score'] < 0].empty): 
                                risk = tests.get_dist_vec(muts_df, True)
                                prot = tests.get_dist_vec(muts_df, False)
                        
                                mw = tests.mannwhitneyu(risk, prot)
                                print(mw.pvalue)
                        
                                perm = tests.run_permutation(muts_df, df, np.mean(risk), np.mean(prot), 1000)
                                print(perm)
                        
                                new_row = {'gene_name': gene_name, 'uniprot_ID': uniprot_ID, 'Mann-Whitney p-val': mw.pvalue, "permutation risk": perm[0], "permutation prot": perm[1]}
                                results_df.append(new_row, ignore_index=True)
                        
                        out_csv = 'parallelized/' + str(protein) + '-pval.csv'
                        results_df.to_csv(out_csv)
        
                except Exception as e:
                        print(e)
                        pass
            
def main():
        parser = argparse.ArgumentParser(description='Get a uniprot ID.')
        parser.add_argument("--uniprot")

        args = vars(parser.parse_args())
        prot = args['uniprot']
        prot_dir = prot[0:2] + '/' + prot[2:4] + '/' + prot[4:6]
    
        #print(os.getcwd())
        #os.chdir('/home/unix/ssadhuka/protein_clusters/Scale_Up')
    
        data = pd.read_csv('all_betas_for_only_miss_vars.t2d.txt', header = None, sep="\t")
        lookup = pd.read_csv('all_genes_symbols_uniprots.txt', header = None, sep="\t")
    
        data.columns = ['mutation', 'effect_size', 'p-value', 'transition', 'gene']
        lookup.columns = ['gene', 'uniprot']
    
        data = pd.merge(data, lookup, on='gene')
    
        data['uniprot_repo'] = data['uniprot'].astype(str).str[0:2] + '/' + data['uniprot'].astype(str).str[2:4] + '/' + data['uniprot'].astype(str).str[4:6]
        # prot_dir = 'Q0/11/13'
        # prot = 'Q01113'
        run_test(prot_dir, prot, data)
        #proteins = set(data['uniprot_repo'])
    
    
if __name__ == '__main__':
        main()
