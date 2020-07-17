import sys
from optparse import OptionParser

import numpy as np
import math
import scipy.stats
import gzip
import subprocess
import math

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("","--ld-file",default=None)
parser.add_option("","--ld-threshold",default=0,type="float")
parser.add_option("","--N",default=100000,type="float") #sample size for sumstat simulations
parser.add_option("","--beta-file",default=None) #specify mapping from SNP to true beta
parser.add_option("","--beta-dist-file",default=None) #specify a distribution of beta values to be assigned to each SNP within a region. No header, columns are chrom, begin, end, p, mean, sd [replicate] where if the optional [1-based] replicate column is present, it will apply the region only to the specified replicate. If multiple lines apply to a SNP, it will take the first one
parser.add_option("","--max-component-size",type="int",default=np.inf) #control the maximum number of SNPs to be included in an LD-block (smaller is faster)
parser.add_option("","--assume-r",action="store_true")
parser.add_option("","--num-sim",type="int",default=1) #generate statistics for multiple indpendent replicates
parser.add_option("","--ldsc-format",action="store_true") #output sumstats in LDSC format
parser.add_option("","--num-causal-snps-out",default=None) #write a file with the number of causal SNPs per iteration

(options, args) = parser.parse_args()

def bail(message):
    sys.stderr.write("%s\n" % message)
    sys.exit(1)
def warn(message):
    sys.stderr.write("Warning: %s\n" % message)
def log(message):
    sys.stderr.write("%s\n" % message)
def open_gz(file, sort=None):
    if sort is None:
        if file[-3:] == ".gz":
            fh = gzip.open(file)
            fh.readline()
        else:
            fh = open(file)
            fh.readline()
    else:
        if file[-3:] == ".gz":
            return subprocess.Popen(["zcat %s | tail -n+2 | sort -grk%s" % (file, sort)], stdout=subprocess.PIPE, shell=True).stdout
        else:
            return subprocess.Popen(["tail -n+2 %s | sort -grk%s" % (file, sort)], stdout=subprocess.PIPE, shell=True).stdout
        
assume_r = options.assume_r

ld_file = options.ld_file

if ld_file is None:
    bail("Need --ld-file")

ld_fh = open_gz(ld_file, 7)
snp_to_component = {}
component_to_snp = {}
index = 0
component = 0
max_component_size = 2
chr_to_snp_pos = {}
snp_to_chr_pos = {}
chr_pos_to_snp = {}

for line in ld_fh:
    cols = line.strip().split()
    if len(cols) != 7:
        bail("Bad line")
    value = float(cols[6])
    if value < 0:
        assume_r = True

    if abs(value) < options.ld_threshold:
        continue
    snp_1 = cols[2]
    snp_2 = cols[5]

    snp_1_chr = cols[0]
    snp_1_pos = int(cols[1])
    snp_to_chr_pos[snp_1] = (snp_1_chr, snp_1_pos)
    snp_2_chr = cols[3]
    snp_2_pos = int(cols[4])
    snp_to_chr_pos[snp_2] = (snp_2_chr, snp_2_pos)

    if snp_1_chr not in chr_to_snp_pos:
        chr_to_snp_pos[snp_1_chr] = set()
    chr_to_snp_pos[snp_1_chr].add(snp_1_pos)
    chr_pos_to_snp[(snp_1_chr, snp_1_pos)] = snp_1
    if snp_2_chr not in chr_to_snp_pos:
        chr_to_snp_pos[snp_2_chr] = set()
    chr_to_snp_pos[snp_2_chr].add(snp_2_pos)
    chr_pos_to_snp[(snp_2_chr, snp_2_pos)] = snp_2

    if snp_1 not in snp_to_component and snp_2 not in snp_to_component:
        component += 1
        snp_to_component[snp_1] = component
        snp_to_component[snp_2] = component
        component_to_snp[component] = set()
        component_to_snp[component].add(snp_1)
        component_to_snp[component].add(snp_2)
    elif snp_1 in snp_to_component and snp_2 not in snp_to_component:
        if len(component_to_snp[snp_to_component[snp_1]]) < options.max_component_size:
            snp_to_component[snp_2] = snp_to_component[snp_1]
            component_to_snp[snp_to_component[snp_1]].add(snp_2)
        else:
            component += 1
            snp_to_component[snp_2] = component
            component_to_snp[component] = set()
            component_to_snp[component].add(snp_2)
    elif snp_2 in snp_to_component and snp_1 not in snp_to_component:
        if len(component_to_snp[snp_to_component[snp_2]]) < options.max_component_size:
            snp_to_component[snp_1] = snp_to_component[snp_2]
            component_to_snp[snp_to_component[snp_2]].add(snp_1)
        else:
            component += 1
            snp_to_component[snp_1] = component
            component_to_snp[component] = set()
            component_to_snp[component].add(snp_1)
    elif snp_2 in snp_to_component and snp_1 in snp_to_component and snp_to_component[snp_1] != snp_to_component[snp_2]:
        if len(component_to_snp[snp_to_component[snp_1]]) + len(component_to_snp[snp_to_component[snp_2]]) <= options.max_component_size:
            component_1 = snp_to_component[snp_1]
            component_2 = snp_to_component[snp_2]
            for snp in component_to_snp[component_2]:
                snp_to_component[snp] = component_1
            component_to_snp[component_1] = component_to_snp[component_1].union(component_to_snp[component_2])
            component_to_snp.pop(component_2)

    if len(component_to_snp[snp_to_component[snp_1]]) > max_component_size:
        max_component_size = len(component_to_snp[snp_to_component[snp_1]])
log("Max component size: %s" % max_component_size)

ld_fh.close()

for chrom in chr_to_snp_pos:
    chr_to_snp_pos[chrom] = sorted(list(chr_to_snp_pos[chrom]))

snp_to_index = {}
component_to_cor = {}
for component in component_to_snp:
    index = 0
    #print "%s %s" % (component, component_to_snp[component])
    for snp in sorted(component_to_snp[component]):
        snp_to_index[snp] = index
        index += 1
    #print("%s %s" % (component, len(component_to_snp[component])))
    component_to_cor[component] = np.identity(len(component_to_snp[component]), dtype=np.float64)

ld_fh = open_gz(ld_file, 7)
for line in ld_fh:
    cols = line.strip().split()
    value = float(cols[6])
    if abs(value) < options.ld_threshold:
        continue
    if not assume_r:
        value = math.sqrt(value)
    snp_1 = cols[2]
    snp_2 = cols[5]
    if snp_to_component[snp_1] == snp_to_component[snp_2]:
        component_to_cor[snp_to_component[snp_1]][snp_to_index[snp_1], snp_to_index[snp_2]] = value
        component_to_cor[snp_to_component[snp_1]][snp_to_index[snp_2], snp_to_index[snp_1]] = value
ld_fh.close()

snp_to_true_beta_params = {}
if options.beta_dist_file:
    import random
    import bisect
    beta_dist_fh = open(options.beta_dist_file)
    for line in beta_dist_fh:
        cols = line.strip().split()
        if len(cols) != 6 and len(cols) != 7:
            warn("Ignoring line without six columns for (chrom, start, end, p, mean, var [replicate]):  %s" % line.strip())
            continue
        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        p = float(cols[3])
        if p < 0 or p > 1:
            log("Error: bad value for bernoulli (%s)" % (p))
            continue
        mean = float(cols[4])
        var = float(cols[5])
        if var < 0:
            log("Error: bad value for var (%s)" % (var))
            continue

        if len(cols) == 7:
            rep = int(cols[6])
            if rep < 1:
                log("Error: bad value for rep (%s)" % (rep))
                continue
        else:
            rep = None

        #find all overlapping snps

        if chrom in chr_to_snp_pos:
            start_ind = bisect.bisect_left(chr_to_snp_pos[chrom], start)
            for ind in range(start_ind, len(chr_to_snp_pos[chrom])):
                if chr_to_snp_pos[chrom][ind] < start:
                    log("There is a bug -- %s is less than %s" % (chr_to_snp_pos[chrom][ind], start))
                    continue
                if chr_to_snp_pos[chrom][ind] > end:
                    break
                cur_snp = chr_pos_to_snp[(chrom, chr_to_snp_pos[chrom][ind])]
                if cur_snp not in snp_to_true_beta_params:
                    snp_to_true_beta_params[cur_snp] = []
                snp_to_true_beta_params[cur_snp].append((p, mean, var, rep-1))

snp_to_true_beta = {}
if options.beta_file:
    beta_fh = open(options.beta_file)
    for line in beta_fh:
        cols = line.strip().split()
        if len(cols) != 2:
            warn("Ignoring line without two columns for (snp, beta):  %s" % line.strip())
            continue
        snp = cols[0]
        if snp in snp_to_index:
            true_beta = cols[1]
            snp_to_true_beta[snp] = true_beta
    beta_fh.close()

snp_to_beta = {}
snp_to_p = {}

if options.ldsc_format:
    sys.stdout.write("SNP\tReplicate\tA1\tA2\tZ\tN\n")
else:
    sys.stdout.write("SNP\tChrom\tPos\tReplicate\tEffect\tStdErr\tP-value\n")
#file1 = open("/humgen/diabetes/UKBB_app27892/UKB_bim_copy/range_files/myfile.txt","w")
#file1.write("SNP\tEffect\tP-value\n")

num_true_causal_snps = {}

for component in component_to_cor:

    cor_matrix = component_to_cor[component]
    try:
        L = np.linalg.cholesky(cor_matrix)
    except np.linalg.linalg.LinAlgError:
        l, Q = np.linalg.eigh(cor_matrix)
        l[l < 0] = 0
        L = np.dot(Q, np.diag(np.sqrt(l)))
        #renormalize diagonal to be 1
        D2 = np.diag(1/np.sqrt(np.diag(np.dot(L, L.transpose()))))
        L = np.dot(D2, L)

    #now do the random sampling
    for it in range(0,options.num_sim):
        if it not in num_true_causal_snps:
            num_true_causal_snps[it] = 0
        cur_true_beta = np.zeros(len(component_to_snp[component]))
        if len(snp_to_true_beta_params) > 0:
            for snp in component_to_snp[component]:
                if snp in snp_to_true_beta_params:
                    for true_beta_params in snp_to_true_beta_params[snp]:
                        if true_beta_params[3] is None or true_beta_params[3] == it:
                            if random.random() < true_beta_params[0]:
                                cur_true_beta[snp_to_index[snp]] = np.random.normal(loc=true_beta_params[1], scale=np.sqrt(true_beta_params[2]), size=1)[0]
                            break

        if len(snp_to_true_beta) > 0:
            for snp in component_to_snp[component]:
                if snp in snp_to_true_beta:
                    cur_true_beta[snp_to_index[snp]] = snp_to_true_beta[snp]

        num_true_causal_snps[it] += np.sum(cur_true_beta != 0)

        u = np.random.normal(loc=0, scale=1, size=len(component_to_snp[component]))
        y = np.array(np.matmul(cor_matrix, cur_true_beta) * np.sqrt(options.N) + np.matmul(L, u))[0:]

        for snp in component_to_snp[component]:
            if len(component_to_snp[component]) == 1:
                beta = y
            else:
                beta = y[snp_to_index[snp]]
            pvalue = 2 * scipy.stats.norm.cdf(-1 * abs(beta))
            beta /= np.sqrt(options.N)
            se = 1 / np.sqrt(options.N)
            if options.ldsc_format:
                sys.stdout.write("%s\t%d\t%s\t%s\t%.3g\t%d\n" % (snp, it+1, "R", "A", beta / se, options.N))
            else:
                sys.stdout.write("%s\t%s\t%d\t%d\t%.3g\t%.3g\t%.3g\n" % (snp, snp_to_chr_pos[snp][0], snp_to_chr_pos[snp][1], it+1, beta, se, pvalue))

            #file1.write("%s\t%.3g\t%.3g\n" % (snp, beta, pvalue))

    #print("Component %s:" % component)
    #print("y:")
    #print(y)

if options.num_causal_snps_out is not None:
    out_fh = open(options.num_causal_snps_out, 'w')
    out_fh.write("Replicate\tNum_Causal\n")
    for it in range(0,options.num_sim):
        out_fh.write("%d\t%d\n" % (it+1, num_true_causal_snps[it]))
    out_fh.close()
