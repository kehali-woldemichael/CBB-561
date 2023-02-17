### This file contains source code for reading sequences from fasta files and weight matrix from matrix.txt file.
### written by Kehali Woldemichael

import sys, re

def read_fasta(filename):
    fasta_file = open(filename,'r')
    lines = fasta_file.readlines()

    hre=re.compile('>(\S+)')
    lre=re.compile('^(\S+)$')

    gene={}

    for line in lines:
        outh = hre.search(line)
        #print(line)
        if outh:
            id = outh.group(1)
            print(outh.group(1))
        else:
            outl = lre.search(line)
            if(id in gene.keys()):
                gene[id] += outl.group(1)
            else:
                gene[id] = outl.group(1)
    return gene

def read_sc(sc_file):
    substitution_scores = open(sc_file)
    lines = substitution_scores.readlines()[1:]
    substitution_scores = {}
    for line in lines:
        if line[0] != '\n':
            id = line[0] 
            substitution_scores[id] = [int(x) for x in line[1:].replace("\n", "").split()]
    return substitution_scores

def nw_score_genes(gene_file_1, gene_file_2, sc_file, scoring_type):
    gene_1 = read_fasta(gene_file_1)
    gene_2 = read_fasta(gene_file_2)
    substitution_scores = read_sc(sc_file)

def return_sc(substitution_scores, n_1, n_2):
    if n_1 == 'A': return substitution_scores[n_2][0]
    if n_1 == 'G': return substitution_scores[n_2][1]
    if n_1 == 'C': return substitution_scores[n_2][2]
    if n_1 == 'T': return substitution_scores[n_2][3]
    if n_1 == 'N': return substitution_scores[n_2][4]