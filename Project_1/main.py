### This file contains source code for running the alignments.
### mainly written by Qi Gao

from linear import Needleman_Wunsch_linear
from affine import Needleman_Wunsch_affine
from read_fasta_and_weight import read_fasta
from read_fasta_and_weight import read_sc
from read_fasta_and_weight import nw_score_genes
from read_fasta_and_weight import return_sc
from generate_output import generate_output


filepath = "/hpc/group/cscbb561s23/projects/project1/"

substitution_scores = read_sc('{}matrix.txt'.format(filepath))  ### This is written by Kehali Woldemichael

# Part 1
gene1 = read_fasta('{}close-first.fasta'.format(filepath))
gene2 = read_fasta('{}close-second.fasta'.format(filepath))
final_res = []

for k in range(len(gene1)):
    key1 = list(gene1.keys())
    key2 = list(gene2.keys())
    align_res = Needleman_Wunsch_linear(gene1[key1[k]], gene2[key2[k]], substitution_scores)
    final_res += generate_output(align_res, key1, key2, k, 0)

outfile = open('linear-close.txt', 'w')
outfile.writelines(final_res)
outfile.close()

gene1 = read_fasta('{}distant-first.fasta'.format(filepath))
gene2 = read_fasta('{}distant-second.fasta'.format(filepath))
final_res = []

for k in range(len(gene1)):
    key1 = list(gene1.keys())
    key2 = list(gene2.keys())
    align_res = Needleman_Wunsch_linear(gene1[key1[k]], gene2[key2[k]], substitution_scores)
    final_res += generate_output(align_res, key1, key2, k, 0)

outfile = open('linear-distant.txt', 'w')
outfile.writelines(final_res)
outfile.close()

# Part 2
gene1 = read_fasta('{}close-first.fasta'.format(filepath))
gene2 = read_fasta('{}close-second.fasta'.format(filepath))
final_res = []

for k in range(len(gene1)):
    key1 = list(gene1.keys())
    key2 = list(gene2.keys())
    align_res = Needleman_Wunsch_affine(gene1[key1[k]], gene2[key2[k]], substitution_scores)
    final_res += generate_output(align_res, key1, key2, k, 0)

outfile = open('affine-close.txt', 'w')
outfile.writelines(final_res)
outfile.close()

gene1 = read_fasta('{}distant-first.fasta'.format(filepath))
gene2 = read_fasta('{}distant-second.fasta'.format(filepath))
final_res = []

for k in range(len(gene1)):
    key1 = list(gene1.keys())
    key2 = list(gene2.keys())
    align_res = Needleman_Wunsch_affine(gene1[key1[k]], gene2[key2[k]], substitution_scores)
    final_res += generate_output(align_res, key1, key2, k, 0)

outfile = open('affine-distant.txt', 'w')
outfile.writelines(final_res)
outfile.close()