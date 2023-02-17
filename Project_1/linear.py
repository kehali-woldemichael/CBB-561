### This file contains source code for Needleman-Wunsch Algorithm with linear scoring
### written by Qi Gao

from read_fasta_and_weight import return_sc

def Needleman_Wunsch_linear(seq1, seq2, substitution_scores, gap_penalty = 1):
    l1 = len(seq1)
    l2 = len(seq2)
    
    # fill in
    F = [([0] * (l1 + 1)) for i in range(l2 + 1)]
    T = [([0] * (l1 + 1)) for i in range(l2 + 1)]
    for i in range(l2):
        for j in range(l1):
            cand = [F[i][j] + return_sc(substitution_scores, seq1[j], seq2[i]),
                F[i][j + 1] - gap_penalty, F[i + 1][j] - gap_penalty]
            idx = cand.index(max(cand))
            T[i + 1][j + 1] = idx
            F[i + 1][j + 1] = cand[idx]
    
    # trace back
    seq1_align = ""
    seq2_align = ""
    align_type = ""
    
    while (i + j >= 0):
        if (T[i + 1][j + 1] == 0):
            seq1_align += seq1[j]
            seq2_align += seq2[i]
            if (seq1[j] == seq2[i]):
                align_type += "|"
            else: 
                align_type += "*"
            i -= 1
            j -= 1
        elif (T[i + 1][j + 1] == 1):
            seq1_align += "-"
            seq2_align += seq2[i]
            align_type += " "
            i -= 1
        else:
            seq1_align += seq1[j]
            seq2_align += "-"
            align_type += " "
            j -= 1
        
    seq1_align = seq1_align[::-1]
    seq2_align = seq2_align[::-1]
    align_type = align_type[::-1]
    
    # calculate summary statistics
    nmatch = align_type.count("|")
    percent_identity = 2 * nmatch / (l1 + l2)
    align_temp = "|" + seq1_align + "|"
    align_temp1 = [i for i in align_temp.split("-") if i != ""]
    align_temp = "|" + seq2_align + "|"
    align_temp2 = [i for i in align_temp.split("-") if i != ""]
    nindels = len(align_temp1) + len(align_temp2) - 2
    mean_indel_length = align_type.count(" ") / nindels
    alignment_length = len(align_type)
    alignment_score = cand[idx]
    
    return seq1_align, seq2_align, align_type, nmatch, percent_identity, nindels, mean_indel_length, alignment_length, alignment_score