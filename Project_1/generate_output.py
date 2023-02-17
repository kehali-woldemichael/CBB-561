### This file contains source code for creating formatted alignment output.
### written by Qi Gao

def generate_output(align_res, seq1_name, seq2_name, i, score_dec):
    output = []
    output.append("Alignment #{}:\n".format(i + 1))
    output.append("\n")
    output.append("Sequence #1: {}\n".format(seq1_name[i]))
    output.append("Sequence #2: {}\n".format(seq2_name[i]))
    output.append("Matches: {}\n".format(align_res[3]))
    output.append("Percent identity: {}%\n".format(round(100 * align_res[4])))
    output.append("Indels: number = {}, mean length = {}\n".format(align_res[5], round(align_res[6], 1)))
    output.append("Alignment length: {}\n".format(align_res[7]))
    output.append("Score = {}\n".format(round(align_res[8], score_dec)))
    output.append("\n")
    for i in range(0, len(align_res[0]), 60):
        output.append("{}\n".format(align_res[0][i:i+60]))
        output.append("{}\n".format(align_res[2][i:i+60]))
        output.append("{}\n".format(align_res[1][i:i+60]))
        output.append("\n")
    return output
