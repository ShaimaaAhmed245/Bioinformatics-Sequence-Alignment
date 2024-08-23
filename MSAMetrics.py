import numpy as np
from math import log

def sum_of_pairs(seq, gap, match, mismatch):
    matrix = np.zeros((len(seq),len(seq)))
    sop = 0
    align_len = len(seq[0].seq)
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            seq1=seq[i].seq 
            seq2=seq[j].seq
            score = getScore(seq1, seq2, gap, match, mismatch)
            matrix[i][j] = score
            sop = sop + score
        matrix[i][i] = align_len*match
    return sop, matrix

def getScore(xSeq, ySeq, gap, match, mismatch):
    score = 0
    for i in range (len(xSeq)):
        if xSeq[i] == '-' or ySeq[i] == '-':
            score = score + gap
        elif xSeq[i] == ySeq[i]:
            score = score + match
        elif xSeq[i] != ySeq[i]:
            score = score + mismatch
    return score

def percent_identity(seq):
    matrix = np.zeros((len(seq),len(seq)))

    total_pairs = len(seq[0].seq)
    # how many identical from the two sequence 
    identical_pairs=0
    percent_Identity = 0
    total_percent_Identity = 0

    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            seq1=seq[i].seq 
            seq2=seq[j].seq
            identical_pairs=0
            for k in range (len(seq1)):
                if (seq1[k] == seq2[k] and (seq1[k]!='-' or seq2[k]!='-' )):
                    identical_pairs += 1
            percent_Identity=identical_pairs/total_pairs*100
            total_percent_Identity +=  percent_Identity
            matrix[i][j] = round(percent_Identity,2)
        matrix[i][i] = 100
    return total_percent_Identity, matrix


def Mutual_Identity(seq):
    residue_frequency = dict()
    matrix = np.zeros((len(seq),len(seq)))
    # frequency of A C G T -
    for i in seq:
        for residue in i:
            if residue in residue_frequency:
                residue_frequency[residue] +=1
            else:
                residue_frequency[residue] = 1
    # frequency of two letters from two sequences ex: AA
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            seq1=seq[i].seq 
            seq2=seq[j].seq      
            for k in range(len(seq[0])):
                residue = seq1[k] + seq2[k]
                if residue in residue_frequency:
                    residue_frequency[residue] +=1
                else:
                    residue_frequency[residue] = 1

    total_residues =sum(residue_frequency.values()) 
    mutual_information = 0
    row_mutual_information = 0

    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            seq1=seq[i].seq 
            seq2=seq[j].seq
            #For percent identity analysis
            row_mutual_information = 0
            for k in range(len(seq1)):
                residue1,residue2=seq1[k], seq2[k]
                p_residue1_residue2 = residue_frequency[residue1+residue2] / total_residues
                p_residue1 = residue_frequency[residue1] / total_residues
                p_residue2 = residue_frequency[residue2] / total_residues
                mutual_information += log(p_residue1_residue2 /(p_residue1*p_residue2))
                row_mutual_information += log(p_residue1_residue2 /(p_residue1*p_residue2))
            matrix[i][j] = round(row_mutual_information,2)
    return mutual_information, matrix
