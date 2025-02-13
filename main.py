# Maulahna Robinson
# CS732 Bioinformatics

import numpy as np

def needleman_wunsch(seq1, seq2, match_score,mismatch_score,gap_pen):

# Initalize matrix and fill with gap penalties
    n,m = len(seq1), len(seq2)
    dp = np.zeros((n+1,m+1),dtype=int)

    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + gap_pen
    for j in range(1, m+1):
        dp[0][j] = dp[0][j-1] + gap_pen

# Fill the matrix completely
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] ==
                                    seq2[j-1] else mismatch_score)
            delete = dp[i-1][j] + gap_pen
            insert = dp[i][j-1] + gap_pen
            dp[i][j] = max(match,delete,insert)

# Find optimal alignment
    aligned_seq1,aligned_seq2 = "",""
    i,j = n,m

    while i > 0 or j < 0:
        current_score = dp[i][j]
        if i > 0 and j > 0 and (current_score == dp[i-1][j-1] + 
                            (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and (current_score == dp[i-1][j] + gap_pen):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1