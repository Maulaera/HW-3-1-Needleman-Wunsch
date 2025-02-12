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
