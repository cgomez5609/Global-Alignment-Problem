"""
Author: Chris E.
Find a Highest-Scoring Alignment of Two Strings
Rosalind problem with the same name as above
"""

# Global Alignment Problem
import numpy as np
import pandas as pd

from Bio.SubsMat import MatrixInfo 

v = 'PLEASANTLY'
w = 'MEANLY'

s_matrix = MatrixInfo.pam250
penalty = 5

def build_matrix(v, w, penalty, score_matrix):
    """
    v and w are two strings: v as our rows and w as our columns.
    Creates our matrix and fills in scores based on the gap penalty named penalty
    and Pam 250 scoring matrix named score_matrix.
    This function also creates the directions for our matrix. If backtracking from
    source the letters represent a specific direction.
    u: up
    l: left
    d: diagonal
    s: source
    """
    arr = np.zeros(dtype=(int), shape=(len(v) +1, len(w) +1))
    arrDir = np.zeros(dtype=(str), shape=(len(v) +1, len(w) +1))
    
    # First columns
    for row in range(0, 1):
        for col in range(1, len(arrDir[row])):
            arr[row][col] = arr[row][col-1] + -penalty
            arrDir[row][col] = 'l'
            
    # First rows
    for row in range(1, len(arrDir)):
        for col in range(0, 1):
            arr[row][col] = arr[row-1][col] + -penalty
            arrDir[row][col] = 'u'
    
    
    arrDir[0][0] = 's'


    for i in range(1, len(arr)):
        for j in range(1, len(arr[i])):
            arr[i][j] = max(arr[i-1][j] - penalty, arr[i][j-1] - penalty) 
            if (w[j - 1], v[i - 1]) in score_matrix:
                arr[i][j] = max(arr[i-1][j-1] + score_matrix[w[j - 1], v[i - 1]], arr[i][j])
            else:
                arr[i][j] = max(arr[i-1][j-1] + score_matrix[v[i - 1], w[j - 1]], arr[i][j])
            
            if arr[i][j] == arr[i-1][j] - penalty:
                arrDir[i][j] = 'u'
            elif arr[i][j] == arr[i][j-1] - penalty:
                arrDir[i][j] = 'l'
            else:
                arrDir[i][j] = 'd'
                
    return arr, arrDir


def backtrack(arrDir):
    """
    Backtracking function to achieve alignment between two strings. Begins
    at the sink and works its way up to the source. Along the way appends 
    letter or gap to one of two lists (va or wa) depending on the directionality 
    of the matrix containing the directions. 
    va: list for v string
    wa: list for w string
    Lists then get reveresed and joined at the end.
    """
    va = []
    wa =[]
    index_v = len(v)-1
    index_w = len(w)-1
    row = len(v)
    col = len(w)
    sink = arrDir[row][col] 
    
    while sink != 's':
        if sink == 'd':
            va.append(v[index_v])
            wa.append(w[index_w])
            index_v -= 1
            index_w -= 1
            row -= 1
            col -= 1
            sink = arrDir[row][col]
        elif sink == 'u':
            wa.append('-')
            va.append(v[index_v])
            index_v -= 1
            row -= 1
            sink = arrDir[row][col]
        elif sink == 'l':
            va.append('-')
            wa.append(w[index_w])
            index_w -= 1
            col -= 1
            sink = arrDir[row][col]

    va = ''.join(reversed(va))
    wa = ''.join(reversed(wa))
    
    return va, wa
    
def main():
    matrix, direc = build_matrix(v, w, penalty, s_matrix)

    # Used pandas dataframe to help viusalize strings in our matrix
    # Not needed for actual problem
    dt_arr = pd.DataFrame(data=matrix, index=[0] + list(v), columns=[0] + list(w))
    dt_dir = pd.DataFrame(data=direc, index=[0] + list(v), columns=[0] + list(w))
    
    va, wa = backtrack(direc)
    score = matrix[len(v)][len(w)]
    
    print(score)
    print(va)
    print(wa)
      
if __name__ == '__main__':
    main()
    

