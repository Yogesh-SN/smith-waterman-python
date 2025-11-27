###Smith-Waterman - List of Needed Functions

# Multiple tracebacks?? Not accounted for :/ A large differnece between match,
# mismatch scores and gap penalties allows for mostly (atleast for short
# sequences) one traceback

#2 strings
#x -> x+1 columns
#y -> y+1 rows

import numpy as np
import pandas as pd


#                                                                                                            '''Initialization '''


def alignment_match(x):
  #Alignment array where the aligned sub string (seq2) can be placed against seq1
  alignment_array = np.array([
      list(x), ['' for i in range(len(x))]
      ])
  return alignment_array

def mat_creater(x_len, y_len):
    # Matrix initialization using np.zeros
    return np.zeros((y_len, x_len), dtype=int)

def sub_matrix(m, mm):
    # Creating substitution matrix: 4x4 matrix for match/mis-match values
    # This function is not used anywhere in the code; its only kept for the sake
    # Of completness of the Smith-Waterman Algorithim
    submat = np.full((4, 4), mm, dtype=int)
    np.fill_diagonal(submat, m)
    return submat

def match_mismatch_diction(y,x,m,mm):
  # Maintaing a dictionary with Keys = (j,i) (Coordinates of Alignment between 2
  # base pairs), Values = match or mis-match score)
  m_mm_dict = dict()
  for j in range(len(y)):
    for i in range(len(x)):
      if y[j] == x[i]:
        m_mm_dict[(j+1),(i+1)] = m
      else:
        m_mm_dict[(j+1),(i+1)] = mm

  return m_mm_dict



#                                                                                                          ''' Matrix Filling '''


def SW_max_function(y,x,i,j,mat,g,m,mm): #row,column,scoringmatrix,gap_penalty
  triplemdict = match_mismatch_diction(y,x,m,mm)


  #The piece-wise Function which picks out the max score between the 4 parameters
  val_list = [0,(mat[i-1][j-1] + triplemdict[(i,j)]),
   (mat[i-1][j] + g), (mat[i][j-1] + g)]

  # Max between -> 0, Diagonal Cell + match/mismatch, Up Cell + gap penalty, Left Cell + gap penalty
  max_val = max(val_list)

  return (max_val,val_list.index(max_val)) #maintaing index of which value from the list is picked; would be used in making the 'max_where dict' -> will be used in traceback

def SW_mat_filling(x,y,m,mm,g):
  max_where_dict = dict() #initalize a max_where dictionary: Having keys as coordinated (i,j) and values as index of which value (max) from the val_list is picked up
  x_len, y_len = len(x),len(y)
  x_list, y_list = list(x),list(y)
  mat_main = mat_creater(x_len + 1,y_len + 1) #Initializing a matrix
  for i in range(1,len(mat_main)):
    for j in range(1, len(mat_main[i])):
      mat_main[i][j], max_where_dict[(i,j)] = SW_max_function(y,x,i,j,mat_main,gap_pen,m,mm) #Filling matrix with values, after being passed through the SW_max_function
      #and Dictionary Filling

  mat_main_df = pd.DataFrame(mat_main,
                             index = ['Seq2 ↓'] + [x for x in seq2], #Nucleotide Lables Row-wise for Seq2
                             columns = ['Seq1 →'] + [y for y in seq1] #Nucleotide Lables Column-wise for Seq1
                             ) #Made main_matrix as df - For labelling of rows and columns, and a overall better visual of the scored matrix

  ind = np.unravel_index(np.argmax(mat_main), mat_main.shape)
  ind = tuple(int(x) for x in ind) #Stores max value's position in the entire scoring matrix

  return mat_main, max_where_dict, mat_main_df, ind


#                                                                                                        ''' Traceback & Recursion '''



#Recursive Function which Traces back 0, starting from the max scoring cell in the matrix

def recursion_traceback(
        scored_matrix,
        traceback_dict,
        seq1,
        seq2,
        m_match,
        m_mismatch,
        gap_pen,
        i,
        j,
        alignment,
        y_subseq,
        score=0
    ):

    # Stop if score hits zero (Smith-Waterman condition)
    if scored_matrix[i][j] == 0:
        return alignment, score

    direction = traceback_dict[(i, j)]

    # 1  -  if score is coming from Diagonal
    if direction == 1:
        alignment[1][j - 1] = y_subseq[-1]
        new_i, new_j = i - 1, j - 1 #Change Coordinate
        score += m_match if seq1[j - 1] == seq2[i - 1] else m_mismatch
        return recursion_traceback(
            scored_matrix, traceback_dict,
            seq1, seq2,
            m_match, m_mismatch, gap_pen,
            new_i, new_j,
            alignment,
            y_subseq[:-1],
            score
        )

    # 2 - if score is coming from Up (gap in seq1)
    if direction == 2:
        alignment[1][j - 1] = '-'
        new_i, new_j = i - 1, j #Change Coordinate
        score += gap_pen
        return recursion_traceback(
            scored_matrix, traceback_dict,
            seq1, seq2,
            m_match, m_mismatch, gap_pen,
            new_i, new_j,
            alignment,
            y_subseq,   # seq2 index stays
            score
        )

    # 3 - if the score is coming from Left (gap in seq2)
    if direction == 3:
        alignment[1][j - 1] = '-'
        new_i, new_j = i, j - 1 #Change Coordinate
        score += gap_pen
        return recursion_traceback(
            scored_matrix, traceback_dict,
            seq1, seq2,
            m_match, m_mismatch, gap_pen,
            new_i, new_j,
            alignment,
            y_subseq,   # seq2 index stays
            score
        )

def smith_waterman(a, b, m, mm, gap):

    alignment_array = alignment_match(a)
    filled_matrix, tracebackdictionary, filled_matrix_df, CoordMax = SW_mat_filling(
        a, b, m, mm, gap
    )

    i_max,j_max = CoordMax[0], CoordMax[1]

    alignment, score = recursion_traceback(filled_matrix, tracebackdictionary, a, b, m_match,
        mm_mismatch,
        gap_pen,
        i_max,
        j_max,
        alignment_array,
        b,
        score=0
    )


    alignment_df = pd.DataFrame(alignment, index=["Seq1", "Seq2"])

    return filled_matrix_df, alignment_df, score

'''
Example Usage - Find output in the README.md file

seq1 = "CGTATCTCATT" #cols #CGTATCTCATT
seq2 = "TTC" #rows #TTCC
m_match,mm_mismatch, gap_pen  = 5,-2,-1

scored_matrix, pairwise_localalignment, scoring,  = smith_waterman(seq1,seq2,m_match,mm_mismatch,gap_pen)

print(scored_matrix,'\n')
print(pairwise_localalignment,'\n')
print(f'Alignment Score = {scoring}')

'''
