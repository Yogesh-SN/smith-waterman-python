# smith-waterman-python
Python implementation of the Smith–Waterman local alignment algorithm, complete with scoring matrix generation and alignment reconstruction.

This project contains a from-scratch implementation of the Smith–Waterman algorithm for local sequence alignment.
Everything - scoring matrix construction, traceback dictionary, and final alignment - is built manually for full transparency and learning value.


**The code computes:**

1. The scoring matrix

2. A traceback dictionary showing how each cell was derived

3. The optimal local alignment between two sequences

4. The maximum alignment score

5. A pandas DataFrame view of the scoring matrix

6. A pandas DataFrame of the final alignment

It’s designed to be readable, hackable, and suitable for bioinformatcians studying python


**Features**

1. Full matrix setup (rows = seq2, columns = seq1)

2. Match/mismatch/gap scoring system

3. Customizable scoring values

4. Clean traceback implementation

5. Works with recursion-based traceback logic

6. Returns both the alignment and the full scoring matrix

**Usage:**

alignment_df, score, matrix_df = smith_waterman(
    seq1="CGTATCTCATT",
    seq2="TTC",
    m=5,
    mm=-2,
    gap=-1
)

print(alignment_df)
print(score)
print(matrix_df)

**Notes / Limitations**:

Local alignment can sometimes produce multiple equally optimal tracebacks.
This implementation focuses on retrieving only **one optimal alignment**.

If you want to discourage branching tracebacks, you can assign:

Match score ≫ Mismatch penalty

Match score ≫ Gap penalty

Using a wide gap between match / mismatch / gap scores usually collapses the ambiguity and yields a single, clean optimal path.
