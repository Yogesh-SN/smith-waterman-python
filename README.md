# smith-waterman-python
Python implementation and recreation of the Smith–Waterman (Smith & Waterman (1981)) local alignment algorithm, complete with scoring matrix generation and alignment reconstruction.

From WikiPedia: "_The Smith–Waterman algorithm performs local sequence alignment; that is, for determining similar regions between two strings of nucleic acid sequences. Instead of looking at the entire sequence, the Smith–Waterman algorithm compares segments of all possible lengths and optimizes the similarity measure._"

"_The Smith-Waterman algorithm is widely used in bioinformatics for identifying similar regions between biological sequences, particularly in scenarios where global alignment is not appropriate.
 Key applications include DNA sequence alignment to compare gene sequences, protein sequence alignment to identify similar functional domains, and database searches to find locally similar sequences within large biological databases_"
 
This project contains a from-scratch implementation of the Smith–Waterman algorithm.

Everything - scoring matrix construction, traceback dictionary, and final alignment - is built manually for full transparency and learning value.
(The only libraries used are NumPy and Pandas)

What the Algorithim is desgined to do:

<img width="300" height="420" alt="image" src="https://github.com/user-attachments/assets/803a64b0-deeb-4737-98f8-ebcad7034914" />                        <img width="301" height="240" alt="image" src="https://github.com/user-attachments/assets/1db76a53-b90e-45f4-a021-432ecd0dd4ef" />




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

2. Match/mismatch/(linear) gap scoring system

3. Customizable scoring values

4. Clean traceback implementation

5. Works with recursion-based traceback logic

6. Returns both the alignment and the full scoring matrix

**Usage Example:**

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

**Output:**

<img width="1201" height="537" alt="image" src="https://github.com/user-attachments/assets/9bb74a2c-cdca-4599-9370-8000c3deba41" />




**Notes / Limitations**:

Local alignment can sometimes produce multiple equally optimal tracebacks.
This implementation focuses on retrieving only **one optimal alignment**.

If you want to discourage branching tracebacks, you can assign:

Match score ≫ Mismatch penalty

Match score ≫ Gap penalty

Using a wide gap between match / mismatch / gap scores usually collapses the ambiguity and yields a single, clean optimal path.
