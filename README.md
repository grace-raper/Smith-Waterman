# The Smith-Waterman Algorithm
The Smith-Waterman algorithm is a dynamic programming algorithm to find the optimal alignment of two sequences.

This repository provides a basic implementation of the Smith-Waterman algorithm that has been tailored to the issue of finding the optimal alinment of two DNA protein sequences. Sequence aligment is scored via the BLOSUM62 scoring matrix. The BLOSUM62 scoring matrix is the default scoring matrix for [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and was developed by analyzing the frequencies of amino acid substitutions in clusters of related proteins.

*note: this algorithm was implemented for an assigment in a computation biology class I took in early 2021.*


### A Short Primer on DNA Sequence Alignment
Sequence alignment is the procedure of comparing two sequences by searching for a series of characters that are in the same order in all sequences. Two sequences can be aligned by writing them across a page in two rows. Identical characters are placed in the same column, and non identical ones can either be placed in the same column as a mismatch or against a gap (-) in the other sequence. Sequences that are aligned in this manner are said to be similar. Sequence alignment is useful for discovering functional, structural, and evolutionary information in biological sequences.

Take for example the following 2 sequences. Sequence A is the 10 character string: "GACGGATTAG". Sequence B is the 11 character string: "GATCGGAATAG". An good alignment for them is found in the table below:

| Index      | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| -----------|---|---|---|---|---|---|---|---|---|----|----|
| Sequence A | G | A | - | C | G | G | A | T | T | A  | G  | 
| Sequence B | G | A | T | C | G | G | A | A | T | A  | G  | 

The only differences in these two sequences can be found at index 3 & 8. Observe that a gap (-) is introduced in the first sequence to let equal bases align perfectly. The goal of the Smith-Waterman algorithm is to present an efficient way to take any two sequences and determine the best alignment between them. 

The best aligment is determined by assigning each column a score based on wheither it is a match, mismatch, or gap. To illustrate this concept we will score the above aligment where a match is worth +1, mismatch is worth -1 & gap is worth -2. As we can see in the table below, the above alignment will give a total score: 9 × 1 + 1 × (-1) + 1 × (-2) = 6.

| Index      | 1  | 2  | 3  | 4  | 5  | 6  | 7  |  8 |  9 | 10 | 11 |
| -----------|---|---|---|---|---|---|---|---|---|----|----|
| SCORE      | +1 | +1 | -2 | +1 | +1 | +1 | +1 | -1 | +1 | +1 | + 1|

In this implementation, the match, mismatch and gap penalty are determined based on the BLOSUM62 scoring matrix.

### Why Dynamic Programming?

One approach to compute similarity between two sequences is to generate all possible alignments and pick the best one. However, the number of alignments between two sequences is exponential and this will result in a slow algorithm so, Dynamic Programming is used as a technique to produce faster alignment algorithm. Dynamic Programming tries to solve an instance of the problem by using already computed solutions for smaller instances of the same problem. Giving two sequences Seq1 and Seq2 instead of determining the similarity between sequences as a whole, dynamic programming tries to build up the solution by determining all similarities between arbitrary prefixes of the two sequences. The algorithm starts with shorter prefixes and uses previously computed results to solve the problem for larger prefixes.
