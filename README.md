Greedy Motif Search Algorithm

Implementation of the Greedy Motif Search algorithm in Python for DNA motif discovery.

Overview

Greedy Motif Search is a bioinformatics algorithm used to identify conserved motifs across multiple DNA sequences.

Features

- Position Weight Matrix (PWM) construction
- Motif scoring
- Probability calculation
- Profile-most probable k-mer selection
- Greedy motif search

Input

- Multiple DNA sequences
- k-mer length (k)

Output

- Best motif set
- Motif score

Example

Input:

3
AGCTAGC
CAGCTAA
TTGCTAG
3

Output:

Motif: ['GCT', 'GCT', 'GCT']
Score: 0

Concepts Used

- Bioinformatics
- DNA Motif Discovery
- Greedy Algorithms
- Probability Matrices
- Python Programming
