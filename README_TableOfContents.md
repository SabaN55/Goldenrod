# Goldenrod

Table of Contents:

N.B. This respository is a small collection of machine learning / bioinformatics python projects I have completed while
in university. It does not represent an exhaustive list of my coding knowledge.



BuildTree.py:

BuildTree is a program that takes sequences and roots them in the least parsimonious tree possible using either
UPGMA or Neighbor-joining.

UltrametricAdditive.py:

This file acts as support to BuildTree.py, it contains the functions is_ultramertric() and is_additive. 
These two receive as input the pairwise distances for a set of sequences and then verify whether the given 
distances satisfy the ultrametricity or additivity criterion, respectively.

generate_HMM_sequence.py:

This function (generate_HMM_sequence) generates a random sequence of characters (in this case, nucleotides), 
while accounting for the random sequence of some hidden states  (which evolve with Markovian transition probabilities).
In this project, our goal was to generate a novel and random genome that would model the C-G / A-T rich segments of
DNA in specific species using a set of given state-dependent emission probabilities.

analyze_sequence.py:

Analyze sequence takes a HMM sequence, such as the one created in generate_HMM_sequence.py and a state-dependent emission
probability matrix to search for a specific sequence query and calculate the probability of different hidden state sequences
for any and all queries found.

simulate.py:

Simulate.py is a computational simulation of a whole-genome shotgun sequencer. Using length of a genome, number of reads, and
the length of the reads, the program is able to statistically model the coverage, number of contigs, length of contigs, etc. of
the shotgun sequence.

assembly_tester.py:

Assembly_tester.py includes functions that work along the same line as simulate.py where .fasta files for contigs and shotgun reads
are passed in and analyzed. The function is able to place each read among the contigs and check to make sure that each paired read
fits into the contigs. If any of the paired reads fail, they are caught by test cases.

viterbi.py:

Viterbi.py is a simulation of the Viterbi algorithm. The dynamic programming function count_segments() determines the number of segments that
exist in each possible hidden state of a Hidden Markov model. In this project, we analyzed the structural patterns of
in the DNA of an ocean-dwelling archaeon in order to predict the locations of structural DNA segments and coding DNA segments.

posterior.py:

These sets of functions serve the same purpose as those in viterbi.py but via a different algorithm. The posterior
dynamic method accounts for all the probabilities that the sequence was generated from each state for every position 
in the genome. 

