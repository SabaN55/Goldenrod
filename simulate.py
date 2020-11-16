import sys
import random
from compsci260lib import *

def run_simulate():
    """
    Simulates the sequencing process then empirically compute and report 
    the quantities from parts a-c.
    """

    iterations = 20
    G = 3000000
    R = 40000
    L = 450

    # Call simulate(G, R, L) `iteration` times. Report the empirical values 
    # from parts a-c for each iteration
    #
    cov = 0
    not_covered = 0
    num_contigs = 0
    avg_contigl = 0
    for i in range(iterations):
        trun = simulate(G,R,L)
        print(f"----- Results for run {i+1} -----")
        print(f"The coverage is {trun[0]}.")
        cov += trun[0]
        print(f"There are {trun[1]} nucleotides that were not covered.")
        not_covered += trun[1]
        print(f"There are {trun[2]} contigs.")
        num_contigs += trun[2]
        print(f"The average length of these contigs is {trun[3]}.")
        avg_contigl += trun[3]
        print("")


    print(f"----- AVERAGE RESULTS -----")
    print(f"The coverage is {cov/iterations}.")
    print(f"There are {not_covered/iterations} nucleotides that were not covered.")
    print(f"There are {num_contigs/iterations} contigs.")
    print(f"The average length of these contigs is {avg_contigl/iterations}.")
    #


def simulate(G, R, L):
    """
    Simulates one iteration of the sequencing process and empirically compute 
    the empirical coverage (average number of times a nucleotide in the genome
    was sequenced), the number of nucleotides not covered by any read, the 
    number of contigs assuming you can use the oracular assembly algorithm to
    assemble all the reads, and the average length of these contigs.

    Args:
        G (int) - the length of the genome
        R (int) - the number of reads
        L (int) - the length of each read

    Returns
        a tuple of floats:

            (Empirical coverage, 
             Number of nucleotides not covered by any read, 
             Number of contigs,
             Average length of these contigs)
    """

    #
    genome = [0 for _ in range(G)]
    #print(genome)
    for i in range(R):
        start = random.randint(0,G-L)
        for l in range(L):
            indx = start + l
            genome[indx] = genome[indx] + 1
    #print(genome)
    C = sum(genome)/G

    nulls = 0
    open = False
    num_contigs = 0
    contig_l = 0
    for nuc in genome:
        if nuc == 0:
            nulls += 1
            open = False
        else:
            if open:
                contig_l += 1
            else:
                contig_l += 1
                num_contigs += 1
                open = True





    #
    return (C, nulls, num_contigs, contig_l/num_contigs)


if __name__ == '__main__':
    """Call run_simulate(), do not modify"""
    run_simulate()
