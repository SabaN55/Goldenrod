
import sys
from math import log
from compsci260lib import *


def run_viterbi():
    hmm_file = 'HMM.methanococcus.txt'
    input_file = 'bacterial.genome.fasta'
    print("Decoding the Viterbi path for %s using %s" % (input_file, hmm_file))

    vit_ret = viterbi_decoding(input_file, hmm_file)

    # Collect the segment counts for each state
    counts = count_segments(vit_ret)

    # Report the first 10 and last 10 segments of your decoding
    #
    print("The first 10 segments are:")
    topsy_turv = vit_ret[::-1]
    for i in range(10):
        print(topsy_turv[i])
    print("")
    print("The last 10 segments are:")
    for i in range(10):
        print(vit_ret[i])
    print("")

    #

    # Then report the number of segments that exist in each state.
    #
    for item in counts.keys():
        print(f"There are {counts[item]} segments in {item}.")
    #


def viterbi_decoding(input_file, hmm_file):
    """Calculate the viterbi decoding of an input sequence

    Arguments:
        input_file (str): path to input fasta file
        hmm_file (str): path to HMM file

    Returns:
        A list of dictionaries of segments in each state (1-indexed).
        An example output may look like:

        [
            {‘start’: 1, ‘end’: 12, ‘state’: ‘state2’},
            {‘start’: 13, ‘end’: 20, ‘state’: ‘state1’},
            ...
        ]
    """

    # open hmm file
    try:
        f_hmm_file = open(hmm_file, 'r')
    except IOError:
        print("IOError: Unable to open HMM file: %s." % hmm_file)
        print("Exiting.")
        sys.exit(1)

    # read the state names
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]
    
    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
    
    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()
    
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1
    
    print("Read a sequence of length", len(emit_str))
    
    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)]

    # Initialize the first column of the matrix
    for i in range(K):
        in_index = get_emit_index(emit_str[0].upper(), emitted_symbols)
        viterbi[i][0] = log(emit_probs[i][in_index]) + log(initial_probs[i])
    
    # Build the matrix column by column
    for j in range(1, len(emit_str)):
        in_index = get_emit_index(emit_str[j].upper(), emitted_symbols)
        
        for i in range(K):
            # Compute the entries viterbi[i][j] and pointers[i][j]
            # Tip: Use float('-inf') for the value of negative infinity
            #
            max_p = float('-inf')
            pointer = 0
            for q in range(K):
                akl = log(transitions[q][i])
                el = log(emit_probs[i][in_index])
                vk = viterbi[q][j-1]
                if akl+el+vk > max_p:
                    max_p = akl+el+vk
                    pointer = q
            viterbi[i][j] = max_p
            pointers[i][j] = pointer
            #

    # Traceback, stored as a list segment lengths in each state as dictionaries
    #
    f_val = float('-inf')
    start = 0
    start_p = 0
    for i in range(K):
        if viterbi[i][len(emit_str)-1] > f_val:
            f_val = viterbi[i][len(emit_str)-1]
            start = i
            start_p = pointers[i][len(emit_str)-1]
    path = []

    c_p = start_p

    #print(pointers)
    sub_from = len(emit_str)-1

    # These are just debugging variables
    c_0 = 0
    c_1 = 0

    retlist = []


    # Here the path is created as we iterate through the input string backwards
    cur_state = states[c_p]
    cur_begin = sub_from +1
    while sub_from >= 0:
        if c_p == 0:
            c_0 += 1
        else:
            c_1 += 1

        if states[c_p] != cur_state:
            end = sub_from + 2
            retlist.append({
                "start": end,
                "end": cur_begin,
                "state": cur_state
            })
            cur_state = states[c_p]
            cur_begin = sub_from +1


        path.append(states[c_p])
        sub_from -= 1
        c_p = pointers[c_p][sub_from]
    retlist.append({
        "start": 1,
        "end": cur_begin,
        "state": cur_state
    })
    return retlist


def count_segments(vit_ret):
    """Calculate the number of segments appearing in each state of
    the viterbi path

    Arguments:
        vit_ret (list of dicts): dictionary of lengths in each state.
            see: return value of viterbi_decoding

    Returns:
        a dictionary of states to number of occurrences in the state. 
        e.g. {'state1': 10, 'state2': 9}
    """

    #
    ret_dict = {}
    for item in vit_ret:
        state = item["state"]
        if state not in ret_dict.keys():
            ret_dict[state] = 1
        else:
            ret_dict[state] += 1
    #
    #print(ret_dict)
    return ret_dict


def get_emit_index(input_val, alphabet):
    """Get the index of the emission value in the alphabet

    Note: This will be a useful function for indexing the emission
    probabilities"""
    return alphabet.index(input_val)


if __name__ == '__main__':
    """Call run_viterbi(), do not modify"""
    run_viterbi()
