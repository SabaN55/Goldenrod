import sys
from math import log, exp
from compsci260lib import get_fasta_dict


def run_posterior():
    input_file = 'bacterial.genome.fasta'
    hmm_file = "HMM.methanococcus.txt"
    posterior = posterior_decoding(input_file, hmm_file)

    # Report the first and last ten segments in your decoding
    #
    print("The first 10 segments are:")
    for i in range(10):
        print(posterior[i])
    print("")
    print("The last 10 segments are:")
    for i in range(len(posterior) - 10, len(posterior)):
        print(posterior[i])
    print("")

    print(count_segments(posterior))
    #


def posterior_decoding(input_file, hmm_file):
    """
    Calculate the posterior decoding and return the decoded segments.

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

    # Read in the input files
    f_hmm_file = open(hmm_file)
    if f_hmm_file is None: sys.exit("Can't open file: " + hmm_file)

    # read the state names and number
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

    print("Done reading sequence of length ", len(emit_str))

    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emitted_symbols,
                          emit_probs, emit_str)

    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions,
                            emitted_symbols, emit_probs, emit_str)

    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]

    # Create the list of decoded segments to return
    #
    # This section creates the state path
    state_path = []
    twos = 0
    for item in posterior:
        tmax = max(item)
        state = states[item.index(tmax)]
        if state == "state2":
            twos += 1
        state_path.append(state)
    #print(len(state_path))
    #print(twos)

    # Using that state path, we can create the return list
    curr_state = state_path[0]
    retlist = []
    start = 1
    i = 0
    while i <= len(state_path)-1:
        if state_path[i] != curr_state:
            retlist.append({
                "start": start,
                "end": i,
                "state": curr_state,
            })
            start = i + 1
            curr_state = state_path[i]
        i += 1
    retlist.append({
        "start": start,
        "end": i,
        "state": curr_state,
    })
    #
    #print(retlist)
    return retlist


def run_forward(states, initial_probs, transitions, emitted_symbols,
                emit_probs, emit_str):
    """Calculates the forward probability matrix.

    Arguments:
        states (list of str): list of states as strings
        initial_probs (list of float): list of initial probabilities for each
            state
        transitions (list of list of float): matrix of transition probabilities
        emission_symbols (list of str): list of emission symbols
        emit_probs (list of list of float): matrix of emission probabilities
            for each state and emission symbol
        emit_str (str):

    Returns:
        (list of list of floats): matrix of forward probabilities
    """

    K = len(states)
    N = len(emit_str)

    forward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    emit_index = get_emit_index(emit_str[0], emitted_symbols)
    for k in range(K):
        forward[0][k] = log(initial_probs[k]) + log(emit_probs[k][emit_index])

    # Iterate
    for i in range(1, N):
        emit_index = get_emit_index(emit_str[i].upper(), emitted_symbols)

        # Compute the forward probabilities for the states
        #
        for l in range(K):
            tosum = []
            for k in range(K):
                el = log(emit_probs[l][emit_index])
                toe = forward[i - 1][k]
                akl = log(transitions[k][l])
                thisk = toe + akl + el
                tosum.append(thisk)

            forward[i][l] = summer(tosum)

        #
    #print(forward)
    return forward


def run_backward(states, initial_probs, transitions, emitted_symbols,
                 emit_probs, emit_str):
    """Calculates the backward probability matrix.

        Arguments:
            states (list of str): list of states as strings
            initial_probs (list of float): list of initial probabilities for 
                each state
            transitions (list of list of float): matrix of transition
                probabilities
            emission_symbols (list of str): list of emission symbols
            emit_probs (list of list of float): matrix of emission
                probabilities for each state and emission symbol
            emit_str (str):

        Returns:
            (list of list of floats): matrix of backwards probabilities
    """

    K = len(states)
    N = len(emit_str)

    backward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    for k in range(K):
        backward[N - 1][k] = log(1)  # which is zero, but just to be explicit...

    # Iterate
    for i in range(N - 2, -1, -1):
        emit_index = get_emit_index(emit_str[i + 1].upper(), emitted_symbols)

        # Compute the backward probabilities for the states
        #
        for k in range(K):
            tosum = []
            for l in range(K):
                el = log(emit_probs[l][emit_index])
                toe = backward[i + 1][l]
                akl = log(transitions[k][l])
                thisl = toe + akl + el
                tosum.append(thisl)

            backward[i][k] = summer(tosum)
        #
    #print(backward)
    return backward


def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


def summer(probs_list):
    '''
    A helper funciton where i is the current pos,
    l is the state that is being transitioned to
    '''
    p = max(probs_list)
    sum = 0
    for q in probs_list:
        sum += exp(q-p)
    return p + log(sum)

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



if __name__ == '__main__':
    """Call run_posterior(), do not modify"""
    run_posterior()
