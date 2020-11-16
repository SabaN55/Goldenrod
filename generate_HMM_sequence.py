from compsci260lib import *
import sys, random, os


def generate_HMM_sequence(hmm_file, seq_length):
    """
    Load an HMM specification from file, and generate state and observed
    sequences through sampling of the HMM.

    The HMM will have states labeled as strings and observations as characters
    which will be useful when we generate an HMM with nucleotide sequence
    observations.

    An example return sequence may look like:
        (['W', 'S', 'W'], ['A', 'T', 'G']) or
        (['F', 'L', 'L'], ['2', '6', '6'])

    Arguments:
        hmm_file (str): path to the HMM file
        seq_length (int): the length of the sequence we will generate

    Returns:
        a tuple of
        (the state sequence as strings,
         observed sequence as single character strings)
    """

    if not os.path.exists(hmm_file):
        print("Can't open HMM parameter file: %s" % hmm_file)
        return -1

    f = open(hmm_file, "r")

    # read the state names
    states = f.readline().strip().split()

    # read the initial probabilities
    initial_probs = f.readline().strip().split()
    initial_probs = [float(p) for p in initial_probs]

    # read the transition matrix
    transitions = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        transitions[state] = [float(p) for p in matrix_row]

    # read the input alphabet
    input_alphabet = f.readline().strip().split()

    # read the emission matrix
    emission = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        emission[state] = [float(p) for p in matrix_row]
        # normalize
        sum_emission = sum(emission[state])
        emission[state] = [x / sum_emission for x in emission[state]]
    f.close()

    # Initialize the return lists

    ret_states = []
    ret_syms = []

    # We begin by choosing one of the hidden states, assuming they are all equal
    all_states = []
    for item in transitions.keys():
        all_states.append(item)
    first_choice = random.randint(0, len(transitions.keys()) - 1)

    h_state = all_states[first_choice]

    # We now have our initial state, we can begin the run

    while seq_length > 0:
        #print(h_state)

        # Start by setting the emission probability based off our current hidden state
        grab_from = []
        chances = emission[h_state]
        for i in range(len(chances)):
            freq = int(chances[i] * 100)
            sym = input_alphabet[i]
            # print(chances[i]*100, input_alphabet[i])
            grab_from.extend([sym for q in range(freq)])
        # print(grab_from)

        # We now have an iterable list that contains the frequency of all symbols
        the_decider = random.randint(0, len(grab_from) - 1)
        next_ch = grab_from[the_decider]

        # We now have our next character, so we need to update our state and sym lists
        ret_states.append(h_state)
        ret_syms.append(next_ch)

        # And now we need to see if we change states
        change_chances = transitions[h_state]
        change_grab = []
        for z in range(len(change_chances)):
            freq = int(change_chances[z] * 100)
            tstate = all_states[z]
            change_grab.extend([tstate for q in range(freq)])
        # print(change_grab)

        # Now we choose again
        change_decider = random.randint(0, len(change_grab) - 1)
        h_state = change_grab[change_decider]

        seq_length = seq_length - 1

    # ACGT
    #print(transitions)
    #print(emission)
    #print(input_alphabet)
    #print(len(ret_states))
    #
    return (ret_states, ret_syms)  # a tuple containing the state and observation sequences
