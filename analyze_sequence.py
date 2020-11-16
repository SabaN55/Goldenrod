from compsci260lib import *
from generate_HMM_sequence import generate_HMM_sequence


def run_analyze_sequence():
    """Load the nucleotide sequence HMM file and generate a 100000 length sequence,
    then do some analysis of the results, including the average length in each state"""
    seq_length = 100000
    (state_sequence, 
    observed_sequence) = generate_HMM_sequence("HMM.parameters.txt",
                                                     seq_length)

    # Write your generated sequence to file. Report the states that prevailed whenever
    # TCGA was observed, and their frequencies when observing TCGA. Finally, report 
    # the average length spent in each state across the entire sequence.
    #
    write_text = True

    #Write the text file:
    if write_text == True:
        filename = "nucleotide_sequence.txt"
        f = open(filename, "w")
        for line in state_sequence:
            f.write(line)
        for line in observed_sequence:
            f.write(line)
        f.close()
        print("Wrote %s" % filename)

    # Report the frequency of each occurance of "TCGA"
    freq_dic = compute_state_frequencies(state_sequence,observed_sequence,"TCGA")

    total_occur = 0
    for item in freq_dic.values():
        total_occur += item

    print(f"The TCGA sequence was found {total_occur} times with these states at these frequencies:")
    for item in sorted(freq_dic):
        print(f"{item}: {freq_dic[item]/total_occur}")

    # Report the frequency of each of the states:
    state_dict = compute_average_length(state_sequence)

    for item in state_dict:
        print(f"The system remains in the {item} state for an average of {state_dict[item]} nucleotides.")
    #


def compute_average_length(state_sequence):
    """
    Given a state sequence, return the average length in each state

    Arguments:
        state_sequence (list of one char strings): generated sequence of states

    Returns:
        a dictionary mapping the state name to the average length in the state.
        Example return:
        {
            'W': 5.5,
            'S': 4.5
        }

    """

    #
    total_run_len = {}
    i = 0
    while i <= len(state_sequence) -1:

        curr_run = 1
        curr_state = state_sequence[i]
        if curr_state not in total_run_len:
            total_run_len[curr_state] = []
        n = 1
        while i+n < len(state_sequence):
            #print(i+n)
            if state_sequence[i+n] == curr_state:
                curr_run += 1
                n += 1
            else:
                break
        total_run_len[curr_state].append(curr_run)
        i += n

    # Now go through this dictionary and generate the avg runs:
    ret_dict = {}
    for item in total_run_len.keys():
        #print(item, len(total_run_len[item]))
        total = 0
        for runl in total_run_len[item]:
            total += runl
        ret_dict[item] = total/len(total_run_len[item])

    return ret_dict
    #


def compute_state_frequencies(state_sequence, observed_sequence, subsequence):
    """Given the state and observed sequences, return the state sequences that
    emitted the query subsequence and frequency in which those state sequences
    sequences emitted the subsequence.

    Arguments:
        state_sequence (list of one char strings): generated sequence of states
        observed_sequencs (list of one char strings): generated observed 
            sequence
        subsequence (list of one char strings): the observed subsequence to
            count

    Returns:
        a dictionary mapping the state name to the frequency of observing the
        provided sequence. Example return:
        {
            'WWWW': 2,
            'WWWS': 1,
            ...
        }
    """

    freq_dic = {}
    skipper = len(subsequence)
    #print(skipper)
    for i in range(len(observed_sequence) - 1):
        seq = "".join(observed_sequence[i:i + skipper])
        if seq == subsequence:
            tstate = "".join(state_sequence[i:i + skipper])
            if tstate not in freq_dic.keys():
                freq_dic[tstate] = 1
            else:
                freq_dic[tstate] += 1


    return freq_dic




if __name__ == '__main__':
    """Main method call, do not modify"""
    run_analyze_sequence()
