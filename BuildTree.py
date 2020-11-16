from compsci260lib import *
from UltrametricAdditive import is_ultrametric, is_additive
from GlobalAlignerPlus import *


def build_tree():
    """
    Read the aligned student sequences, compute the dictionary of distances,
    construct the appropriate phylogenetic tree.
    """

    # Load the aligned student sequences and create the dictionary of distances
    # in the same format as found in UltrametricAdditive.py. You will need to
    # convert the student names to integers ("Student_1" -> 1)
    #
    fseqs = get_fasta_dict("students.fasta")
    rsequences = get_fasta_dict("students.aligned.fasta")
    sequences = {}
    dist = {}
    numlist = []

    for seqkey in rsequences.keys():
        num = seqkey[-1]
        numlist.append(num)
        sequences[int(num)] = rsequences[seqkey]

    # Now we generate the list of pairs that will be used to calculate distances
    tocomp = []

    for num1 in sequences.keys():
        for num2 in sequences.keys():
            if num1 != num2:
                if sorted([num1, num2]) not in tocomp:
                    tocomp.append([num1, num2])
    # print(tocomp)

    for item in tocomp:
        str = f"{item[0]},{item[1]}"
        d = compute_dist(sequences.get(item[0]), sequences.get(item[1]))
        # print(str,d)
        dist[str] = d

    # Uncomment the following line to print out the distance dictionary
    print_dist(dist)

    # Check if the distances are ultrametric
    threshold = 0.003  # problem-specific distance threshold for this problem
    if is_ultrametric(dist, threshold=threshold):
        print("\nThe distance is ultrametric.")
    else:
        print("\nThe distance is not ultrametric.")

    # Check if the distances are additive
    if is_additive(dist, threshold=threshold):
        print("\nThe distance is additive.\n")
    else:
        print("\nThe distance is not additive.\n")

    # Report the Newick representation of the phylogenetic tree    
    #
    print("The Newick representation of the phylogenetic tree is:")
    print(f"{compute_nj_newick(numlist, dist)}")
    #

    # Call `summarize_alignment` to report the differences between the two
    # most similar and two most different sequence alignments
    #
    # CLosest Covid-19:
    print(summarize_alignment(rsequences["Student_3"], rsequences["Student_8"]))

    # Furthest Covid-19
    print(summarize_alignment(rsequences["Student_3"], rsequences["Student_6"]))

    #Furthest overall
    print(summarize_alignment(rsequences["Student_1"], rsequences["Student_6"]))
    #


def compute_nj_newick(seq_names, dist):
    """
    Performs neighbor joining to construct the phylogenetic tree from a
    distance table.

    Args:
        seq_names (list of ints): representing the sequence names.
                e.g. [1, 2, ..] for ['Student_1', 'Student_2', ...]

        dist (dict of str to float): distance table mapping pairs of students 
            to float distances. Refer to UltrametricAdditive.py for examples of
            distance tables.

    Returns:
        the Newick representation of the phylogenetic tree as a string
    """

    # ------------- Implementation of the Neighbor Joining algorithm ----------
    # Keeping track of variable names:

    # dist              - dictionary containing the computed pair-wise
    #                     distances between nodes
    # node_list         - array containing the list of current nodes (L), which
    #                     gradually decreases in size while iterating to build 
    #                     the tree
    # avg_dist          - dictionary containing the averaged distances 
    #                     (r values) for all of the current nodes (those in L)
    # D                 - dictionary containing the 'adjusted' pairwise
    #                     distances between nodes
    # newick            - dictionary to maintain the Newick representation of
    #                     each leaf node or subtree

    # ------------- Initialization: -------------------------------------------

    node_list = []
    for item in seq_names:# L = the list of current nodes
        node_list.append(int(item))
    newick = {}
    for i in range(1, len(node_list) + 1):
        newick[i] = "" + str(i)

    avg_dist = {}  # averaged distances (r values) for all current nodes

    for i in range(1, len(node_list) + 1):

        avg_dist[i] = 0
        for j in range(1, len(node_list) + 1):

            if i != j:
                avg_dist[i] += dist["%d,%d" % (i, j)] if i < j else \
                    dist["%d,%d" % (j, i)]

        avg_dist[i] = avg_dist[i] / (len(node_list) - 2)

    max_node = len(node_list)  # the maximum key used for a node

    #print(avg_dist)

    # -------------- Iteration to build the tree --------------------

    # As long as the list contains more than 2 nodes, iterate 
    while len(node_list) > 2:

        # ---------- Begin your code -------------

        # Compute the 'adjusted' distances between nodes using the original
        # distances (from dist)
        # and averaged distances (from avg_dist)

        # Let D be the dict to contain 'adjusted' distances
        D = {}

        # Loop through each pair of nodes as entered in the dist dict
        # Use the entries from the avg_dist and calculate entries for the D
        # dict

        for item in dist.keys():
            pair = item.split(",")
            D[item] = dist[item] - avg_dist[int(pair[0])] - avg_dist[int(pair[1])]
            # print(item)
            # print(dist[item])
            # print(D[item])

        # Pick the pair i,j in node_list for which adjusted distance D_ij is
        # minimal.
        # You may find the function two_min_in_dict helpful.
        (i, j) = two_min_in_dict(D)  # Replace with your pair

        # Define a new node k and set dist[m,k] for all nodes m in node_list
        # as (dist[i,m] + dist[j,m] - dist[i,j]) / 2

        # max_node had been earlier set to the largest key used for a node
        k = max_node + 1
        max_node += 1
        m = 0
        for ind in range(len(node_list)):
            m = int(node_list[ind])
            comp = {i, j, k, m}
            #print(comp)
            #print(len(comp))
            if len(comp) == 4:
                k1 = sorted((m,k))
                k2 = sorted((m,i))
                k3 = sorted((m,j))
                k4 = sorted((j,i))
                key1 = f"{k1[0]},{k1[1]}"
                key2 = f"{k2[0]},{k2[1]}"
                key3 = f"{k3[0]},{k3[1]}"
                key4 = f"{k4[0]},{k4[1]}"

                dist[key1] = (dist[key2] + dist[key3] - dist[key4]) / 2
        #print(node_list)

        # ---------- End your code -------------

        # Add the new node k to the Newick format representation of the tree
        # with edges of lengths 
        # dik = (dist[i,j] + avg_dist[i] - avg_dist[j])/2
        # djk = dist[i,j]-d[i,k], 
        # joining k to i and j

        d_ik = (dist["%d,%d" % (i, j)] + avg_dist[i] - avg_dist[j]) / 2
        d_jk = dist["%d,%d" % (i, j)] - d_ik
        newick[k] = "(" + newick[i] + ":" + "%.7f" % (d_ik) + "," \
                    + newick[j] + ":" + "%.7f" % (d_jk) + ")"

        # Remove i and j from node_list and add k
        temp = []
        for idx in range(0, len(node_list)):

            if node_list[idx] != i and node_list[idx] != j:
                temp.append(node_list[idx])

        temp.append(k)

        node_list = list(temp)

        # Update the r terms
        if len(node_list) > 2:
            avg_dist[k] = 0
            #print(avg_dist.keys())

            for ind in range(0, len(node_list) - 1):
                m = node_list[ind]
                avg_dist[m] = avg_dist[m] * (len(node_list) - 1)
                avg_dist[m] -= dist["%d,%d" % (m, i)] if m < i \
                    else dist["%d,%d" % (i, m)]
                avg_dist[m] -= dist["%d,%d" % (m, j)] if m < j \
                    else dist["%d,%d" % (j, m)]
                avg_dist[m] += dist["%d,%d" % (m, k)]

                avg_dist[m] /= (len(node_list) - 2)
                avg_dist[k] += dist["%d,%d" % (m, k)]

            avg_dist[k] = avg_dist[k] / (len(node_list) - 2)

        # Remove any elements from the dict that contain nodes i or j
        delete_from_dict(dist, i, j)
        delete_from_dict(D, i, j)
        delete_from_dict(newick, i, j)
        delete_from_dict(avg_dist, i, j)

    # Return the Newick representation
    return ("(" + newick[node_list[0]] + ":" +
            "%.7f" % (list(dist.values())[0]) +
            "," + newick[node_list[1]] + ":0);\n")


def summarize_alignment(seq1, seq2):
    """
    Summarize the alignment between two sequences by computing the number of
    matches, mismatches, and gaps. This code will contain similar logic to 
    the provided `compute_dist` function.

    Note: that we performed multiple sequence alignment to obtain the 
    aligned sequences in students.aligned.fasta. So, for any pair of sequences, 
    you may find a gap at the same place in the two aligned sequences. Gaps 
    should only be counted if they are matched with a non-gap character 
    (ignore a gap aligned to a gap).

    Args:
        seq1 (str): the first sequence, extracted from a multiple sequence
                    alignment
        seq2 (str): the second sequence, extracted from a multiple sequence
                    alignment

    Returns:
        a tuple of the number of (matches, mismatches, gaps) between the two
        sequences
    """

    #
    matches = 0
    mismatches = 0
    gaps = 0
    for i in range(len(seq1)):
        nuc1 = seq1[i]
        nuc2 = seq2[i]
        if nuc1 != "-" or nuc2 != "-":
            if nuc1 != "-" and nuc2 != "-":
                if nuc1 == nuc2:
                    matches += 1
                else:
                    mismatches += 1
            else:
                gaps += 1

    #

    # return the number of matches, mismatches and gaps as a tuple
    return (matches, mismatches, gaps)


########################################################################
# Provided functions for this problem
########################################################################

def compute_dist(seq1, seq2):
    """Returns the distance between two sequences. The distance is computed
    as the ratio between the number of mismatches and the total number
    of matches or mismatches.

    Args:
        seq1 (string) - first sequence
        seq2 (string) - second sequence

    Returns:
        the ratio of mismatches over the total number of matches or mismatches
        as a float
    """

    mismatch = 0
    match_or_mismatch = 0

    for i in range(0, len(seq1)):
        if seq1[i] == "-":
            continue
        elif seq2[i] == "-":
            continue
        elif seq1[i] == seq2[i]:
            match_or_mismatch += 1
        else:
            mismatch += 1
            match_or_mismatch += 1

    return float(mismatch) / match_or_mismatch


def print_dist(dist):
    """
    Print the distance table
    """

    idx = [int(i.split(',')[0]) for i in list(dist.keys())]
    idx.extend([int(i.split(',')[1]) for i in list(dist.keys())])
    max_idx = max(idx)

    print("\t", end=' ')
    for col in range(2, max_idx + 1):
        print("{:>7}".format(col), end=' ')
    print()

    for row in range(1, max_idx):
        print('%d\t' % row, end=' ')
        for col in range(2, row + 1):
            print('       ', end=' ')
        for col in range(row + 1, max_idx + 1):
            print('{:>7.4f}'.format(dist[str(row) + ',' + str(col)]), end=' ')
        print()

    ########################################################################


# Functions used for Neighbor Joining Algorithm
########################################################################

def delete_from_dict(dictionary, i, j):
    """Deletes the dict entries with keys that contain i or j."""

    for k in list(dictionary.keys()):
        ks = [int(_) for _ in str(k).split(',')]
        if i in ks or j in ks:
            del dictionary[k]


def min_in_dict(wiki):
    """Returns the key associated with the minimum value in the dict."""
    import operator
    return min(iter(wiki.items()), key=operator.itemgetter(1))[0]


def two_min_in_dict(dictionary):
    import operator
    sorted_dict = sorted(iter(dictionary.items()), key=operator.itemgetter(1))
    element = sorted_dict[0][0]  # get the first element of the tuple
    (i, j) = element.split(",")
    return (int(i), int(j))


if __name__ == '__main__':
    '''
    dist_1 = {"1,2": 0.3, "1,3": 0.7, "1,4": 0.9,
              "2,3": 0.6, "2,4": 0.8,
              "3,4": 0.6}
    studs = [1,2,3,4]
    compute_nj_newick(dist_1,studs)
    '''
    build_tree()
