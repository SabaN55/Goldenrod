from compsci260lib import *


def run_assembly_tester():
    """
    Read in the contig files, verify the reads in each contig, and estimate
    the gap length between the contigs.
    """
    reads_file = 'paired.reads.fasta'
    contig0_file = 'contig0.fasta'
    contig1_file = 'contig1.fasta'

    # load the fasta files
    reads = get_fasta_dict(reads_file)
    contig0 = get_fasta_dict(contig0_file)
    contig1 = get_fasta_dict(contig1_file)

    # Determine how the reads map to the contigs
    contig_reads = find_contig_reads(reads, contig0, contig1)

    print(contig_reads["seq37"])

    # Report the reads whose ends map to both contig0 and contig1
    # and their estimated gap lengths
    #
    # Case1: sequence doesnt match with a contig
    case1 = True
    for i in contig_reads:
        if contig_reads[i]['contig_a'] == None and contig_reads[i]['contig_b'] == None:
            case1 = False
            print(f'{i} FAILS C1: does not match up with any contigs.')
    if case1:
        print("PASSED C1: Every read appears in the correct orientation somewhere in the set of contigs.")
    print("")
    #



    case2 = True
    for i in contig_reads:
        if contig_reads[i]['contig_a'] == contig_reads[i]['contig_b'] and contig_reads[i]['contig_a'] != None:
            if contig_reads[i]['start_a'] > contig_reads[i]['start_b']:
                if contig_reads[i]['end_a'] - contig_reads[i]['start_b'] > 2010 or contig_reads[i]['end_a'] - \
                        contig_reads[i]['start_b'] < 1990:
                    case2 = False
                    print(f'{i} FAILS C2: Paired reads are not 2000 +- 10 long.')
            if contig_reads[i]['start_b'] > contig_reads[i]['start_a']:
                if contig_reads[i]['end_b'] - contig_reads[i]['start_a'] > 2010 or contig_reads[i]['end_b'] - \
                        contig_reads[i]['start_a'] < 1990:
                    case2 = False
                    print(f'{i} FAILS C2: Paired reads are not 2000 +- 10 long.')
    if case2:
        print(
            "PASSED C2: Whenever a mated pair of reads appears in the same contig, the distance from the beginning of the first read to the end of the last read is 2000 ± 10.")
    print("")
    #
    # Num_multi_contig is the count of sequences where each read appears in a seperate contig.
    # This test is also important because we will be determining gap length for all cases where the reads are in two different contigs

    multi_contigs = []
    case3 = True
    for i in contig_reads:
        if contig_reads[i]['contig_a'] != contig_reads[i]['contig_b'] and contig_reads[i]['contig_a'] != None and \
                contig_reads[i]['contig_b'] != None:
            if contig_reads[i]['contig_a'] != 'contig1':
                case3 = False
                print(f'{i} FAILS C3: First read is not in the first contig.')
            else:
                # Here is where we calculate the gap lengths

                readl = len(reads[i+'a']) + len(reads[i+'b'])
                up = (len(contig1['contig1']) - contig_reads[i]['end_a'])
                down = contig_reads[i]['start_b']

                multi_contigs.append((i, 2000-readl-(up+down)))

    if case3:
        print(
            'PASSED C3: Whenever a mated pair of reads appears in different contigs (but in a single supercontig) that the first read in the pair appears before the second read in the pair.')
    print("")

    print(f'There were {len(multi_contigs)} mated pairs of reads where one read maps to contig0 and the other maps to contig1:')
    for i in multi_contigs:
        print(f'There is a gap between {i[0]}a and {i[0]}b: it is {i[1]} +/- 10 nucleotides long.')


def find_contig_reads(reads, contig0, contig1):
    """
    Determine whether the sequencing reads map to contig0/contig1/both/neither
    and where in the contig it matches. `reads` will contain both ends 'a' and
    'b', but you will return a dictionary using the read name as a whole
    (without 'a' and 'b').

    It should return a dictionary mapping the name of the mated pair of reads
    to:
        - 'contig_a' (str): the contig in which read 'a' was found, as 
          `contig0', `contig1', or None.
        - start_a (int): the start position (1-indexed) read end 'a' mapped to
          within its respective contig (None if not found in any contig)
        - end_a (int): the end position (1-indexed) read end 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'contig_b' (str): the contig in which read 'b' was found, as
          `contig0', `contig1', or None.
        - start_b (int): the start position (1-indexed) read end 'b' mapped to
          within its respective contig (None if not found in any contig)
        - end_b (int): the end position (1-indexed) read end 'b' mapped to
          within its respective contig (None if not found in any contig)

    The returned dictionary should look something like:
    {
        'seq1': {
            'contig_a': 'contig0',
            'start_a': 301,
            'end_a': 800,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        'seq2': {
            'contig_a': 'contig1',
            'start_a': 1101,
            'end_a': 1600,
            'contig_b': 'contig0'
            'start_b': 201,
            'end_b': 700,
        },
        'seq3' : {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        ...
    }

    Arguments:
        reads (dict str to str): dictionary of sequencing reads
        contig0 (dict str to str): dictionary of reads in contig0
        contig1 (dict str to str): dictionary of reads in contig1

        see: get_fasta_dict

    Returns:
        Dictionary mapping reads to information about their locations in contigs.
    """

    #
    retdict = {}

    sequences = set()
    for i in reads:
        sequences.add(i[0:len(i) - 1])

    for read in sequences:
        father_seq = read
        read_seq_a = reads[father_seq + 'a']
        read_seq_b = reads[father_seq + 'b']
        # if len(read_seq_a) != len(read_seq_b): print("oh n0")
        # frameshift = len(read_seq_a)
        retdict[read] = {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None,
            'start_b': None,
            'end_b': None
        }

        # START WITH THE A
        frameshift = len(read_seq_a)
        # FIRST CONTIG
        first_contig = contig0['contig0']

        for i in range(len(first_contig) - frameshift):
            segment = first_contig[i:i + frameshift]
            if segment == read_seq_a:
                retdict[read]['contig_a'] = 'contig0'
                retdict[read]['start_a'] = i + 1
                retdict[read]['end_a'] = i + frameshift + 1

        sec_contig = contig1['contig1']

        for i in range(len(sec_contig) - frameshift):
            segment = sec_contig[i:i + frameshift]
            if segment == read_seq_a:
                retdict[read]['contig_a'] = 'contig1'
                retdict[read]['start_a'] = i + 1
                retdict[read]['end_a'] = i + frameshift + 1

        # NOW B
        frameshift = len(read_seq_b)
        # FIRST CONTIG
        first_contig = contig0['contig0']

        for i in range(len(first_contig) - frameshift):
            segment = first_contig[i:i + frameshift]
            if segment == read_seq_b:
                retdict[read]['contig_b'] = 'contig0'
                retdict[read]['start_b'] = i + 1
                retdict[read]['end_b'] = i + frameshift + 1

        sec_contig = contig1['contig1']

        for i in range(len(sec_contig) - frameshift):
            segment = sec_contig[i:i + frameshift]
            if segment == read_seq_b:
                retdict[read]['contig_b'] = 'contig1'
                retdict[read]['start_b'] = i + 1
                retdict[read]['end_b'] = i + frameshift + 1
    #

    return retdict


if __name__ == '__main__':
    """Call run_assembly_tester(). Do not modify this block."""
    run_assembly_tester()
    
