from Bio.Seq import Seq
import Bio.SeqIO
import numpy as np

# file_name = 'dna.example.fasta'
file_name = 'dna2.fasta'

# Read the file
data = list(Bio.SeqIO.parse(file_name, 'fasta'))

print(f'{len(data)} records')

def sequence_lengths(data):

    return np.array([len(d) for d in data])

def longest_sequence(data):

    return sequence_lengths(data).max()

def mask_sequence_length(data, sequence_length):

    # Get sequence lengths
    _l = sequence_lengths(data)

    return [data[i] for i in range(len(_l)) if _l[i] == sequence_length]

def get_longest_sequence(data):

    # Get the length
    _max = longest_sequence(data)

    # Mask and return value
    return mask_sequence_length(data, _max)

def get_shortest_sequence(data):

    _min = shortest_sequence(data)

    return mask_sequence_length(data, _min)

def shortest_sequence(data):

    return sequence_lengths(data).min()

def n_shortest_sequence(data):

    _l = sequence_lengths(data)

    _min = shortest_sequence(data)

    return (_l == _min).sum()

def n_longest_sequence(data):

    _l = sequence_lengths(data)

    _max = longest_sequence(data)

    return (_l == _max).sum()

def set_reading_frame(seq, offset):

    return seq[offset:]

def is_start_codon(codon):

    if codon.lower() in ['atg']:

        return True

    else:

        return False

def is_stop_codon(codon):

    if codon.lower() in ['taa', 'tag', 'tga']:

        return True

    else:

        return False


def find_codon(seq, codon, n=0, codon_index=[]):
    """
    Recursively searches a sequence and returns the
    start location of the specified codon
    This can be generalized to search recursively for
    any codon sequence ... should make that happen

    Args:
        seq (sequence): nucleotide sequence
        codon (str): triplet codon to search for
        n (int): starting index for search (0)
        codon_index (list): list to store start locations
                            Used as a passthrough for index
                            tracking/appending.

    Returns:
        codon_index (list): list of start locations
    """

    # Find the first instance of a start codon
    _seq = seq[n:].upper()

    if n == 0:

        codon_index = []

        # Case standardization
        codon = codon.upper()

    # Find the next instance of the specified codon
    index = _seq.find(codon)

    if index != -1:

        codon_index.append(n + index)

        return find_codon(seq, codon, n + index+1, codon_index)

    elif index == -1:

        return codon_index

def find_start_codons(seq):
    """
    Find the start location of start codons

    Args:
        seq (sequence): nucleotide sequence

    Returns:
        codon_index (list): list of start codon start locations
    """

    return find_codon(seq, 'ATG')

def find_stop_codons(seq):
    """
    Return start location of stop codons

    Args:
        seq (sequence):

    Returns:
        codon_index (list): list of start codon locations
    """

    # Need to iterate over stop codons
    codon_index = []

    for codon in ['TAA', 'TAG', 'TGA']:

        codon_index += find_codon(seq, codon)

    return codon_index

def is_orf(start_index, stop_index):
    """
    Needed a way to determine if a start and stop codon
    are in frame or not.

    Args:
        start_index (int): start position of start_codon
        stop_index (int): stop position of stop_codon

    Returns:
        boolean: True if in frame, false if not
    """

    # Will be an ORF if the stop codon is in frame
    # and the stop_index is after the start_index
    return (start_index < stop_index) & (((stop_index - start_index) % 3) == 0)

def orf_length(start_index, stop_index):
    """
    Returns length of open reading frame (ORF)

    Args:
        start_index (int): start location
        stop_index (int): stop location

    Returns:
        l (int): length of ORF (-1 if invalid ORF)
    """

    if is_orf(start_index, stop_index):

        return stop_index - start_index

    else:

        return -1

def orf_reading_frame(start_index, stop_index):
    """
    Return reading frame of ORF

    Args:
        start_index (int): location of start_codon
        stop_index (int): location of stop_codon

    Returns:
        reading_frame (int): [1, 2, 3] if valid ORF, -1 otherwise
                             Note: 1, 2, 3 to match indexing in quiz
    """

    if is_orf(start_index, stop_index):

        # +1 to reflect indexing in quiz and lectures
        # Otherwise indexing would be [0, 1, 2] which
        # will lead to errors in responses
        return (start_index % 3) + 1

    else:

        return -1

def get_orf_reading_frame(seq):
    """
    Needed a wrapper to get ORF reading frames
    """

    orfs = get_orfs(seq)

    reading_frame = [orf_reading_frame(start, stop) for start, stop in orfs]

    return reading_frame

def get_orfs(seq):
    """
    Return a list of valid open reading frames (ORFs).

    Args:
        seq (sequence): nucleotide sequence

    Returns
        orf_index (list of lists): each element contains a list of [start, stop]
    """

    start_codons = find_start_codons(seq)

    stop_codons = find_stop_codons(seq)

    # Need to create comprehensive comparison of start/stop codons
    # Added in -1 filter so we only return valid ORFs
    # orf_index = [[start, stop] for start in start_codons for stop in stop_codons if orf_length(start, stop) != -1]

    # After looking at the quiz, I think the instructors want us to assume
    # absolute termination at the nearest stop codon. This is not strictly
    # the case in reality, but ... fine.
    def nearest_stop_codon(start):

        stop = np.array(stop_codons)

        # Only include stop codons that are past
        # the start codon
        stop = stop[stop > start+3]

        # We are only interested in stop locations
        # that are in frame. Apply mask so only
        # that subset is available.
        stop = stop[is_orf(start, stop)]

        if len(stop) > 0:
            # Need to add 3 here so we don't find
            # stop codons that are part of the
            # start codon
            return stop.min()
        else:
            return None

    orf_index = []

    for start in start_codons:

        stop = nearest_stop_codon(start)

        if stop is not None:

            if orf_length(start, stop) != -1:

                orf_index.append([start, stop])

    return orf_index

def get_orf_lengths(seq):
    """
    Needed a way to get all lengths of valid ORFs in a
    specified sequence

    Args:
        seq (sequence): nucleotide sequence

    Returns:
        l (list): List of all ORF lengths
    """

    orf_index = get_orfs(seq)

    # The lengths must include the stop codon too, so add 3
    # Recall that the stop index reflects the beginning of the
    # stop codon.
    l = [x[1] - x[0] + 3 for x in orf_index]

    return l

def get_orf_complete(seq):
    """
    Needed a way to package up ORF information into a dictionary

    Key is (start, stop)
    This includes the following:
        - length (np.array): ORF lengths
        - reading_frame (int): ORF reading frame
    """

    # Get the ORFs
    orfs = get_orfs(seq)

    # Get lengths
    l = get_orf_lengths(seq)

    # Get reading frame
    reading_frame = get_orf_reading_frame(seq)

    # Create dictionary output
    orf_complete = {}

    for i, index in enumerate(iter(orfs)):
        start, stop = index
        orf_complete[(start, stop)] = {}
        orf_complete[(start, stop)]['length'] = l[i]
        orf_complete[(start, stop)]['reading_frame'] = reading_frame[i]

    return orf_complete

def get_fasta_orfs(data):
    """
    Simple wrapper to process multiple sequences read in
    from fasta (or other) file

    Args:
        data (list): list of sequence objects

    Returns:
        orfs (dict): dictionary containing info about all ORFs
                     in fasta or other file
    """

    fasta_orfs = {}

    for d in data:

        _id = d.id

        fasta_orfs[_id] = get_orf_complete(d.seq)

    return fasta_orfs


def fragment_seq(seq, n, offset=0, fragments=[]):
    """
    Breaks a sequence up into fragments of length n
    """

    if offset == 0:

        fragments = []

    # Converting to an immutable object type
    fragment = str(seq[offset:offset+n])

    # print('"' + fragment + '"')
    new_offset = offset + 1
    # print(new_offset)
    # Only include the fragment if it is correctly sized
    if (new_offset <= len(seq)) & (len(fragment) == n):

        fragments.append(fragment)

        return fragment_seq(seq, n, new_offset, fragments)

    else:

        return fragments

def get_fragment_counts(seq, n):
    """
    Needed to count the number of occurances
    of each repeat in a sequence.

    Args:
        seq (sequence): nucleotide sequence
        n (int): repeat length

    Returns:
        fragment_counts (dict): dictionary of {'<fragment>': count}
    """

    fragment_counts = {}

    # Get the repeats
    fragments = fragment_seq(seq, n)

    # Get unique fragments
    #  This will be used to initialize dictionary counts
    #  to 0
    for f in iter(set(fragments)):
        fragment_counts[f] = 0

    # Now, count the fragments
    for f in fragments:
        fragment_counts[f] += 1

    return fragment_counts

def get_fragment_counts_complete(data):
    """
    Wrapper to iterate through all elements in FASTA file

    Args:
        data (list): list of sequences

    Returns:
        fragment_counts (dict):
    """

    pass
