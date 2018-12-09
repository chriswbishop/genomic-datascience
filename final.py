from Bio.Seq import Seq
import Bio.SeqIO
import numpy as np

file_name = 'dna.example.fasta'

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

