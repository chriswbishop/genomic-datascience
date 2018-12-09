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

def find_start_codon(seq, n, start_codon):
    """
    Recursively find start codons in a sequence

    This can be generalized to search recursively for
    any codon sequence ... should make that happen
    """
    # Find the first instance of a start codon
    _seq = seq.seq[n:].upper()

    if n == 0:
        start_codon = []

    index = _seq.find('ATG')

    if index != -1:

        start_codon.append(n + index)

        return find_start_codon(seq, n + index+1, start_codon)

    elif index == -1:

        return start_codon

