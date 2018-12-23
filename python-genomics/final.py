from Bio.Seq import Seq
import Bio.SeqIO
import numpy as np

# file_name = 'dna.example.fasta'
file_name = 'dna2.fasta'

# Read the file
data = list(Bio.SeqIO.parse(file_name, 'fasta'))


def sequence_lengths(data):
    """
    Compute sequence lengths

    Required for final.

    Args:
        data (list): list of sequence objects

    Returns:
        l (np.array): array of lengths
    """

    return np.array([len(d) for d in data])


def longest_sequence(data):
    """
    Compute max length of sequences.

    Required for final.

    Args:
        data (list): list of sequence objects

    Returns:
        max_length (int): maximum observed sequence length
    """

    return sequence_lengths(data).max()


def mask_sequence_length(data, sequence_length):
    """
    Return sequences of specified sequence length.

    Turned out to be a useful function for the final.

    Args:
        data (list): list of sequence objects
        sequence_length (int): target sequence length

    Returns:
        data (list): subsetted list of sequences of
                     desired length
    """

    # Get sequence lengths
    _l = sequence_lengths(data)

    return [data[i] for i in range(len(_l)) if _l[i] == sequence_length]


def get_longest_sequence(data):
    """
    Get ... longest ... sequence

    Args:
        data (list): list of sequence objects

    Returns:
        data (list): subsetted data with longest length
    """

    # Get the length
    _max = longest_sequence(data)

    # Mask and return value
    return mask_sequence_length(data, _max)


def get_shortest_sequence(data):
    """
    Get ... shortest ... sequence

    Args:
        data (list): list of sequence objects

    Returns:
        data (list): subsetted data with shortest length
    """

    _min = shortest_sequence(data)

    return mask_sequence_length(data, _min)


def shortest_sequence(data):
    """
    Compute min length of sequences.

    Required for final.

    Args:
        data (list): list of sequence objects

    Returns:
        min_length (int): minimum observed sequence length
    """

    return sequence_lengths(data).min()


def find_codon(seq, codon, n=0, codon_index=[]):
    """
    Recursively searches a sequence and returns the
    start location of the specified codon.

    Note: code supports arbitrary sequences and lengths.
          So, "codon" here is perhaps too application
          specific.

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
    are in the same frame or not. Stop codon must follow
    start codon.

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
        # Special value indicating an invald ORF
        # Used elsewhere for filtering, etc.
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
    Needed a wrapper to get all ORFs in a sequence

    Args:
        seq (sequence): nucleotide sequence

    Returns:
        reading_frame (list): list of reading frames for
                              all ORFs in sequence
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

    # After looking at the quiz, I think the instructors want us to assume
    # absolute termination at the nearest stop codon. This is not strictly
    # the case in reality, but ... fine.
    def nearest_stop_codon(start):

        # Convert to numpy array to facilitate filtering
        stop = np.array(stop_codons)

        # Only include stop codons that are past
        # the start codon. +3 here is a magic number
        # to advance past the triplet codon.
        stop = stop[stop > start+3]

        # We are only interested in stop locations
        # that are in frame. Apply mask so only
        # that subset is available.
        stop = stop[is_orf(start, stop)]

        # Possible that we will not find an in-frame
        # stop code (i.e., no ORF)
        # Conditional statement here protects against
        # this error and facilitates filtering return
        if len(stop) > 0:
            return stop.min()
        else:
            return None

    # Will be a list of lists of values
    orf_index = []

    # Iterate through all start codons and find the
    # nearest, in-frame stop codon
    for start in start_codons:

        stop = nearest_stop_codon(start)

        if stop is not None:

            if orf_length(start, stop) != -1:

                orf_index.append([start, stop])

    return orf_index


def get_orf_lengths(seq):
    """
    Needed a way to get all lengths of valid ORFs in a
    specified sequence.

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
    Needed a way to package up ORF information into a dictionary.
    This high-level data structurew as easier to manipulate for
    the final.

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
    Fragments a sequence into fragments of length n.
    At its core, this Recursively traverses the
    DNA sequence.

    Note: ran into max recursion depth errors,
    had to set sys.setrecursionlimit(5000)

    Args:
        seq (sequence): nucleotide sequence
        n (int): fragment size
        offset (int): starting point offset
        fragments (list): list of DNA fragments.
                          This is primarily used as
                          a passthrough during recursion
    """

    if offset == 0:

        fragments = []

    # Converting to an immutable object string
    #   (Hashing sequence in dictionaries later was
    #    problematic)
    fragment = str(seq[offset:offset+n])

    # Advance to next location in sequence
    new_offset = offset + 1

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


def get_fragment_counts_complete(data, n):
    """
    Wrapper to iterate through all elements in FASTA file

    Args:
        data (list): list of sequences

    Returns:
        fragment_counts (dict):
    """

    fragment_counts = {}

    for i, d in enumerate(iter(data)):
        print(i)
        fragment_counts[d.id] = get_fragment_counts(d.seq, n)

    return fragment_counts


def aggregate_fragment_count(data, n):
    """
    Final questions often required me to aggregate information
    and counts over all sequences. This wrapper function
    helped me do that quickly and reproducibly.

    Args:
        data (list): list of sequences
        n (int): fragment size

    Returns:
        _agg (dict): dictionary of 'seq': count
    """

    fragment_counts = get_fragment_counts_complete(data, n)

    _agg = {}

    # Aggregate counts
    for seq_id, fc in fragment_counts.items():
        for fc_id, count in fc.items():

            if fc_id not in _agg:
                _agg[fc_id] = count
            else:
                _agg[fc_id] += count

    return _agg
