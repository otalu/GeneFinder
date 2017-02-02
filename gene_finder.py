# -*- coding: utf-8 -*-
"""
Genez
@author: Onur, the Incompetent

"""

import random


from amino_acids import aa, codons, aa_table   # you may find these useful


from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    # pass

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    # pass
    reverse = []
    index = len(dna)-1
    while index >= 0:
        nucleotide = get_complement(dna[index])
        reverse.append(nucleotide)
        index = index - 1
    return "".join(reverse)

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    # pass

    for item in range(0, len(dna)-2, 3):
        look = dna[item] + dna[item+1] + dna[item+2]
        if look == "TAG" or look == "TAA" or look == "TGA":
            return dna[:(item)]
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    # pass
    item = 0
    non_nested = []
    while item < len(dna):
        if dna[item:item+3] == 'ATG':
            new = rest_of_ORF(dna[item:])
            non_nested.append(new)
            item = item + len(new)
        else:
            item = item + 3
    return (non_nested)


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    # pass
    part1 = find_all_ORFs_oneframe(dna)
    part2 = find_all_ORFs_oneframe(dna[1:])
    part3 = find_all_ORFs_oneframe(dna[2:])
    return part1+part2+part3


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    # pass

    forwards = find_all_ORFs(dna)
    backwards = find_all_ORFs(get_reverse_complement(dna))
    both = forwards + backwards
    return both


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    # pass

    all_ORF = list(find_all_ORFs_both_strands(dna))
    longest = 0
    item = 0
    while item < (len(all_ORF)-1):
        if all_ORF[item] > all_ORF[item + 1]:
            longest = all_ORF[item]
            item += 1
        else:
            longest = all_ORF[item + 1]
            item += 1
        return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    # pass

    lonc = []  # lonc = longest ORF non coding
    for i in range(num_trials):
        dna_shuffled = shuffle_string(dna)
        new_dna = list(dna_shuffled)
        length = len(longest_ORF(new_dna))
        lonc.append(length)
    lonc.sort()
    return lonc(-1)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    # pass
    AA = ''  # empty string
    # dna = list(dna)
    for item in range(0, len(dna) - len(dna) % 3, 3):  # runs loop for every codon
        amino_acid = aa_table[dna[item:item+3]]  # assigns aminoacid the protein
        AA += amino_acid  # adds AA to the list
    return AA


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    # pass
    threshold = longest_ORF_noncoding(dna, 1500)

if __name__ == "__main__":
    import doctest
    # doctest.testmod()
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
