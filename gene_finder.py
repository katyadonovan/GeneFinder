# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Katya Donovan
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

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
    if nucleotide == 'A':
	       return 'T'
    elif nucleotide == 'T':
	       return 'A'
    elif nucleotide == 'G':
	       return 'C'
    elif nucleotide == 'C':
	       return 'G'


    # TODO: implement this



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
    lengthofdna = len(dna)
    a= ''
    for i in range(lengthofdna):
        a = a + get_complement(dna[lengthofdna-1-i])

    return a


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
    >>> rest_of_ORF("")
    ''
    """

    for i in range(0, len(dna), 3):
        dna1 = dna[i:i+3]
        if dna1 ==  'TAG':
            return dna[:i]
        elif dna1 == 'TAA':
            return dna[:i]
        elif dna1 == 'TGA':
            return dna[:i]
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
    >>>
    """
    orf_list = []
    for i in range(0, len(dna),3):
        codon = dna[i:i+3]
        if codon == 'ATG':
            a = rest_of_ORF(dna[i:])
            orf_list.append(a)
            dna = dna[i+len(a):]
    return orf_list


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
    a = []
    a= a+ find_all_ORFs_oneframe(dna)
    a= a+ find_all_ORFs_oneframe(dna[1:])
    a= a+ find_all_ORFs_oneframe(dna[2:])
    return a


"""    for i in range(0,len(dna),3):
          dna1 = dna[i:i+3]
          if dna1 == 'ATG':
              a.appen(dna(i:))
    for i in range(0,len(dna),3):
        dna1 = dna[i:i+2]
        if dna1 == 'ATG':
            a.append(dna(i:))
"""
def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    a=[]
    a = find_all_ORFs(dna)
    b=  get_reverse_complement(dna)
    c = find_all_ORFs(b)
    a = a+c
    return a

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    non_nested_orfs = find_all_ORFs_both_strands(dna)
    size = 0
    i=0
    a= 0
    while i < len(non_nested_orfs):
        if len(non_nested_orfs[i]) > size:
            size = len(non_nested_orfs[i])
            a= i
        i+=1
    #print("a: ", a)
    #print("non_nested_orfs: ", non_nested_orfs)
    return non_nested_orfs[a]

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the 9specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i=0
    new_ORF_length=0
    while i < num_trials:
        shuffled_dna = shuffle_string(dna)
        real_longest_ORF = len(longest_ORF(shuffled_dna))
        if real_longest_ORF > new_ORF_length:
            new_ORF_length = real_longest_ORF
        i+=1
    #print(new_ORF_length)
    return new_ORF_length

#longest_ORF_noncoding('ATGCGAATGTAGCATCAAA',120)


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
    amino = []
    remainder = len(dna)%3
    new_length = len(dna)-remainder
    print(dna)
    for i in range(0, new_length,3):
        codon = dna[i:i+3]
        amino_acid = aa_table[codon]
        amino.append(amino_acid)
    amino_str=''.join(amino)
    print(amino_str)
    return amino_str


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = 10#longest_ORF_noncoding(dna,1500)
    sequence = []
    protein = []
    for i in find_all_ORFs_both_strands(dna):
        if len(find_all_ORFs_both_strands(i)) > threshold:
            sequence.append(i)
    print(sequence)
    for j in sequence:
        #print ("j: ", j)
        protein.append(coding_strand_to_AA(j))
    return protein


if __name__ == "__main__":
    import doctest
    dna_seq = load_seq('data/X73525.fa')
    print(gene_finder(dna_seq))
    #doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
