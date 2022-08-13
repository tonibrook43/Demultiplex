#!/usr/bin/env python
# Author: Toni Brooks tonib@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.2"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from asyncio import base_subprocess
from email.mime import base

DNA_bases = "A, C, T, G, N"
RNA_bases = "A, C, G, U, N"

def rev_comp(DNA): #we are looping over the dna string so we are looking at each position in our string and saving it to base
    cDNA=""
    for base in DNA:
        cDNA += base_dict[base] #need to make a list to store my complemented bases
    return cDNA[::-1]

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    phred_score_int = (ord(letter) - 33)
    return phred_score_int

def qual_score(phred_score: str) -> float:
    """Write your own doc string"""
    sum=0
    for letter in phred_score:
        sum += (convert_phred(letter))
    return sum/len(phred_score)

def validate_DNA_seq(DNA):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    #here is a comment
    DNA = DNA.upper()
    return len(DNA) == DNA.count("A") + DNA.count("T") + DNA.count("G") + DNA.count("C")

def validate_base_sequence(base_sequence, RNAflag):
    """Checks and returns if the string base_sequence contains only upper- DNA (T, C, A, and G) characters"""
    seq = base_sequence.upper()
    return len(seq) == (seq.count("U" if RNAflag else "T") + seq.count('C') + seq.count('A') + seq.count('G') + seq.count('N'))
    pass

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def onelinefasta(fasta_file, output_file):
    """Takes lines from a fasta sequence and concatanates it into one line and outputs it to a new file"""
    with open('fasta_file') as fasta_input, open('output_file', 'w') as fasta_output:
        x = []
    for line in fasta_input:
        if line.startswith('>'): 
            if x:
                fasta_output.write(''.join(x) + '\n') 
                x = []
            fasta_output.write(line)
        else:
            x.append(line.strip())
    if x:
        fasta_output.write(''.join(x) + '\n')
    pass

def (base_seq):
    """"Return the percentage of G and C characters in base_seq"""
    seq = base_seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

if __name__ == "__main__":
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    assert gc_content("gaaa") == 0.25
    print("correctly calculated GC content")
    assert validate_base_seq("AATAGAT", False) == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False
    assert validate_DNA_seq("aaaaa") == True, "DNA string not recognized"
    print("Correctly identified a DNA string")
    assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"
    print("Correctly determined non-DNA")
pass

#Most functions pulled from lecture notes#
