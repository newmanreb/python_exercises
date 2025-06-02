# # 2.1: Refer back to exercise 1.4, where we printed a DNA string in blocks, with a space between each
# # block. Now further develop your code so that it displays a DNA string in the style used in GenBank
# # records. So given the DNA sequence:
#
# sequence = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAG"
#
# def genbank_format(sequence, block_size=10, blocks_per_row=6):
#     seq = sequence.lower()
#     row_size = block_size * blocks_per_row
#
#     for start in range(0, len(seq), row_size):
#         row_seq = seq[start:start + row_size]
#
#         blocks = []
#         for i in range(0, len(row_seq), block_size):
#             block = row_seq[i:i + block_size]
#             blocks.append(block)
#
#         formatted_row = ' '.join(blocks)
#         print(f"{start + 1:<9}{formatted_row}")
#
# genbank_format(sequence)

# 2.2: Write a function to translate a DNA sequence into an amino acid sequence (without importing modules).
# Find the standard genetic code table online, and create a dictionary to hold a translational table.
# Use single letter amino-acid codes, and assume coding starts from the first base only.

# translation_sequence = "aggagtaagcccttgcaactggaaatacacccattg"
#
# def translate_dna(sequence):
#     codons = {
#         "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I",
#         "ATC": "I",
#         "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S",
#         "TCG": "S",
#         "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A",
#         "GCC": "A",
#         "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q",
#         "CAG": "Q",
#         "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
#         "TGC": "C",
#         "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R",
#         "AGG": "R",
#         "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
#     }
#     sequence = sequence.upper()
#     protein = ""
#     for each_codon in range(0, len(sequence)-2, 3):
#         codon = sequence[each_codon:each_codon+3]
#         amino_acid = codons.get(codon, "X")
#         protein += amino_acid
#
#     return protein
#
# print(translate_dna(translation_sequence))
#
# challenge_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAA"
# print(translate_dna(challenge_sequence))

## 2.3: Write a function which generates the reverse complement of a sequence. Bonus points for dealing with gaps or
## incorrect base letters.

def reverse_complement(sequence):
    """
    Takes a sequence and returns a reverse complement, which means complement bases and reversed in order."
    N will return N; gaps as - or as a space will return -. Unknown bases will return X.
    """
    complements = {
        "T": "A",
        "A": "T",
        "G": "C",
        "C": "G",
        "N": "N", # ambiguous base
        "X": "X", # unknown base
        "-": "-", # gap
        " ": "-", # treat spaces as gaps
    }
    sequence = sequence.upper()
    complement = ""

    for each_base in sequence:
        base = complements.get(each_base, "X")
        complement += base

    return complement[::-1]

seq_for_rev_complement = "aggagtaagcccttgcaactggaaatacacccattg"
print(reverse_complement(seq_for_rev_complement))