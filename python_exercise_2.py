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
#         "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
#         "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T",
#         "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
#         "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D",
#         "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
#         "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
#         "GGG": "G",
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

# ## 2.3: Write a function which generates the reverse complement of a sequence. Bonus points for dealing with gaps or
# ## incorrect base letters.
#
# def reverse_complement(sequence):
#     """
#     Takes a sequence and returns a reverse complement, which means complement bases and reversed in order."
#     N will return N; gaps as - or as a space will return -. Unknown bases will return X.
#     """
#     complements = {
#         "T": "A",
#         "A": "T",
#         "G": "C",
#         "C": "G",
#         "N": "N", # ambiguous base
#         "X": "X", # unknown base
#         "-": "-", # gap
#         " ": "-", # treat spaces as gaps
#     }
#     sequence = sequence.upper()
#     complement = ""
#
#     for each_base in sequence:
#         base = complements.get(each_base, "X")
#         complement += base
#
#     return complement[::-1]
#
# seq_for_rev_complement = "aggagtaagcccttgcaactggaaatacacccattg"
# print(reverse_complement(seq_for_rev_complement))

## 2.4  Combine translation and reverse complement functions to generate a six-frame translation of
## a DNA sequence. This means you should translate three forward reading frames starting at the
## first, second, and third base of the first codon of the forward sequence, and three reverse
## reading frames starting at the first, second, and third base of the first codon of the reverse
## complement of the sequence.

# def forward_and_reverse(sequence):
#     codons = {
#         "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I",
#         "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
#         "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T",
#         "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
#         "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D",
#         "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
#         "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
#         "GGG": "G",
#     }
#
#     complements = {
#         "T": "A",
#         "A": "T",
#         "G": "C",
#         "C": "G",
#         "N": "N", # ambiguous base
#         "X": "X", # unknown base
#         "-": "-", # gap
#         " ": "-", # treat spaces as gaps
#     }
#
#     sequence = sequence.upper()
#
#     complement = ""
#     for each_base in sequence:
#         base = complements.get(each_base, "X")
#         complement += base
#
#     revcomp = complement[::-1]
#
#     protein1 = ""
#     protein2 = ""
#     protein3 = ""
#     protein4 = ""
#     protein5 = ""
#     protein6 = ""
#
#     for each_base in range(len(sequence) - 2):
#         codon = sequence[each_base:each_base+3]
#         amino_acid = codons.get(codon, "X")
#         frame = each_base % 3
#         if frame == 0:
#             protein1 += amino_acid
#         elif frame == 1:
#             protein2 += amino_acid
#         else:
#             protein3 += amino_acid
#
#     for each_base in range(len(revcomp) - 2):
#         codon = revcomp[each_base:each_base+3]
#         amino_acid = codons.get(codon, "X")
#         frame = each_base % 3
#         if frame == 0:
#             protein4 += amino_acid
#         elif frame == 1:
#             protein5 += amino_acid
#         else:
#             protein6 += amino_acid
#
#     return protein1, protein2, protein3, protein4, protein5, protein6
#
#
# seq = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCC"
#
# p1, p2, p3, p4, p5, p6 = forward_and_reverse(seq)
# print("Forward")
# print("1", p1)
# print("2", p2)
# print("3", p3)
# print("Reverse")
# print("4", p4)
# print("5", p5)
# print("6", p6)

## 2.5. Count single, di-nucleotide and tri-nucleotides in a sequence.

def count_bases(sequence):
    adenosine = 0
    cytosine = 0
    guanine = 0
    thymine = 0
    undefined = 0

    for each_base in sequence:
        if each_base == "a":
            adenosine += 1
        elif each_base == "c":
            cytosine += 1
        elif each_base == "g":
            guanine += 1
        elif each_base == "t":
            thymine += 1
        else:
            undefined += 1
            
    print("A:", adenosine)
    print("C:", cytosine)
    print("G:", guanine)
    print("T:", thymine)
    print("U:", undefined)

count_bases("aggagtaagcccttgcaactggaaatacacccattg")

