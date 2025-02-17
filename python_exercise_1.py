## 1.1
a = 5
b = 6
c = a**2 + b**2
print(c)
# What if a = 10 and b = 12?
d = 10
e = 12
f = d**2 + e**2
print(f)

## 1.2
string = "TheUniversityOfManchesterFacultyofBiologyMedicineAndHealth"
print(string[3:13], string[15:25])

dna_seq = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTGAT"
print(dna_seq[137:143], dna_seq[274:280])

## 1.3
# Given two integers a and b, return sum of all odd integers from a to be inclusively.
a, b = 50, 100  # example values
total = 0   # initialise sum

for i in range(a,b + 1):    # loop through range of a and b + 1 so it's inclusive
    if i % 2 != 0:  # check if the number is odd
        total += i  # add the odd number to the total

print(total)    # output = 1875

# What would the output be if a = 10 and b = 25?
a, b = 10, 25  # example values
total = 0   # initialise sum

for i in range(a,b + 1):    # loop through range of a and b + 1 so it's inclusive
    if i % 2 != 0:  # check if the number is odd
        total += i  # add the odd number to the total

print(total)    # output = 144

## 1.4
# Given a string representing a DNA sequence, print in blocks, with gaps every so many spaces
# Given the following sequence print with block size 10:
sequence = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTGAT"
block_size = 10 # length of the blocks we want
sequence_length = len(sequence) # get the length of the input sequence
index = 0   # start index at 0

while index < sequence_length:  # loop through the length of the sequence
    block = sequence[index:index + block_size]  # extract a block of characters by block_size
    print(block)    # print that block
    index += block_size # increase the index by the block size to keep going

# Another way
sequence = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTGAT"
block_size = 10
sequence_length = len(sequence)

for i in range(0, len(sequence), block_size):
    print(sequence[i:i + block_size])

## 1.5
# Write code to transcribe a DNA sequence to RNA.
