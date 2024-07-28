# -*- coding: utf-8 -*-
"""alntomatrix.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1vnXMc1GsEAyaUkmiZwkp7vWSz8khPJoT

# **Now for aligment format without matrices ", "**
"""

# Upload a CLUSTAL alignment file
with open('CLUSTALcons.aln', 'r') as file:
    # Move the file pointer to the 5th character (index 4, since it's 0-based)
    file.seek(8)
    # Read the rest of the file from the current position
    content = file.read()
alignment_data = content[:-1]

print(alignment_data)



print(alignment_data)
# OR define the alignment data
#alignment_data = """
#ENSGGOT00000004685.3/1-4167 ------------------------------------------------------------
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 --AAAAGGAACTTTTTCGTAACAGTTGTACAACCAGGCTAAAGAGCACACCCTCGCCTTC
#ENSPTRT00000078486.1/1-4419 GAAAAAGGAACTTTCTCGTAACAGTTGTACAACCAGGCTAAAGAGCACACCCTCGTCTTC
#ENSPPAT00000042227.1/1-4166 ------------------------------------------------------------
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#
#ENSGGOT00000004685.3/1-4167 ------------------------------------------------------------
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 CCTGAAGGGGAGGAAACCTGCAATAATGTGAGGTAATTAAATTTAGCCTTAGAAAAGTCT
#ENSPTRT00000078486.1/1-4419 CCTGAAGGGGAGGAAACATGCAATAATGTGAGGTAATTAAATTTAGCCTTAGAAAAGTCT
#ENSPPAT00000042227.1/1-4166 ------------------------------------------------------------
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#
#ENSGGOT00000004685.3/1-4167 -------------------------------------------GCTTTATTTTTCAATTC
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 TTGAAAAGTTCTAATCAAGACTTTCAGCTTTCCGATTGTGTGCGCTTTATTTTTCAATTC
#ENSPTRT00000078486.1/1-4419 TTGAAACGTTCTAATCAAGACTTTCAGCTTTCCGATTGTGTGCGCTTTATTTTTCAATTC
#ENSPPAT00000042227.1/1-4166 -------------------------------------------GCTTTATTTTTCAATTC
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#"""



#if seqnumber was wrong, write it manually
#seqnumber = ""
#seqint = int(seqnumber)
#line_count = alignment_data.count('\n')
#x= linecount/(seqint + 1)

# Define the number of sequences
space_count = alignment_data.count('ENS')
line_count = alignment_data.count('\n')
last_newline_index = alignment_data.rfind('\n')
x= line_count-space_count
seqnumber = space_count/x
a = int(seqnumber)
print("seqnumber, as number of sequences is ",seqnumber, "... NOT: If seqnumber was wrong, write it manually and control the output")
print()





# Splitting the alignment data by newline character
lines = alignment_data.split('\n')

# Extracting data after spaces in each line
extracted_data = [line.split(' ')[-1] for line in lines]

# Joining the extracted data back to a string
result = '\n'.join(extracted_data)

# Split the new alignment data by lines

alignment_lines = result.strip().split('\n')



#Define the indices

indices = [(i*a+i)+c for c in range(0, a)for i in range(0, x) ]
print("Indices:", indices)

# Initialize an empty list to store the combined strings
combined_strings = []

# Iterate over the indices
for index in indices:
    # Extract the string at the specified index and append it to the list
    combined_strings.append(alignment_lines[index])

# Combine the strings into one
result = ''.join(combined_strings)

# Print the result
print(result)

character_count = len(result)
print("Character count:", character_count)
int_count = int(character_count)
seq_lenght = int(int_count/a)
print("Sequence lenght:", seq_lenght)


def add_x_every_60(text):
    resultmod2 = ''
    for i in range(0, len(text), seq_lenght):
        resultmod2 += text[i:i+seq_lenght] + '\n'
    return resultmod2



modified_textwo2 = add_x_every_60(result)
print(modified_textwo2)

from os import read
f = open("alignment.txt", "w")
f.write(modified_textwo2)
f.close()

read_file = open("alignment.txt", "r")
print(read_file.read())

"""# **Alignment with Matrices ", "**"""

# Upload a CLUSTAL alignment file
with open('CLUSTALcons.aln', 'r') as file:
    # Move the file pointer to the 5th character (index 4, since it's 0-based)
    file.seek(8)
    # Read the rest of the file from the current position
    content = file.read()
alignment_data = content[:-1]
# OR define the alignment data
#alignment_data = """
#ENSGGOT00000004685.3/1-4167 ------------------------------------------------------------
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 --AAAAGGAACTTTTTCGTAACAGTTGTACAACCAGGCTAAAGAGCACACCCTCGCCTTC
#ENSPTRT00000078486.1/1-4419 GAAAAAGGAACTTTCTCGTAACAGTTGTACAACCAGGCTAAAGAGCACACCCTCGTCTTC
#ENSPPAT00000042227.1/1-4166 ------------------------------------------------------------
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#
#ENSGGOT00000004685.3/1-4167 ------------------------------------------------------------
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 CCTGAAGGGGAGGAAACCTGCAATAATGTGAGGTAATTAAATTTAGCCTTAGAAAAGTCT
#ENSPTRT00000078486.1/1-4419 CCTGAAGGGGAGGAAACATGCAATAATGTGAGGTAATTAAATTTAGCCTTAGAAAAGTCT
#ENSPPAT00000042227.1/1-4166 ------------------------------------------------------------
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#
#ENSGGOT00000004685.3/1-4167 -------------------------------------------GCTTTATTTTTCAATTC
#ENST00000265113.9/1-3919    ------------------------------------------------------------
#ENSNLET00000035538.1/1-4464 TTGAAAAGTTCTAATCAAGACTTTCAGCTTTCCGATTGTGTGCGCTTTATTTTTCAATTC
#ENSPTRT00000078486.1/1-4419 TTGAAACGTTCTAATCAAGACTTTCAGCTTTCCGATTGTGTGCGCTTTATTTTTCAATTC
#ENSPPAT00000042227.1/1-4166 -------------------------------------------GCTTTATTTTTCAATTC
#ENSPPYT00000017907.3/1-3902 ------------------------------------------------------------
#"""


#if seqnumber was wrong, write it manually
#seqnumber = ""
#seqint = int(seqnumber)
#line_count = alignment_data.count('\n')
#x= linecount/(seqint + 1)

# Define the number of sequences
space_count = alignment_data.count('ENS')
line_count = alignment_data.count('\n')
x= line_count-space_count
seqnumber = space_count/x
a = int(seqnumber)
print("seqnumber, as number of sequences is ",seqnumber, "... NOT: If seqnumber was wrong, write it manually and control the output")






# Splitting the alignment data by newline character
lines = alignment_data.split('\n')

# Extracting data after spaces in each line
extracted_data = [line.split(' ')[-1] for line in lines]

# Joining the extracted data back to a string
result = '\n'.join(extracted_data)

# Split the new alignment data by lines

alignment_lines = result.strip().split('\n')



#Define the indices

indices = [(i*a+i)+c for c in range(0, a)for i in range(0, x) ]
print(indices)

# Initialize an empty list to store the combined strings
combined_strings = []

# Iterate over the indices
for index in indices:
    # Extract the string at the specified index and append it to the list
    combined_strings.append(alignment_lines[index])

# Combine the strings into one
result = ''.join(combined_strings)

# Print the result
print(result)

character_count = len(result)
print("Character count:", character_count)
int_count = int(character_count)
seq_lenght = int(int_count/a)
print("Sequence lenght:", seq_lenght)

def add_comma_every_char(text):
    resultmod = ''
    for char in text:
        resultmod += char + '", "'
    return resultmod

def modify_text(text):
    # Add 'x' to the first character
    text = '"' + text
    return text

def add_x_every_60(text):
    resultmod2 = ''
    for i in range(0, len(text), seq_lenght*5):
        resultmod2 += text[i:i+seq_lenght*5] + '\n'
    return resultmod2

def removelast(text):
    resultmod3 = text[:-5]
    return resultmod3








modified_text = add_comma_every_char(result)
init_text = modify_text(modified_text)
modified_text12 = add_x_every_60(init_text)
modified_text2 = removelast(modified_text12)
print(modified_text2)

f = open("Alignmentmatric.txt", "w")
f.write(modified_text2)
f.close()