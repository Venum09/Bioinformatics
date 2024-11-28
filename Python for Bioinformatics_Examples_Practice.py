#Logical Operations
T=25
A=25
C=40
G=40

T=25

T==A
#Strings
Seq1 = 'ATGCATGCATATAATATATATGC'
Seq2 = 'TATAGCGCUTTAA'
print(Seq1)
print(Seq2)

#Converting sequences in Upper/Lower case

Seq1.lower()

Seq2.lower()


#Count the length of the sequence

len(Seq1)

'A' in Seq1
'U' in Seq1 + Seq2

# Count characters

Seq1.count('A')
A_count=Seq1.count('A')
T_count=Seq1.count('T')
G_count=Seq1.count('G')
C_count=Seq1.count('C')
 
print(A_count)
print(T_count)
print(C_count)


# Now get GC context

GC=G_count + C_count
length=len(Seq1)
GC_Percentage=GC/length*100
print(GC_Percentage)

# Combine sequences

WS=Seq1 +' '+ Seq2.lower()
print(WS)

# Add " " to add space or "%s %s"%

#round function to roundoff at certain decimals

# GC method 2
Seq1_length=len(Seq1)
GC=Seq1.count('G')+Seq1.count('C')
GC_percentage=GC/Seq1_length*100
print('GC % is = ', GC_percentage)
GC_rounded=round(GC_percentage,2)
print('GC % is = ', GC_rounded)


# Number of Purines and Pyrimidines
#variable.()count format to count data in variable
purines=Seq1.count('A') + Seq1.count('G')
print('Number of Purines is', purines)

pyrmdn=Seq1.count('C') + Seq1.count('T')
print('Number of Pyrimidines is',pyrmdn)

percentage_purines=purines/Seq1_length*100
print('Percentage of Purines = ', percentage_purines)

#Extracting Sub-sequences using Indexing and Slicing

# Indexing starts form '0; mathematically (n-1)

# Example
print('Nucleotide number form  3 to 5 is ', Seq1[8])

       # To start form the right end use negative integers; it starts from '-1'

print('Nucleotide number from 5 to 3 is', Seq1[-7])


# SLicing Example

Seq1[0:-1]

# Adding third Index

#Format Variable[m:n:s]; m = start, n=length covered, s=interval

print('Third index is ', Seq1[0:8:2])
Seq1[0:9:3]

# Short hand Indexing

# Exapmple | variable[:n] slice from the begining just before n, varible[m:] starts slicing from index m to the end of the sequences, variable[:] returns the copy of the entire sequence

Seq1[::]
Seq1[:5]
Seq1[2:]


# Data Type-Lists

list1=['A','T','G','C']
subseq=['TTA','CGCG','GAAT']

# Convert String to the list

string='ATGCTAGG'
stl=list(string)
print(stl)

noi=len(subseq)
print(noi)

# Adding New Element

list1.append("TGA")
print(list1)

# Indexing

# Example- First element of the list

list1[0:5]

       # Sort a list

list1.sort()
print(list1)

       # Reverse the list

list1.reverse()
print(list1)

# Insertion and Deletion

list1.insert(0, "TTT")
print(list1)

del list1[0]
print(list1)

       # Delete item value
list1.remove("T")
print(list1)

       # Change nucleotides

list1[0]= 'AT'
print(list1)

       # Convert list to string

Seq="".join(list1)
print(Seq)


# Reverse Sequqences

    # When sequence is a string

print(Seq1)
Seq1[::-1]

rSeq1=Seq1[::-1]
print(rSeq1)

    # When sequence is list

# Approach1

forSeq=''.join(list1)
print(forSeq)
# Reverse the sequence
revSeq=forSeq[::-1]
print(revSeq)

#Approach2 

revSeq=''.join(revSeq)
print(revSeq)

# Compare

print(''.join(list1))
print(revSeq)

# Convert Reverse list to string

Scopy=['A','T','T','G']
Scopy.extend(list1)

print(len(Scopy),Scopy)


Scopy.reverse()
print(Scopy)

rSeq=''.join(Scopy)
print('total is',rSeq)

print(len(rSeq))
print(len(list1))
print(len(Scopy))






# Repeat Programming Tasks using functions

    # 1st word 'def' then parenthesis; add parameters in parenthesis() to add new values and add colon(:). Add return value to use output outside the function

def GC_percentage():
    Seq1_length=len(Seq1)
    GC=Seq1.count('G')+Seq1.count('C')
    GC_percentage=GC/Seq1_length*100
    print('GC % is = ', GC_percentage)
    GC_rounded=round(GC_percentage,2)
    print('GC % is = ', GC_rounded)

# Calling the function


GC_percentage()

# Reusing the values
def GC_percentage():
    Seq1_length=len(Seq1)
    GC=Seq1.count('G')+Seq1.count('C')
    GC_percentage=GC/Seq1_length*100
    print('GC % is = ', GC_percentage)
    GC_rounded=round(GC_percentage,2)
    print('GC % is = ', GC_rounded)
    return GC_rounded

GC_value=GC_percentage()
print(GC_value)


# Modifying the functions

def GC_percentage(sequence):
    Seq1=(sequence)
    GC=Seq1.count('G')+Seq1.count('C')
    GC_percentage=GC/Seq1_length*100
    GC_rounded=round(GC_percentage,2)
    print('GC % is = ', GC_rounded)
    return GC_rounded

S1='ATGTCGTAG'
S2='GCTGAGTAG'

# Utilising existing function (GC_percentage)

GC1=GC_percentage(S1)
GC2=GC_percentage(S2)

# Counting features in multiple genebank files

# 1 Get Gene Bank Files
# 2 Make the same directory
# 3 Count features using python
# 4 Save results for individual files

import glob #Library used to find files that match pattern
from Bio import SeqIO # Biopython lib. for reading and writing sequences file format
from collections import Counter # Specialised dictionary for counting hashable objects
import pandas as pd # Data manipulation and analysis library
import os # Used for interacting with OS, e.g reading & writing

# Assigns the directory path where your gene bank files are located to the variable file_dir.

file_dir=basedir="c:/Users/Chintansinh Rathod/Desktop/Python/GBF"

#Uses the glob function to search within file_dir for all files ending with .gb and stores the list of these file paths in gfiles.

gfiles=glob.glob("%s/*.gb"%file_dir)

gfiles
print(gfiles)
print(len(gfiles))


'C:/Users/Chintansinh Rathod/Desktop/Python/GBF/S1.gb'

def count_features(gfile):
    genebank_object=SeqIO.read(gfile,"gb")
    features=genebank_object.features
    feature_types=[feature.type for feature in features]
    feature_count=Counter(feature_types)
    print('features have been counted')

    dataframe=pd.DataFrame(feature_count.items(), columns=['feature','count'])

    directory,filename=os.path.split(gfile)
    filename=filename.strip('.gb')

    basedit='c:/Users/Chintansinh Rathod/Desktop/Python/GBF'

    outputfile='%s/%s.csv'%(basedir,filename)
    dataframe.to_csv(outputfile,index=False)

    print('Count data has been saved')

for gfile in gfiles:
    count_features(gfile)


# Count features in multiple Gene Bank Data


# import libraries

import glob #Library used to find files that match pattern
from Bio import SeqIO # Biopython lib. for reading and writing sequences file format
from collections import Counter # Specialised dictionary for counting hashable objects
import pandas as pd # Data manipulation and analysis library
import os # Used for interacting with OS, e.g reading & writing

# Set a path

file_dir='c:/Users/Chintansinh Rathod/Desktop/Python/GBF'

# Get file names

gfiles=glob.glob('%s/*.gb'%file_dir)
print(gfiles)

gfiles[0]


# Function to read genebank files and get the features

def read_file(gfile):
    genebank_object=SeqIO.read(gfile,"gb")
    features=genebank_object.features
    feature_types=[feature.type for feature in features]
    return feature_types

# Count features of the genebank files

def count_features(feature_types):
    feature_count=Counter(feature_types)
    print("Features have been counted")
    return feature_count

# Function to identify features in all the files

def scan_all_features(files):
    allfeatures=[]
    for gfile in gfiles:
        feature_types=read_file(gfile)
        allfeatures.extend(feature_types)

    allfeatures=set(allfeatures)
    allfeatures=list(allfeatures)
    print('All features have been identified')
    return allfeatures


# Get all the features

allfeatures=scan_all_features(gfiles)


print(allfeatures)

# Create empty list to count features

allfeature_count=[]

# Use for loop to iterate over the files and count features

for gfile in gfiles:
    directory,filename=os.path.split(gfile)
    filename=filename.strip('.gb')
    feature_types=read_file(gfile)
    feature_count=count_features(feature_types)
    temp_count=[]

    temp_count.append(filename)

    for feature in allfeatures:
        if feature in feature_count.keys():
            temp_count.append(feature_count[feature])
        else:
            temp_count.append(0)
    allfeature_count.append(temp_count)

    print(len(allfeature_count))

print(allfeature_count)


# Columns for data frame

columns=[]
columns.append('File')
columns.extend(allfeatures)

# Data frame using pandas

dataframe=pd.DataFrame(allfeature_count,columns=columns)
print(dataframe)



# Finding Reverse Sequence

Seq1[::-1]
reverse_seq1=Seq1[::-1]
print(Seq1)
print(reverse_seq1)


# Dealing with List data

# Approach 1- Convert to string and find reverse using slicing

f_seq= ''.join(list1)
print(f_seq)

# Approach 2- List slicing and convert to string
print(list1)
list_1=''.join(list1)
print(list1)
rlist1=list_1[::-1]

#Regular example to try
A=['A','C','G','T','A']
B=''.join(A)
print(B)
rB=B[::-1]
print(rB)


# Approach 3

eSeq = ['A','T','G','G','A','T','G','G','A','T','G','G']

exSeq = ['A','T','G','G','A','T','G','G','A','T','G','G']

print(exSeq)

exSeqS=''.join(exSeq)
print(exSeqS)

exSeqSR=exSeqS[::-1]
print(exSeqSR)




# Generate own DNA Sequences

length=10
#Supply the list of nucletides

Ntds=['A','T','G','C']

# Importing Random Module

from random import choice

# Call empty sequence

seqE=''

# Randomly create a sequence by adding bases

for i in range(length):
    base=choice(Ntds)
    seqE+=base

print(seqE)

# Approach 2 - List Comprehension

seqE=[choice(Ntds) for i in range (length)]
print(seqE)

SeqE=''.join(seqE)
print(SeqE)



#Genrerating Multiple Sequences

    # Import random module

def genrate_seq(length=10):
    bases=['A','T','G','C']
    seq=[choice(bases) for i in range(length)]
    seq=''.join(seq)
    return seq

    # Generate sequences of the same length Scenario 1

num_of_seq=5
seq=[genrate_seq(5) for i in range(num_of_seq)]
print(seq)

    # Generate sequences of the same length Scenario 2: Different length

lengths=[5,10,15,20]
seq=[genrate_seq(length) for length in lengths]
print(seq)

len(seq[0])


# Random seq generations

from random import choice
from random import randint

# Define the number of sequences

num_of_seq=5


lengths=[randint(5,10) for i in range(num_of_seq)]
print(lengths)

seq=[genrate_seq(length) for length in lengths]
print(seq)




# How to find the complements of nucleotides sequences

Mseq='GTGATGACGGCGGTATCGGC'
print(Mseq)

len(Mseq)

Bdicts={'A':'T','G':'C','T':'A','C':'G'}
print(Bdicts)

# Find complements of the nuclotides apporach 1

seqcomp=''

for base in Mseq:
    Bdict=Bdicts[base]
    seqcomp+=Bdict

print(' ',Mseq)
print(' ',seqcomp)

# Approach 2: Using list comprehensions

seqcomp1=[Bdicts[base] for base in Mseq]
seqcomp1=''.join(seqcomp1)
print(seqcomp1)



# Find complement of multiple sequences 

Dseq=['ATTCA', 'TATACGCAGG', 'CCTTATAGCAGAGCA', 'TAAAGTTTGATTTACGCTAC']
len(Dseq)

def scf (sequence):
    Bdicts={'A':'T','G':'C','T':'A','C':'G'}
    sc=[Bdicts[base] for base in sequence]
    sc=''.join(sc)
    return sc


# Use for loop to iterate over the list of sequences and their respoective complementry seqs.
print('')
 
for sequence in Dseq:
    sc=scf(sequence)
    Dseq=''.join(Dseq)
    
    print('Original Sequence    :    ',sequence)
    print('Complementry Seqence :    ',sc)


# Reverse Complement of Multiple Sequences

def reverse_sequence(sequence):
    reverse=sequence[::-1]
    rev_comp=scf(Dseq)
    return rev_comp


# Call the sequence to find the reverse of the sequence

for sequence in Dseq:
    revS=reverse_sequence(sequence)

print('Forward:    ',Dseq)
print('Reverse:    ',revS)



    
