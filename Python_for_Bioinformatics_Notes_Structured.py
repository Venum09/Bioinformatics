# Python for Bioinformatics - Comprehensive Notes with Examples

# ======================================
# Section 1: Introduction to Python
# ======================================

# Python is a versatile programming language widely used in bioinformatics for data analysis, visualization, and automating biological computations.

# Example: Print a simple message
print("Welcome to Python for Bioinformatics!")

# ======================================
# Section 2: Variables and Data Types
# ======================================

# Variables store data, and Python supports different data types like integers, floats, strings, and lists.

# Example 1: Basic variable assignment
gene_name = "BRCA1"
mutation_count = 23

# Example 2: Using different data types
gc_content = 0.59  # Float for GC content
is_expressed = True  # Boolean for gene expression status

# "f"-strings allow you to embed variables and expressions directly inside strings using `{}`.

print(f"Gene: {gene_name}, Mutations: {mutation_count}, GC Content: {gc_content}, Expressed: {is_expressed}")

# ======================================
# Section 3: Control Structures
# ======================================

# Control structures like loops and conditionals help in decision-making and iterating over data.

# Example 1: If-Else for decision-making
sequence_length = 500
if sequence_length > 1000:
    print("Long sequence")
else:
    print("Short sequence")

# Example 2: For loop to iterate over a list
genes = ["BRCA1", "TP53", "EGFR"]
for gene in genes:
    print(gene)

# Example 3: While loop
counter = 0
while counter < 3:
    print("Processing...")
    counter += 1

# ======================================
# Section 4: Functions
# ======================================

# Functions encapsulate reusable blocks of code.

# Example: Calculate GC content
def calculate_gc(sequence):
    g = sequence.count("G")
    c = sequence.count("C")
    return (g + c) / len(sequence)

dna_sequence = "ATGCGCGATCGT"
gc = calculate_gc(dna_sequence)
print(f"GC Content: {gc:.2f}")

# ======================================
# Section 5: File Handling
# ======================================

# File handling is essential for working with biological datasets.

# Example: Read a FASTA file
with open("example.fasta", "r") as file:
    for line in file:
        print(line.strip())

# Example: Writing to a file
with open("output.txt", "w") as file:
    file.write("Analysis Complete\n")

# ======================================
# Section 6: Libraries for Bioinformatics
# ======================================

# Python has powerful libraries for bioinformatics like Biopython.

# Example: Parse a FASTA file with Biopython
from Bio import SeqIO

for record in SeqIO.parse("example.fasta", "fasta"):
    print(f"ID: {record.id}, Sequence: {record.seq}")

# Example: Perform pairwise alignment
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
alignment = aligner.align("ATCG", "ATG")
print(alignment)

# ======================================
# Section 7: Data Analysis with Pandas
# ======================================

# Pandas is used for handling structured data like gene expression matrices.

import pandas as pd

# Example: Create and manipulate a DataFrame
data = {
    "Gene": ["BRCA1", "TP53", "EGFR"],
    "Expression": [10.5, 8.3, 15.2]
}
df = pd.DataFrame(data)
print(df)

# Example: Filter genes based on expression
high_expression = df[df["Expression"] > 10]
print(high_expression)

# ======================================
# Section 8: Visualization with Matplotlib
# ======================================

# Visualization helps in interpreting biological data.

import matplotlib.pyplot as plt

# Example: Plot gene expression levels
genes = ["BRCA1", "TP53", "EGFR"]
expression = [10.5, 8.3, 15.2]

plt.bar(genes, expression)
plt.title("Gene Expression Levels")
plt.xlabel("Genes")
plt.ylabel("Expression")
plt.show()

# ======================================
# Section 9: Advanced Topics
# ======================================

# Explore advanced topics like machine learning for genomics.

# Example: K-means clustering using Scikit-learn
from sklearn.cluster import KMeans
import numpy as np

# Simulated gene expression data
data = np.array([[1.1, 2.2], [1.3, 2.1], [10.1, 10.2], [10.3, 10.5]])
kmeans = KMeans(n_clusters=2)
kmeans.fit(data)

print(f"Cluster centers: {kmeans.cluster_centers_}")
