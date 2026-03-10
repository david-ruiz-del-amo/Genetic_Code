"""
This is a basic project which principal purpose is to translate a sequence in spite of not containing start codons or stop codons
"""
#INPUTS
#Receive the sequence
seq = input(f"Enter a sequence:").upper().replace(" ", "") #Clean the sequence deleting all the spaces which the sequence may contain    

#FUNCTIONS
#Create a function which returns the complementary strand of DNA
def complementary_seq(seq):
    complementary_strand = ""
    for base in seq:
        if base == "A":
            complementary_strand += "T"
        elif base == "T":
            complementary_strand += "A"
        elif base == "G":
            complementary_strand += "C"
        elif base == "C":
            complementary_strand += "G"
        else:
            complementary_strand += base #In case there is an indetermination or a misspelled base
    
    return complementary_strand

#Create a function which returns a mRNA sequence from a DNA sequence (it will make a more readable code at the output)
def mrna(seq):
    return seq.replace("T", "U")

#Create a function which returns a DNA sequence from a mRNA sequence
def dna(seq):
    return seq.replace("U", "T")

#Create a function which returns the polypeptide codified by the sequence
def polypeptide(mrna_seq):
    #A dictionary is used to make a relation between the amino acid (key) and the codon (value)
    amino_acids = {
        "Phe": ("UUU", "UUC"),
        "Leu": ("UUA", "UUG", "CUU", "CUC", "CUA", "CUG"),
        "Ile": ("AUU", "AUC", "AUA"),
        "Met": ("AUG"),
        "Val": ("GUU", "GUC", "GUA", "GUG"),
        "Ser": ("UCU", "UCC", "UCA", "UCG", "AGU", "AGC"),
        "Pro": ("CCU", "CCC", "CCA", "CCG"),
        "Thr": ("ACU", "ACC", "ACA", "ACG"),
        "Ala": ("GCU", "GCC", "GCA", "GCG"),
        "Tyr": ("UAU", "UAC"),
        "Stop": ("UAA", "UAG", "UGA"),
        "His": ("CAU", "CAC"),
        "Gln": ("CAA", "CAG"),
        "Asn": ("AAU", "AAC"),
        "Lys": ("AAA", "AAG"),
        "Asp": ("GAU", "GAC"),
        "Glu": ("GAA", "GAG"),
        "Trp": ("UGG"),
        "Cys": ("UGU", "UGC"),
        "Arg": ("CGU", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "Gly": ("GGU", "GGC", "GGA", "GGG")
    }
    
    polypeptide_seq = [] #Amino acids are being stored on a list in order to format the output easily
    for i in range(0, len(mrna_seq), 3):
        codon = mrna_seq[i : i + 3]
        exists = None #It is used to detect unknown codons (e.g. CUN)
        for key, value in amino_acids.items():
            if codon in value:
                polypeptide_seq.append(key)
                exists = True
                break #Escape the loop so as to go faster
        if exists == None:
            polypeptide_seq.append("[?]")
        if polypeptide_seq[-1] == "Stop":
            break #Stop transalting
    
    return " - ".join(polypeptide_seq) #Returns a decorated string

#Create a function which will be used to format the outputs
def spaced(seq):
    codon_list = []

    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        codon_list.append(codon)

    return " ".join(codon_list) #Returns a formatted string which has the same sequence, but it is separated to read better the sequence

#OUTPUT
"""
Regardless of the acid nucleic, the output are merged. 
It may be redundant or inefficient, but to complete its purpose 
(tiny tasks like highschool/university homework) it is enough.
"""
#Print the main information regardless which nucleic acid is
print(50*"-") #Printing a line to make the output tidier and prettier
print(f"Given sequence:    5' {spaced(seq)} 3'")
print(f"Inverted sequence: 3' {spaced(seq[::-1])} 5'")
print(f"Length: {len(seq)} bp")
print(50*"-")

#Print the dsDNA
print("DNA:")
print(f"5' {spaced(dna(seq))} 3'")
print(f"3' {spaced(complementary_seq(dna(seq)))} 5'")
print(50*"-")

#Print the corresponding mRNA
print("mRNA:")
print(f"5' {spaced(mrna(seq))} 3'")
print(50*"-")

#Print the polypeptide information
print("Polypeptide sequence:")
polypeptide_seq = polypeptide(mrna(seq)) #Store the sequence in a variable so as not to call the function twice
print(f"NH2 - ... - {polypeptide_seq} - ... - COOH")
print(f"Polypeptide length: {len(polypeptide_seq.split(" - ")) - polypeptide_seq.count("Stop")} amino acids") #It must split again in order to count the amino acids
print(50*"-")