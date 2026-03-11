"""
This is an advanced project which purpose is to detect in a sequence a start codon and from that point, look for the biggest ORF and translate it
However, this is not as advanced as it could be because do not take into account the mRNA splicing and other kinds of processings because it would be 
tough to predict which splicing would occur
"""

#INPUTS
#Receive the nucleic acid and check if it is misspelled
while True:
    nucleic_acid = input("Please, enter a nucleic acid (DNA or RNA): ").upper().replace(" ", "") #Format the input
    if nucleic_acid not in ("DNA", "RNA"):
        print("Not valid input.")
    else:
        break

#Receive the sequence wich is going to be analyzed
while True:
    seq = input(f"Enter a {nucleic_acid} sequence:").upper().replace(" ", "") #Clean the sequence deleting all the spaces which the sequence may contain
    
    #Check if there are a confusion with the nucleic acid selected
    if nucleic_acid == "DNA":
        if "U" not in seq:
            break
        else: 
            print("A DNA sequence must not contain U")
    elif nucleic_acid == "RNA":
        if "T" not in nucleic_acid:
            break
        else: 
            print("A RNA sequence must not contain T")

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
        "Met": ("AUG",),
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
        "Trp": ("UGG",),
        "Cys": ("UGU", "UGC"),
        "Arg": ("CGU", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "Gly": ("GGU", "GGC", "GGA", "GGG")
    }
    
    polypeptide_seq = [] #Amino acids are being stored on a list in order to format the output easily

    longer_seq = ""
    
    #Search for the longer codifing sequence
    for base_start in range(len(mrna_seq)):
        newer_seq = ""
        if mrna_seq[base_start : base_start + 3] == "AUG":
            for base_stop in range(base_start + 3, len(mrna_seq), 3):
                if mrna_seq[base_stop : base_stop + 3] in amino_acids["Stop"]:
                    newer_seq = mrna_seq[base_start: base_stop + 3] #Assign as the new sequence the fragment between a start codon and the stop codon
                    break #Stop searching for further stop codons
            if newer_seq == "": #In case a stop codon is not found
                newer_seq = mrna_seq[base_start:]
            if len(newer_seq) > len(longer_seq):
                longer_seq = newer_seq
    if longer_seq == "": #If no ORF found, take mRNA_seq as the correct ORF
        longer_seq = mrna_seq
    
    for i in range(0, len(longer_seq), 3):
        codon = longer_seq[i : i + 3]
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
#Print the main information regardless which nucleic acid is
print(50*"-") #Printing a line to make the output tidier and prettier
print(f"Given sequence:    5' {spaced(seq)} 3'")
print(f"Inverted sequence: 3' {spaced(seq[::-1])} 5'")
print(f"Length: {len(seq)} bp")
print(50*"-")

if nucleic_acid == "DNA":
    #Print the dsDNA
    print("DNA:")
    print(f"5' {spaced(seq)} 3'")
    print(f"3' {spaced(complementary_seq(seq))} 5'")
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
    
elif nucleic_acid == "RNA":
    #Print the dsDNA
    print("DNA:")
    print(f"5' {spaced(dna(seq))} 3'")
    print(f"3' {spaced(complementary_seq(dna(seq)))} 5'")
    print(50*"-")

    #Print the corresponding mRNA
    print("mRNA:")
    print(f"5' {spaced(seq)} 3'")
    print(50*"-")

    #Print the polypeptide information
    print("Polypeptide sequence:")
    polypeptide_seq = polypeptide(seq) #Store the sequence in a variable so as not to call the function twice
    print(f"NH2 - ... - {polypeptide_seq} - ... - COOH")
    print(f"Polypeptide length: {len(polypeptide_seq.split(" - ")) - polypeptide_seq.count("Stop")} amino acids") #It must split again in order to count the amino acids
    print(50*"-")
