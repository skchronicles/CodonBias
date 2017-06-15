from __future__ import print_function
import random


class GenerateRandom(object):
    def __init__(self, seq_len=999):

        self.sequence = ''
        self.seq_len = seq_len
        self.nt_composition = {  # Note that this takes the form of a cumulative distribution
            'A': 0.29,  # 0.29
            'C': 0.5,   # 0.21
            'G': 0.71,  # 0.21
            'T': 1      # 0.29
        }
        self.protein_sequence = ''

    def generate_sequence(self):
        """Method to Generate our random sequences based off of Nucleotide Composition"""
        self.sequence = 'ATG'
        for i in range(self.seq_len-6):  # Length minus 6 bc we are adding start and stop codons
            rand_num = random.random()   # returns a random number from 0.0 to 1.0
            for nt, probability in sorted(self.nt_composition.items()):
                if rand_num < probability:
                    self.sequence += nt
                    break
        self.sequence += str(random.choice(["TGA","TAA","TAG"]))
        return self.sequence

    def translate_dna(self, sequence):
        """Method to convert the DNA Sequence into an Amino Acid Sequence"""
        codons = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
        for base in range(0, self.seq_len, 3):
            if sequence[base:base+3] in codons:
                self.protein_sequence += codons[sequence[base:base+3]]
        return self.protein_sequence

    def codon_counter(self,dnaseq):
        """create a dictionary to count the occcurence of each codon in the dna string"""
        pass


def main():
    ###############################################################
    """General Usage"""
    # sequenceLength = 33                          # by default this is 999
    # seq2 = GenerateRandom(sequenceLength)        # create object
    # sequence = seq2.generate_sequence()          # creates dna sequence
    # print(sequence)
    # prot = GenerateRandom().translate_dna(sequence) # creates protein sequence
    # print(prot)
    ###############################################################
    """Creating a List containing 500 random sequences"""
    # NumberofSeq = 500
    # seqs = [GenerateRandom() for n in range(NumberofSeq)]
    # generateSeqs = [g.generate_sequence() for g in seqs]
    # print(generateSeqs)
    ################################################################
    """Creating Two Lists (proteins, DNA) containg a specified number of sequences"""
    NumberofSeq = 6000  # CHANGE THIS TO 6000 LATER
    seqobjectlist = []
    protobjectlist = []
    for obj in range(NumberofSeq):
        seqobjectlist.append(GenerateRandom().generate_sequence())
    for dnaobj in seqobjectlist:
        protobjectlist.append(GenerateRandom().translate_dna(dnaobj))
    # print(seqobjectlist)
    # print(protobjectlist)
    #################################################################
    """Basic Dictionary Append Example"""
    # items = "Whats the easiest way to add the list items to a dictionary "
    # stats = {}
    # for i in items:
    #     if i in stats:
    #         stats[i] += 1
    #     else:
    #         stats[i] = 1
    # print(stats)
    ##################################################################
    """Counting the occurences of each Amino Acid in all of the randomly generated protein sequences"""
    aminoAcid_count_dict = {}
    for protseq in protobjectlist:
        for acid in protseq:
            if acid in aminoAcid_count_dict:
                aminoAcid_count_dict[acid] += 1
            else:
                aminoAcid_count_dict[acid] = 1
    # print(aminoAcid_count_dict)
    ###################################################################
    """Counting the occurences of each codon triplet in all of the randomly generated DNA sequences"""
    codon_count_dict = {}
    for dnaseq in seqobjectlist:
        for codonpostion in range(0,len(dnaseq),3):
            if dnaseq[codonpostion:codonpostion+3] in codon_count_dict:
                codon_count_dict[dnaseq[codonpostion:codonpostion+3]] += 1
            else:
                codon_count_dict[dnaseq[codonpostion:codonpostion + 3]] = 1
    # print(codon_count_dict)
    ####################################################################
    """Finding the frequency of each Amino Acid in all of the randomly generated protein sequences"""
    aminoAcid_frquencies_dict = {}
    for k,v in aminoAcid_count_dict.items():
        frequeny = round((v/float(NumberofSeq*333)),2)  # (MULTIPLY NumberofSeq * (seq_len/3))
        aminoAcid_frquencies_dict[k] = frequeny
    # print(aminoAcid_frquencies_dict)
    freqcount = 0
    for v in aminoAcid_frquencies_dict.values():
        freqcount += v
    # print("All the indivdual amino acid frequencies add up to " + str(freqcount))
    ####################################################################
    """Finding the frequency of each codon triplet in all of the randomly generated DNA sequences"""
    codon_frequencies_dict = {}
    for k,v in codon_count_dict.items():
        frequeny = round((v/float(NumberofSeq*333)), 2)  # (MULTIPLY NumberofSeq * (seq_len/3))
        codon_frequencies_dict[k] = frequeny
    # print(codon_frequencies_dict)
    freqcount2=0
    for v in codon_frequencies_dict.values():
        freqcount2 += v
    # print("All the codon triplet frequencies add up to " + str(freqcount2))
    ####################################################################
    """Finding the fraction of each codon encoding for the same Amino Acid"""
    """******Printing this off later!*******"""
    codon_sameAA_frequency_dict = {}
    aa2codon_dict = {
        'I': ('ATT', 'ATC', 'ATA'),
        'L': ('CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG'),
        'F': ('TTT', 'TTC'),
        'M': ('ATG'),
        'C': ('TGT', 'TGC'),
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'Y': ('TAT', 'TAC'),
        'W': ('TGG'),
        'Q': ('CAA', 'CAG'),
        'N': ('AAT', 'AAC'),
        'H': ('CAT', 'CAC'),
        'E': ('GAA', 'GAG'),
        'D': ('GAT', 'GAC'),
        'K': ('AAA', 'AAG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        '*': ('TAA', 'TAG', 'TGA')}
    for k, v in aa2codon_dict.items():
        if v == "ATG" or v == "TGG":
            codonAAfrequency = round((codon_count_dict[v] / float(aminoAcid_count_dict[k])), 2)
            codon_sameAA_frequency_dict[v] = codonAAfrequency
        else:
            for codon in v:
                codonAAfrequency = round((codon_count_dict[codon] / float(aminoAcid_count_dict[k])), 2)
                codon_sameAA_frequency_dict[codon] = codonAAfrequency
    print(codon_sameAA_frequency_dict)
    ##########
    # Testing: Works!
    # if "TTT" and "TTC" in codon_sameAA_frequency_dict:
    #    print("TTT " + str(codon_sameAA_frequency_dict["TTT"]))
    #    print("TTC " + str(codon_sameAA_frequency_dict["TTC"]))
    ####################################################################
    """Finding the frequency of each codon per 1000 kilobases"""
    """******Printing this off later!*******"""
    codon_occurence_perKB_dict = {}
    for k, v in codon_count_dict.items():
        frequencyPerKB = round((3*v)/float(NumberofSeq), 1)  # 3 (bc total needs to equal 1000, not 333)
        codon_occurence_perKB_dict[k] = frequencyPerKB            # NumberofSeq = 6000
    print(codon_occurence_perKB_dict)
    ##########
    # Testing: Works! # fcount = 999.9 (not 1000 due compounding rounding)
    # fcount = 0
    # for v in codon_occurence_perKB_dict.values():
    #    fcount += v
    # print("total: " + str(fcount))
    ####################################################################
    """Printing the resulting 6000 Randomly Generated Proteins Sequences to a FASTA file"""
    """******Outputing this off!*******"""
    file_name = "skylerKuhn_translated_sequences.txt"
    with open(file_name, "w") as fastafile:
        for num, sequence in enumerate(protobjectlist):
            fastafile.write(">AA sequence " + str(num)+ "\n")
            fastafile.write(str(sequence) + "\n")
    fastafile.close()
    ####################################################################
    """Comparing and the frequencies of this Program to the Actual frequencies from PCC1720_codons.txt"""
    """******Outputing this off!*******"""
    import re
    CODON_TABlE_RE = re.compile("(\w+)\s+([\w?|\*?]\s+\d.\d+\s+\d+.\d)\s+\d+")  # (UGA) (* 0.21  0.7)  1345
    codonTableDict = {}
    file_name2 = "PCC1720_codons.txt"
    with open(file_name2) as infile:
        for line in infile:
            linelist = line.strip().replace("(", "").split(")")
            for element in linelist:
                codonline = element.strip()
                codon_match = CODON_TABlE_RE.search(codonline)
                if codon_match:
                    codonValueList = CODON_TABlE_RE.search(codonline).group(2).split()
                    codonKey = CODON_TABlE_RE.search(codonline).group(1).replace("U", "T")  # Convert RNA to DNA
                    codonTableDict[codonKey] = codonValueList
    print(codonTableDict)
    infile.close()
    file_name3 = "skylerKuhn_bias_output.txt"
    with open(file_name3, "w") as bias_out:

        for k,v in codon_occurence_perKB_dict.items():
            actual_list = codonTableDict[k]
            actual_freq = actual_list[2]
            diff = round(float(codon_occurence_perKB_dict[k]) - float(actual_freq), 1)
            bias_out.write("Simulated freq. for " + str(k) + " (" + str(actual_list[0]) + "): " + \
                           str(codon_occurence_perKB_dict[k]) + " per 1000 \t\t" + "Actual Frequency: " + \
                           str(actual_freq) + "\t" + "Diff: " + str(diff) + "\n")
    bias_out.close()
    with open(file_name3, "a") as bias_out2:
        bias_out2.write("\n")
        for k, v in codon_occurence_perKB_dict.items():
            actual_list = codonTableDict[k]
            actual_fraction = actual_list[1]
            diff_fract = round(float(codon_sameAA_frequency_dict[k]) - float(actual_fraction), 2)
            bias_out2.write("Simulated fraction for " + str(k) + " (" + str(actual_list[0]) + "): " + \
                           str(codon_sameAA_frequency_dict[k]) + " \t\t" + "Actual Fraction: " + \
                           str(actual_fraction) + "\t" + "Diff: " + str(diff_fract) + "\n")
    bias_out2.close()
    ####################################################################
    # That's all folks! - Skyler A. Kuhn

if __name__ == "__main__":
    main()