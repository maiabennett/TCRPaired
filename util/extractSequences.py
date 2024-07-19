
# Based on code developed by Austin Seaman, Dario Ghersi, and Ryan Ehrlich for the TRain and SwarmTCR projects
# source: 

import os
import pickle, argparse
import sys
import csv
from venv import create 
import pandas as pd 
import numpy as np
from Bio import Align


# The following indices correspond to their respective CDR Loops
c1 = [26,38]
c2 = [55,65]
c25 = [80,86]
c3 = [104,116]

class ref_CDRAASeq():
    # intialize extract_CDRaa
    def __init__(self, file):
        self.file = file

    # getSeqs returns entire AA sequences for each TCR gene family
    def getSeqs(self, lines, indices, a, flag):
        if flag == 0:
            start, end = indices[a], indices[a + 1]
            seq = lines[start:end]
        if flag == 1:
            start = indices[a]
            seq = lines[start:]
        i2 = [i for i, x in enumerate(seq[0]) if x == '|']
        gene = seq[0][i2[0] + 1:i2[1]]
        aa = seq[1:]
        aaSeq = '%s%s' % (aa[0][:-1], aa[1][:-1])
        return(gene, aaSeq)

    # getLoops pulls uses the CDR indices above to return a AA sequence
    def getLoops(self, seq_):
        cdr1 = seq_[c1[0]:c1[1]].replace('.', '')
        cdr2 = seq_[c2[0]:c2[1]].replace('.', '')
        cdr25 = seq_[c25[0]:c25[1]].replace('.', '')
        cdr3 = seq_[c3[0]:c3[1]].replace('.', '')
        whole = seq_.replace('.', '')
        return (cdr1, cdr2, cdr25, cdr3, whole)

    # storeSeqs produces a dictionary for TRAV and TRBV genes. Dictionary keys are
    # TR genes, the values are the AA sequences for each loop (c1, c2, c25, c3).
    def storeSeqs(self):
        imgtFile = open(self.file, 'r')
        cdrDict = {}
        lines = imgtFile.readlines()
        indices = [i for i, x in enumerate(lines) if x[0] == '>']
        for a in range(len(indices)):
            if a != len(indices) - 1:
                gene, aaSeq = self.getSeqs(lines, indices, a, 0)
                cd1, cd2, cd25, cd3, whole = self.getLoops(aaSeq)
                cdrDict[gene] = []
                cdrDict[gene].append(cd1)
                cdrDict[gene].append(cd2)
                cdrDict[gene].append(cd25)
                cdrDict[gene].append(cd3)
                cdrDict[gene].append(whole)
            else:
                gene, aaSeq = self.getSeqs(lines, indices, a, 1)
                cd1, cd2, cd25, cd3, whole = self.getLoops(aaSeq)
                cdrDict[gene] = []
                cdrDict[gene].append(cd1)
                cdrDict[gene].append(cd2)
                cdrDict[gene].append(cd25)
                cdrDict[gene].append(cd3)
                cdrDict[gene].append(whole)
        return(cdrDict)

    # Write dictionary to pickle file to be used in downstream analysis
    def write_dictionary(self, dictionary, out_name):
        outfile = '%s.pickle' % (out_name)
        with open(outfile, 'wb') as handle:
            pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Write dictionary to csv file to be used in downstream analysis    
    def write_csv(self, dictionary, out_name):
        outfile = "%s.csv" % (out_name)
        with open(outfile, "w") as f:  # Open file, name is set by default but can be changed
            f.write("Gene, CDR1, CDR2, CDR2.5, CDR3, whole   \n")
            for each in dictionary:  # Loop through each clone ID
                chains = dictionary[each]
                f.write(each + "," + ",".join(chains) + "\n")


# Convert new paired sequencing data to CDR loop format
def convert_data_paired(new, gene_column_alpha, gene_column_beta):
    ref = pd.read_csv('./util/ref/paired_genes.csv')
    new['CDR1a'] = new[gene_column_alpha].map(ref.set_index('Gene')['CDR1'])
    new['CDR2a'] = new[gene_column_alpha].map(ref.set_index('Gene')['CDR2'])
    new['CDR2.5a'] = new[gene_column_alpha].map(ref.set_index('Gene')['CDR2.5'])
    new['CDR1b'] = new[gene_column_beta].map(ref.set_index('Gene')['CDR1'])
    new['CDR2b'] = new[gene_column_beta].map(ref.set_index('Gene')['CDR2'])
    new['CDR2.5b'] = new[gene_column_beta].map(ref.set_index('Gene')['CDR2.5'])
    return new
        


# Data --> full, aligned sequence methods
# Taken directly from SeqConductor.py and adapted for use with dataframes

gene_dic = {}  # Gene dictionary, global

#################
#    Methods    #
#################
def create_gene_dic(fasta_file):
    """
    Populate the gene dictionary from IMGT data.

    Parameters
    __________
    fasta_file : str
        Fasta file containing IMGT data of gene segments of TCRs
    Returns
    _______
    gene_dic : dict
        Gene segment name -- AA sequence
    """
    with open(fasta_file, "r") as family:
        temp_fam = ""  # Holds current gene segment name
        for line in family:
            if line[0] == ">":  # Find header
                header = line.split("|")
                temp_fam = header[1]  # gene segment name, creates temp header name
                gene_dic[temp_fam] = ""  # prepares location in dictionary
            # Collecting sequence data from fasta file per header
            else:
                if line.__contains__('*'):
                    gene_dic[temp_fam] += line[:-1].replace('*', "")
                else:
                    gene_dic[temp_fam] += line[:-1]
    return gene_dic


def get_tcr_info(data):
    """
    Create dictionary of components of TCR chains from table

    Parameters
    __________
    data : pandas.DataFrame
        DataFrame containing AV, AJ, BV, BJ, CDR3A, CDR3B columns

    Returns
    _______
    tcr_dic : dict
        clone ID Positions: 0 - TRAV, 1 - CDR3A, 2 - TRAJ, 3 - TRBV, 4 - CDR3B, 5 - TRBJ, clone ID - 6
    """
    tcr_dic = {}  # Creates output dictionary
    clone_id_count = {}  # Keeps track of repeat clone ids
    previous_id = ""
    positions = list(range(7))
    # if not clone_id provided
    for key, value in data.iterrows():  # Reading over table with pandas
        values = value.array
        # Finds each clone_id.
        if isinstance(values[positions[1]], str):
            # If no clone id, update with previous name or no clone id provided
            if not isinstance(values[positions[0]], str) or positions[0] == -1:
                if previous_id != "":
                    values[positions[0]] = previous_id
                else:  # When not using clone id
                    value[positions[0]] = "id-000"
            # Remove second class (ex.  TRAV12-2*01-*03 and TRBV20-1*01-*07: TRAV12-2*01 and TRBV20-1*01)
            # .split(";")[0] splits for when record have two gene segments listed, typically identical segments in aa
            if values[positions[1]].count("-") >= 2:
                values[positions[1]] = "-".join(values[positions[1]].split("-")[:-1])
            if values[positions[4]].split(";")[0].count("-") >= 2:
                values[positions[4]] = "-".join(values[positions[4]].split(";")[0].split("-")[:-1])
            # Confirm each gene segment name is in-fact a name in gene_dic
            if confirm_segment([values[positions[1]], values[positions[3]], values[positions[4]].split(";")[0],
                                values[positions[6]]]):
                if values[positions[0]] in clone_id_count.keys():  # Catch duplicate clone id names and rename "-01"
                    clone_id_count[values[positions[0]]] += 1
                    values[positions[0]] = \
                        f"{values[positions[0]]}-{str(clone_id_count[values[positions[0]]]).zfill(3)}"
                else:
                    clone_id_count[values[positions[0]].split("-")[0]] = 0
                previous_id = values[positions[0]].split("-")[0]
                tcr_dic[values[positions[0]]] = [values[positions[1]], values[positions[2]], values[positions[3]],
                                                 values[positions[4]].split(";")[0], values[positions[5]],
                                                 values[positions[6]]]
                
    return tcr_dic


def confirm_segment(seg_list):
    """
    Loop through gene segments of a clone id and return if a segment is not found in gene_dic

    Parameters
    __________
    seg_list : lst
        List of gene segments w/ clone id position 0

    Returns
    _______
    confirm segment : boolean
        True if gene segment discovered in gene_dic
    """
    for seg in seg_list:
        if seg not in gene_dic.keys():
            return False
    return True


def make_tcr_seq(tcr_parts):
    """
    Sends tcr parts into align_overlap to create resulting Alpha and Beta chain sequences.

    Parameters
    __________
    tcr_parts : dict
        dictionary containing each aa segment of TCR as constructed from get_tcr_info

    Returns
    _______
    output_segs : dict
        a dictionary with clone ID key and a list of Alpha and Beta chain sequences.
    """
    output_seqs = {}
    # For each clone-id
    for each in tcr_parts:
        # Alignment construction for TRAV, and CDR3A
        temp_seq_alpha = align_overlap(gene_dic[tcr_parts[each][0]], tcr_parts[each][1], "V", each)
        # Alignment construction for TRAV + CDR3A and TRAJ
        temp_seq_alpha = align_overlap(temp_seq_alpha, gene_dic[tcr_parts[each][2]], "J", each)
        # # Alignment construction for TRBV, and CDR3B
        temp_seq_beta = align_overlap(gene_dic[tcr_parts[each][3]], tcr_parts[each][4], "V", each)
        # Alignment construction for TRBV + CDR3B and TRBJ
        temp_seq_beta = align_overlap(temp_seq_beta, gene_dic[tcr_parts[each][5]], "J", each)
        parts = tcr_parts[each]
        # 0: Alpha_chain, 1: Beta chain, 2: aV segment, 3: CDR3a, 4: aJ segment, 5: bV segment, 6: CDR3b, 7: bJ segment
        # 8: aV seq, 9: aJ seq, 10: bV seq, 11: bJ seq
        output_seqs[each] = [temp_seq_alpha, temp_seq_beta, parts[0], parts[1], parts[2], parts[3], parts[4], parts[5],
                            gene_dic[tcr_parts[each][0]], gene_dic[tcr_parts[each][2]], gene_dic[tcr_parts[each][3]],
                            gene_dic[tcr_parts[each][5]]]
    return output_seqs


# Method: 05 - align_overlap()
# Input:
#   front - start of seq
#   end - end of seq
#   part - If overlapping V or J segment
# Goal: Looks for overlap between front and end AAs besides when single AA overlap
# Output: Resulting overlapped segments -- if == J: full TCR sequence
def align_overlap(front, end, part, clone_id):
    """
    Looks for overlap between front and end AAs besides when single AA overlap

    Parameters
    __________
    front : str
        start of seq
    end : str
        end of seq
    part : str
        "V" or "J" depending on what segments are being overlapped to the cdr
    clone_id : str
        reference to what constant region to append to end of variable region sequence

    Returns
    _______
    output : str
        Resulting overlapped segments -- if == J and clone_id: full TCR sequence
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    best_align = []
    # Finding shortest alignment with best score to append regions
    # overlap cut modified for v or j based on testing
    if part == "V":
        # Alignment Weights
        aligner.match_score = 2.5
        aligner.mismatch_score = -1
        aligner.gap_score = -1.5
        aligner.extend_gap_score = -0.5
        # Perform alignment
        # When an alignment conforms to standard of at least one of the last 3 aa's of the V segment included in align
        good = False
        if len(end) > 7:
            overlap = (len(end) * -1) + 5  # Will decrease overlap if amino acids don't align properly (add from len)
        else:
            overlap = (len(end) * -1)  # If short CDR provided
        while not good:
            temp_align = aligner.align(front[overlap:], end)
            # Reference first alignment - catch if no alignment found at all.
            try:
                align_ = temp_align[0]
            except:
                output = front + end
                print("Poor alignment possible: " + clone_id)
                return output
            # Collect alignment info
            front_seq = temp_align.alignment.format().split("\n")[0]
            pairs_info = temp_align.alignment.format().split("\n")[1]
            end_seq = temp_align.alignment.format().split("\n")[2]
            # If last match is within the range needed
            if end_seq[0] == " " and front_seq[0] != " ":
                best_align = []
                best_align.append(overlap)  # overlap
                best_align.append(temp_align.alignment.format().split("\n")[1].count(" "))  # Count start position
                good = True
            else:
                if pairs_info.count("|") < 2:  # Only a single aa match - force append CDR3 and warn
                    best_align = "BAD"
                    good = True
                elif len(front_seq) != 2:  # Bump overlap (trim front_seq) by one if not a good alignment
                    overlap += 1

    elif part == "J":
        offset = 6  # Helps with alignments to avoid too spaced out alignment
        # Alignment Weights
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.gap_score = -1
        aligner.extend_gap_score = -0.5
        # Perform alignment
        temp_align = aligner.align(front[(len(end) - offset) * -1:], end)
        try:
            align_ = temp_align[0]  # Reference first align - catch if no alignment found at all.
        except:
            output = front + end
            print("Poor alignment possible: " + clone_id)
            return output
        # Collect alignment info to print
        front_seq = temp_align.alignment.format().split("\n")[0]
        pairs_info = temp_align.alignment.format().split("\n")[1]
        end_seq = temp_align.alignment.format().split("\n")[2]
        # Collected needed information from alignment
        start_seq = temp_align.alignment.format().split("\n")[0]
        end_seq = temp_align.alignment.format().split("\n")[2]
        end_pos = len(temp_align.alignment.format().split("\n")[1])
        best_align = [start_seq, end_seq, end_pos]
    if part == "V":
        if best_align != "BAD":
            end_front = len(front) + best_align[0] + best_align[1]
            output = front[:end_front] + end
        else:
            output = front + end
        return output
    elif part == "J":  # If joining region overlap
        output = front + best_align[1][best_align[2]:]  # Overlaps front with remainder of best_align end seq
        return output


def if_gap(front, end, region, gap_count):  # Check for cdr len in previous method, send in if V or J region.
    """
    Determine if there are gaps in the alignments and where in the sequence it's located

    Parameters
    __________
    front : str
        start of seq
    end : str
        end of seq
    region : str
        "V" or "J" depending on what segments are being overlapped to the cdr
    gap_count : int
        location of gap

    Returns
    _______
    seq : str
        Adjusted aligned sequence
    """
    temp_seq = ""
    if region == "V":  # When the variable region and cdr are being overlapped.
        for letter in front:
            if letter != "-":
                temp_seq += letter
            else:
                return temp_seq + end[len(temp_seq):]
    elif region == "J":  # When the joining region and cdr are being overlapped.
        return front + end[gap_count:]


def append_constant(tcr_seq_dic, alpha, beta):
    """
    Take in the variable region sequence and append the constant region based on the organism name provided

    Parameters
    __________
    tcr_seq_dic : dict
        dictionary of tcr components clone_id: alpha, beta, [other parts]
    alpha : str
        alpha chain constant region id
    beta : str
        beta chain constant region id

    Return
    ______
    tcr_seq_dic : dict
        updated dictionary with all tcr sequences appended with constant regions
    """
    constant_dic = create_constant_dic("Homo sapiens")
    for tcr_id in tcr_seq_dic:
        for segment in constant_dic:
            if constant_dic[segment][0].startswith(alpha) and constant_dic[segment][1] == "EX1":
                tcr_seq_dic[tcr_id][0] = tcr_seq_dic[tcr_id][0] + constant_dic[segment][2].replace("X", "")
            if constant_dic[segment][0].startswith(beta) and constant_dic[segment][1] == "EX1":
                tcr_seq_dic[tcr_id][1] = tcr_seq_dic[tcr_id][1] + constant_dic[segment][2].replace("X", "")
    return tcr_seq_dic


def create_constant_dic(organism):
    """
    Populate the gene dictionary from IMGT data.

    Parameters
    __________
    organism : str
        name of organism provided by user

    Returns
    _______
    constant_dic : dict
        dictionary of the constant regions for the organism selected
    """
    constant_dic = {}
    constant_region_file = './util/ref/constant_regions.fasta' 
    with open(constant_region_file, "r") as f:
        temp_imgt_id = ""
        flag = False  # Flag to know when to collect sequence information for dic
        for line in f:
            if line[0] == ">":  # Determine if header
                if line.split("|")[2].upper() == organism.upper():  # Determine if correct organism
                    header = line.split("|")
                    temp_imgt_id = header[0][1:] + "-" + header[4]
                    segment_name = header[1]  # ex. TRAC*01
                    exon_num = header[4]  # ex. EX1
                    constant_dic[temp_imgt_id] = [segment_name, exon_num]  # prepares location in dictionary
                    flag = True
                else:
                    flag = False
            elif flag:
                # Collecting sequence data from fasta file per header
                if len(constant_dic[temp_imgt_id]) == 2:  # Catches when multiple lines of seq
                    if line.__contains__('*'):
                        constant_dic[temp_imgt_id].append(line[:-1].replace('*', ""))
                    elif line.__contains__('X'):
                        constant_dic[temp_imgt_id].append(line[:-1].replace('X', ""))
                    else:
                        constant_dic[temp_imgt_id].append(line[:-1])
                elif len(constant_dic[temp_imgt_id]) == 3:
                    if line.__contains__('*'):
                        constant_dic[temp_imgt_id][2] += (line[:-1].replace('*', ""))
                    elif line.__contains__('X'):
                        constant_dic[temp_imgt_id][2] += (line[:-1].replace('X', ""))
                    else:
                        constant_dic[temp_imgt_id][2] += (line[:-1])
    return constant_dic


def make_fasta_files(alpha_file, beta_file, tcr_dic):
    """
    Create alpha.fasta and beta.fasta - containing alpha and beta chain sequences with id as header

    Parameters
    __________
    alpha_file : str
        name of alpha chain file
    beta_file : str
        name of beta chain file
    tcr_dic : dict
        dictionary of all tcrs
    character_list : str
        list of characters to not include in fasta header id
    """
    with open(alpha_file, "w") as file_a:
        with open(beta_file, "w") as file_b:
            for clone_id in tcr_dic:
                clone_id_print = clone_id
                tcr_alpha_chain = tcr_dic[clone_id][0]
                tcr_beta_chain = tcr_dic[clone_id][1]
                count_1 = 1
                count_2 = 1
                if tcr_alpha_chain != '':
                    file_a.write('>' + clone_id_print + "\n")
                    for aa in tcr_alpha_chain:
                        if count_1 % 81 != 0:
                            file_a.write(aa)
                            count_1 += 1
                        else:
                            file_a.write('\n' + aa)
                            count_1 = 2
                    file_a.write('\n')
                if tcr_beta_chain != '':
                    file_b.write('>' + clone_id_print + '\n')
                    for aa in tcr_beta_chain:
                        if count_2 % 81 != 0:
                            file_b.write(aa)
                            count_2 += 1
                        else:
                            file_b.write('\n' + aa)
                            count_2 = 2
                    file_b.write('\n')


def return_sequences(tcr_seq_dic):
    """
    Convert the tcr_seq_dic dictionary to a pandas DataFrame.

    Parameters
    __________
    tcr_seq_dic : dict
        dictionary containing clone IDs as keys and alpha and beta sequences as values

    Returns
    _______
    df : pd.DataFrame
        DataFrame containing clone ID, alpha sequence, and beta sequence
    """
    # 0: Alpha_chain, 1: Beta_chain, 2: aV segment, 3: CDR3a, 4: aJ segment, 5: bV segment, 6: CDR3b, 7: bJ segment
    # 8: aV seq, 9: aJ seq, 10: bV seq, 11: bJ seq
    new = pd.DataFrame.from_dict(tcr_seq_dic, orient = "index") 
    new.columns = ['Alpha', 'Beta', 'AV', 'CDR3A', 'AJ', 'BV', 'CDR3B', 'BJ', 'AV_seq', 'AJ_seq', 'BV_seq', 'BJ_seq']
    new['clone.id'] = new.index
    new = new.reset_index(drop=True)
    new = new[['clone.id', 'Alpha', 'Beta']]
    return new

# all functionalities taken from SeqConductor.py
def manage_sequences(data, alpha_file=None, beta_file=None, constant_regions=False, append=True):  # main function for binding data
    create_gene_dic("./util/ref/family_seq.fasta")
    tcr_dic = get_tcr_info(data)
    # create full tcr sequences
    tcr_seq_dic = make_tcr_seq(tcr_dic)
    if constant_regions:
        tcr_seq_dic = append_constant(tcr_seq_dic, "TRAC*01", "TRBC1*01")
    #return tcr_seq_dic
    # join with input dataframe
    if append:
        data = return_sequences(tcr_seq_dic)
        return data
        #print(type(data))
        #print(data.dtypes)
    
    if alpha_file is not None and beta_file is not None:
        make_fasta_files(alpha_file, beta_file, tcr_seq_dic)
