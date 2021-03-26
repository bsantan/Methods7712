##Libraries
import argparse
import gzip
import numpy as np
import copy
from datetime import datetime
from difflib import SequenceMatcher
from collections import defaultdict
import json


#Define arguments for each required and optional input
def defineArguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sequencer-file",  dest="SequencerFile",required=True,help="SequencerFile")
    parser.add_argument("--query-file",  dest="QueryFile",required=True,help="QueryFile")
    parser.add_argument("--min-match-length", dest="MinMatchLength",required=False,default=0,help="MinMatchLength",type=int)

    return parser


#Create a dictionary for given fasta file with ID and Reads, assumes that format for ID includes > at beginning, and read exists for each ID
def parse_file(filename):

    fasta_list = []
    if filename.endswith(".gz"):
        f = gzip.open(filename,"rt")
    else:
        f = open(filename)
    lines = f.readlines()
    for i in range(0, len(lines)-1):
        fasta_dict = {}
        if lines[i].startswith(">"):
            fasta_dict["Seq_ID"] = lines[i].replace(">","")
            fasta_dict["Read"] = lines[i+1]
            fasta_list.append(fasta_dict)

    return fasta_list

def reverse_comp(sequence):

    sequence = sequence.replace("A", "%temp%").replace("T", "A").replace("%temp%","T").replace("G", "%temp%").replace("C", "G").replace("%temp%","C")
    reverse = ''.join(reversed(sequence))

    return reverse 

#Reverse complement given sequence.
def generate_reverse_sequence(reads):

    fasta_list_reverse = []

    for i in range(0,len(reads)):
        rev_dict = {}
        rev_dict["Read"] = reverse_comp(reads[i]["Read"])
        rev_dict["Seq_ID"] = reads[i]["Seq_ID"]
        fasta_list_reverse.append(rev_dict)

    print(reads[0])
    print(fasta_list_reverse[0])

    return fasta_list_reverse


#Compare sequence kmers to query kmer using sliding window approach
def sliding_window(query_dict,seq_reads,min_match_length,rev: bool = False):

    match = []

    for i in range(0,len(seq_reads)):
    #for i in range(0,500):
        #Check if part of sequencing read matches part of query read
        num_seq_comparisons = len(query_dict[0]["Read"])-len(seq_reads[i]["Read"])+1 + 2*(len(seq_reads[i]["Read"])-1)
        startcount = 1
        endcount = 1
        for j in range(1,num_seq_comparisons-min_match_length):
            matching_kmer = {}
            if (j < len(seq_reads[i]["Read"])): 
                query_match_start = 0
                seq_match_start = len(seq_reads[i]["Read"])-1-j
                seq_match_end = len(seq_reads[i]["Read"])-1
            elif (len(seq_reads[i]["Read"]) <= j < (num_seq_comparisons - len(seq_reads[i]["Read"])-1)):
                query_match_start = startcount
                seq_match_start = 0
                seq_match_end = len(seq_reads[i]["Read"])-1
                startcount += 1
            elif (j >= num_seq_comparisons - len(seq_reads[i]["Read"])-1) and (len(seq_reads[i]["Read"])-endcount) > 0:
                query_match_start = startcount
                seq_match_start = 0
                seq_match_end = len(seq_reads[i]["Read"])-endcount
                endcount += 1
                startcount += 1
            seq_match = seq_reads[i]["Read"][seq_match_start:seq_match_end]
            query_match = query_dict[0]["Read"][query_match_start:query_match_start+len(seq_match)]

            if (seq_match == query_match and len(seq_match) > min_match_length and rev == False):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["ID"] = seq_reads[i]["Seq_ID"]
                matching_kmer["sstart"] = seq_match_start
                matching_kmer["send"] = seq_match_end
                matching_kmer["qstart"] = query_match_start
                matching_kmer["qend"] = query_match_start + len(seq_match)
                match.append(matching_kmer)
            if (seq_match == query_match and len(seq_match) > min_match_length and rev == True):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["ID"] = seq_reads[i]["Seq_ID"]
                matching_kmer["sstart"] = seq_match_end
                matching_kmer["send"] = seq_match_start
                matching_kmer["qstart"] = query_match_start
                matching_kmer["qend"] = query_match_start + len(seq_match)
                match.append(matching_kmer)

    return match

def combine_matches(match_dict):

    new_dict = []

    for i in range(0,len(match_dict)):
        for j in range(1,len(match_dict)):
            new_match = defaultdict(list)
            #Check for extension in reverse direction
            if (match_dict[j]["qstart"] - match_dict[i]["qend"] == 1):
                new_match["Sequence"].append(match_dict[i]["Sequence"])
                new_match["Sequence"].append(str(match_dict[j]["Sequence"]))
                new_match["ID"].append(match_dict[i]["ID"])
                new_match["ID"].append(match_dict[j]["ID"])
                new_match["sstart"].append(match_dict[i]["sstart"])
                new_match["send"].append(match_dict[i]["send"])
                new_match["sstart"].append(match_dict[j]["sstart"])
                new_match["send"].append(match_dict[j]["send"])
                new_match["qstart"].append(match_dict[i]["qstart"])
                new_match["qend"].append(match_dict[j]["qend"])
                new_dict.append(new_match)

    return new_dict

def output_json(input_dict,output_filename):

    with open(output_filename, 'w') as json_file:
        json.dump(input_dict,json_file)


def main():

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    SequencerFile = args.SequencerFile
    QueryFile = args.QueryFile
    MinMatchLength = args.MinMatchLength

    #Command expected to execute the script
    cmd = "--sequencer-file %s --query-file %s --min-match-length %i" % (SequencerFile,QueryFile,MinMatchLength)

    ##Parse query and sequencing reads file and generate kmers
    seq_reads = parse_file(SequencerFile)
    query = parse_file(QueryFile)

    #Generate reverse complement of sequencing reads
    reverse_seq_reads = generate_reverse_sequence(seq_reads)

    #Find all matching sequencing and query reads and output into dictionary
    seq_hash_forward = sliding_window(query,seq_reads,MinMatchLength,rev = False)
    #print("F: ",seq_hash_forward)
    seq_hash_reverse = sliding_window(query,reverse_seq_reads,MinMatchLength,rev = True)
    #print("R: ",seq_hash_reverse)

    #Combine forward and reverse matches
    seq_hash_all = seq_hash_forward + seq_hash_reverse

    #Output all matches to a json file
    output_json(seq_hash_all,"/Users/brooksantangelo/Documents/Methods7712/Module2/SANTANGELO_Module2Day3LikeHW_Part1_Code/sequence_match.json")

    forward_matches = combine_matches(seq_hash_forward)
    #print("Combined F: ",forward_matches)
    reverse_matches = combine_matches(seq_hash_reverse)
    #print("Combined R: ",reverse_matches)

    print("Length of matching kmers forward: ",len(seq_hash_forward),"Combined: ",len(forward_matches))
    print("Length of matching kmers reverse: ",len(seq_hash_reverse),"Combined: ",len(reverse_matches))

if __name__ == '__main__':
    main()
