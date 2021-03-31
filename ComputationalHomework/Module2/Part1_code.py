######################################################################
##
##
## 
######################################################################
#
#The Query_Extend_Part1.py is the first part of an algorithm that takes a given set of Next Generation sequencing reads and aligns them to a given query sequence. This algorithm is expected to be run in tandem with Part2 in the future, which will output a human readable form of the result. On its own, Part1 will produce an intermediate file of all seuqence-query matches.

##Libraries
import argparse
import sys
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
    #Allow user to specify whether to output intermediate file of sequence mathces
    parser.add_argument("--intermediate-file-output", dest="IntermediateFileOutput",required=False,action="store_true",help="IntermediateFileOutput")
    parser.add_argument("--output-directory",  dest="OutputDirectory",required='--intermediate-file-output' in sys.argv,help="OutputDirectory")

    return parser


#Create a dictionary for given fasta file with ID and Reads, assumes that format for ID includes '>' at beginning, and read exists for each ID
def parse_file(filename):

    fasta_list = []
    if filename.endswith(".gz"):
        f = gzip.open(filename,"rt")
    else:
        f = open(filename)
    lines = f.readlines()
    for i in range(0, len(lines)-1):
        fasta_dict = {}
        #Separate each read based on fasta file format of ">" character
        if lines[i].startswith(">"):
            fasta_dict["Seq_ID"] = lines[i].replace(">","")
            #Remove whitespace after each read to get correct length
            fasta_dict["Read"] = lines[i+1].strip()
            fasta_list.append(fasta_dict)

    return fasta_list

#Replace A <-> T and C <-> G for a given sequence, and reverse the order.
def reverse_comp(sequence):

    sequence = sequence.replace("A", "%temp%").replace("T", "A").replace("%temp%","T").replace("G", "%temp%").replace("C", "G").replace("%temp%","C")
    reverse = ''.join(reversed(sequence))

    return reverse 

#Reverse complement all sequences given a dictionary of reads and ID from a fasta file. 
def generate_reverse_sequence(reads):

    fasta_list_reverse = []

    #Loop through every read from a given dictionary
    for i in range(0,len(reads)):
        rev_dict = {}
        rev_dict["Read"] = reverse_comp(reads[i]["Read"])
        rev_dict["Seq_ID"] = reads[i]["Seq_ID"]
        #Append the reverse complemented read and original ID to a new dictionary
        fasta_list_reverse.append(rev_dict)

    return fasta_list_reverse


#Compare sequence kmers to query kmer using sliding window approach where the end of each sequencing read is compared to the beginning of the query and moved along the query for larger comparisons.
def sliding_window(query_dict,seq_reads,min_match_length,rev: bool = False):

    match = []

    for i in range(0,len(seq_reads)):
    #for i in range(0,50):
        #Check if part of sequencing read matches part of query read
        num_seq_comparisons = len(query_dict[0]["Read"])-len(seq_reads[i]["Read"])+1 + 2*(len(seq_reads[i]["Read"])-1)
        #Initialize counts which will define the window of comparison
        startcount = 1
        endcount = 1
        for j in range(0,num_seq_comparisons-min_match_length):
            matching_kmer = {}
            #Evaluation of match of end of sequencing read through first full comparison of sequencing read to beginning of query.
            if (j < len(seq_reads[i]["Read"])): 
                query_match_start = 0
                #Starting index of sequence read will decrease with each iteration
                seq_match_start = len(seq_reads[i]["Read"])-1-j
                #Ending index of sequence read will always be the final base pair
                seq_match_end = len(seq_reads[i]["Read"])-1
            #Evaluation of match of entire sequence read with a portion of query read.
            elif (len(seq_reads[i]["Read"]) <= j < (num_seq_comparisons - len(seq_reads[i]["Read"])+1)):
                #Starting index of query read will increase with each iteration
                query_match_start = startcount
                #Starting index of sequence read will always be the first base pair
                seq_match_start = 0
                #Ending index of sequence read will always be the final base pair
                seq_match_end = len(seq_reads[i]["Read"])-1
                startcount += 1
            #Evaluation of match of start of sequencing read through comparison of first base pair of sequence read to last base pair of query.
            elif (j >= num_seq_comparisons - len(seq_reads[i]["Read"])+1) and (len(seq_reads[i]["Read"])-endcount) > 0:
                #Starting index of query read will increase with each iteration
                query_match_start = startcount
                #Starting index of sequence read will always be the first base pair
                seq_match_start = 0
                #Ending index of sequence read will decrease with each iteration
                seq_match_end = len(seq_reads[i]["Read"])-endcount
                endcount += 1
                startcount += 1
            seq_match = seq_reads[i]["Read"][seq_match_start:seq_match_end]
            query_match = query_dict[0]["Read"][query_match_start:query_match_start+len(seq_match)]

            #Check whether sequence and query segments match and are above defined minimum match length
            #For forward orientation
            if (seq_match == query_match and len(seq_match) > min_match_length and rev == False):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["ID"] = seq_reads[i]["Seq_ID"]
                matching_kmer["sstart"] = seq_match_start
                matching_kmer["send"] = seq_match_end
                matching_kmer["qstart"] = query_match_start
                matching_kmer["qend"] = query_match_start + len(seq_match)
                match.append(matching_kmer)
            #For reverse complement, reverse the coordinates of the sequence read that are recorded
            if (seq_match == query_match and len(seq_match) > min_match_length and rev == True):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["ID"] = seq_reads[i]["Seq_ID"]
                matching_kmer["sstart"] = seq_match_end
                matching_kmer["send"] = seq_match_start
                matching_kmer["qstart"] = query_match_start
                matching_kmer["qend"] = query_match_start + len(seq_match)
                match.append(matching_kmer)

    return match

#Output json file of dictionary generated with all sequence-query matches and coordinates as intermediate.
def output_json(input_dict,output_dir):

    with open(output_dir+"/sequence_match.json", 'w') as json_file:
        json.dump(input_dict,json_file)


def main():

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    SequencerFile = args.SequencerFile
    QueryFile = args.QueryFile
    MinMatchLength = args.MinMatchLength
    IntermediateFileOutput = args.IntermediateFileOutput
    OutputDirectory = args.OutputDirectory

    #Command expected to execute the script
    cmd = "--sequencer-file %s --query-file %s --min-match-length %i --output-directory %s" % (SequencerFile,QueryFile,MinMatchLength,OutputDirectory)
    
    starttime = datetime.now()
    print("Minimum Match Length: ",MinMatchLength,"Intermediate File Output: ",IntermediateFileOutput)

    print("Start Time: ",starttime)

    ##Parse query and sequencing reads file and generate kmers
    seq_reads = parse_file(SequencerFile)
    query = parse_file(QueryFile)

    #Generate reverse complement of sequencing reads
    reverse_seq_reads = generate_reverse_sequence(seq_reads)

    #Find all matching sequencing and query reads and output into dictionary
    seq_hash_forward = sliding_window(query,seq_reads,MinMatchLength,rev = False)
    seq_hash_reverse = sliding_window(query,reverse_seq_reads,MinMatchLength,rev = True)

    #Combine forward and reverse matches
    seq_hash_all = seq_hash_forward + seq_hash_reverse

    #Output all matches to a json file if specified.
    if args.IntermediateFileOutput:
        output_json(seq_hash_all,OutputDirectory)


    endtime = datetime.now()
    tdelta = endtime-starttime
    print("Total computation time: ",tdelta)

    print("Length of matching segments forward: ",len(seq_hash_forward))
    print("Length of matching segments reverse: ",len(seq_hash_reverse))

if __name__ == '__main__':
    main()
