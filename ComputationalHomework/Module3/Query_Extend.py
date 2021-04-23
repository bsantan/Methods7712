######################################################################
##
##
## 
######################################################################
#
#The Query_Extend.py is an algorithm that takes a given set of Next Generation sequencing reads and aligns them to a given query sequence. This algorithm can optionally produce an intermediate file of all seuqence-query matches. The algorithm will also output a fasta file with the extended query, and a table of all sequencing reads which aligned to the query. A histogram of the length of reads which aligned will also be produced. It can be specified whether to extend the query based on the longest read possible or the longest read which aligned with other sequencing segments. 

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
import csv
import matplotlib.pyplot as plt


#Define arguments for each required and optional input
def defineArguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sequencer-file",  dest="SequencerFile",required=True,help="SequencerFile")
    parser.add_argument("--query-file",  dest="QueryFile",required=True,help="QueryFile")
    parser.add_argument("--min-match-length", dest="MinMatchLength",required=False,default=0,help="MinMatchLength",type=int)
    #Allow user to specify whether to output intermediate file of sequence mathces
    parser.add_argument("--intermediate-file-output", dest="IntermediateFileOutput",required=False,action="store_true",help="IntermediateFileOutput")
    parser.add_argument("--output-directory",  dest="OutputDirectory",required=True,help="OutputDirectory")
    parser.add_argument("--match-algorithm", dest="MatchAlgorithm",required=False,default="dense",help="MatchAlgorithm")

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
        read = ''

        if lines[i].startswith(">"):
            fasta_dict["Seq_ID"] = lines[i].replace(">","")
            i += 1
            while not lines[i].startswith(">") and i <= len(lines)-1:
                #Remove whitespace after each read to get correct length
                read += lines[i].strip()
                #Break when last index of file is reached, important for when file is only 2 lines long
                if (i == len(lines)-1): break
                i += 1
            fasta_dict["Read"] = read
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
    #for i in range(0,5000):
        #Check if part of sequencing read matches part of query read
        num_seq_comparisons = len(query_dict[0]["Read"])-len(seq_reads[i]["Read"])+1 + 2*(len(seq_reads[i]["Read"])-1)
        #Initialize counts which will define the window of comparison
        startcount = 1
        endcount = 1
        for j in range(min_match_length,num_seq_comparisons-min_match_length):
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

            #Evaluate if the next seqeunce segment also matches query, and skip the current one if so for beginning of sequence read
            if (j < len(seq_reads[i]["Read"])-1) and j > 0: 
                seq_match_next_start = seq_match_start - 1
                seq_match_next = seq_reads[i]["Read"][seq_match_next_start:seq_match_end]
                query_match_next = query_dict[0]["Read"][query_match_start:query_match_start+len(seq_match_next)]
                if (seq_match_next == query_match_next):
                    continue

            #Evaluate if the next seqeunce segment also matches query, and skip the current one if so for end of sequence read
            if (j > num_seq_comparisons - len(seq_reads[i]["Read"])+1) and (len(seq_reads[i]["Read"])-endcount+1) > 0: 
                seq_match_next_end = seq_match_start - 1
                seq_match_next = seq_reads[i]["Read"][seq_match_start:seq_match_next_end]
                query_match_next = query_dict[0]["Read"][query_match_start:query_match_start+len(seq_match_next)]
                if (seq_match_next == query_match_next):
                    continue

            #Check whether sequence and query segments match and are above defined minimum match length
            #For forward orientation
            if (seq_match == query_match and len(seq_match) > min_match_length and rev == False):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["qseqid"] = query_dict[0]["Seq_ID"].rstrip("\n")
                matching_kmer["sseqid"] = seq_reads[i]["Seq_ID"].rstrip("\n")
                matching_kmer["sstart"] = seq_match_start
                matching_kmer["send"] = seq_match_end
                matching_kmer["qstart"] = query_match_start
                matching_kmer["qend"] = query_match_start + len(seq_match)
                match.append(matching_kmer)
            #For reverse complement, reverse the coordinates of the sequence read that are recorded and record the sequence in reverse complement orientation
            if (seq_match == query_match and len(seq_match) > min_match_length and rev == True):
                matching_kmer["Sequence"] = seq_reads[i]["Read"]
                matching_kmer["qseqid"] = query_dict[0]["Seq_ID"].rstrip("\n")
                matching_kmer["sseqid"] = seq_reads[i]["Seq_ID"].rstrip("\n")
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

def verify_query_presence(input_dict,query_dict):

    index_match = []
    #Initialize boolean at true, assuming query will be found
    boolean = True

    #Find all indices of query that have match
    for i in range(0,len(input_dict)):
        for j in range(input_dict[i]["qstart"],input_dict[i]["qend"]):
            if (j not in index_match):
                index_match.append(j)

    for i in range(0,len(query_dict[0]["Read"])):
        if (i not in index_match):
            print("Part of query not found: index ",i,", ",query_dict[0]["Read"][i])
            #Set boolean to false if any part of query is not found
            boolean = False

    return boolean

def extend_query(input_dict,query_dict,output_dir):

    extension_segments_start = []
    extension_segments_end = []

    #Find all segments that overlap with beginning or end of query to extend
    for i in range(0,len(input_dict)):
        #When sequence read overlaps with beginning of query
        if (input_dict[i]["qstart"] == 0):
            #Check if read is reverse complemented and append beginning part of matching sequence in that orientation
            if (input_dict[i]["sstart"] > input_dict[i]["send"]):
                extension_segments_start.append(input_dict[i]["Sequence"][0:input_dict[i]["send"]])
            #When read is not reverse complement, append beginning part of mathcing sequence in normal orientation
            else:
                extension_segments_start.append(input_dict[i]["Sequence"][0:input_dict[i]["sstart"]])

        #When sequence read overlaps with end of query
        if (input_dict[i]["qend"] == len(query_dict[0]["Read"])):
            #Check if read is reverse complemented and append end part of matching sequence in that orientation
            if (input_dict[i]["sstart"] > input_dict[i]["send"]):
                extension_segments_end.append(input_dict[i]["Sequence"][input_dict[i]["sstart"]+1:len(input_dict[i]["Sequence"])])
            #When read is not reverse complement, append end part of matching sequence in normal orientation
            else:
                extension_segments_end.append(input_dict[i]["Sequence"][input_dict[i]["send"]+1:len(input_dict[i]["Sequence"])])

    #Generate a list of all extension to evaluate common lengths
    start_lengths = []
    end_lengths = []

    for i in range(len(extension_segments_start)):
        start_lengths.append(len(extension_segments_start[i]))

    for i in range(len(extension_segments_end)):
       end_lengths.append(len(extension_segments_end[i]))       

    #Output a histogram for lengths of all extensions at start and end of query if matches were found
    if(len(start_lengths) > 0):
        generate_histogram(start_lengths,"Start",output_dir)
    else: print("No matches for start of query found.")
    if(len(end_lengths) > 0):
        generate_histogram(end_lengths,"End",output_dir)
    else: print("No matches for end of query found.")

    return extension_segments_start,extension_segments_end

#Evaluate the extended query by identifying the longest segment at beginning and end of query to append
def extend_with_longest_segment(query_dict,extension_segments_start,extension_segments_end):

    #If no matches were found, append empty string to query
    if(len(extension_segments_start) == 0): extension_segments_start.append("")
    if(len(extension_segments_end) == 0): extension_segments_end.append("")

    #Find the longest extension in beginning and end
    start_extension = max(extension_segments_start,key=len)
    end_extension = max(extension_segments_end,key=len)

    final_query = start_extension+query_dict[0]["Read"]+end_extension

    if (query_dict[0]["Read"] in final_query): print("It Worked!")

    return final_query

#Evaluate the extended query by identifying the longest segment that is aligned with another segment
def extend_with_dense_match(query_dict,extension_segments_start,extension_segments_end):

    start_extension = find_longest_sequence(extension_segments_start,"start")
    end_extension = find_longest_sequence(extension_segments_end,"end")

    final_query = start_extension+query_dict[0]["Read"]+end_extension

    return final_query

#Find longest sequence that has at least one other segment aligned
def find_longest_sequence(extension_segments,segment):

    #If no matches to query were found, append empty string to query
    if(len(extension_segments) == 0): sequence = ""
    else:
        print("There is length for ext seg: ",len(extension_segments))
        #Find which segments match each other and store as dictionary of matching index with length
        matching_indices = []
        for i in range(len(extension_segments)):
            for j in range(len(extension_segments)):
                match_dict = {}
                min_length = min(len(extension_segments[i]),len(extension_segments[j]))
                if (segment == "end"):
                    first = extension_segments[i][0:min_length]
                    second = extension_segments[j][0:min_length]
                elif (segment == "start"):
                    first = extension_segments[i][len(extension_segments[i]) - min_length:len(extension_segments[i])]
                    second = extension_segments[j][len(extension_segments[j]) - min_length:len(extension_segments[j])]
                if first == second and i != j:
                    match_dict["Ref_Index"] = i
                    match_dict["Match_Index"] = j
                    match_dict["Length"] = min_length
                    matching_indices.append(match_dict)

        #If no matches between segments were found, append empty string to query
        if(len(matching_indices) == 0): sequence = ""
        else:
            print("There is length for indices: ",len(matching_indices))
            #Get max length of all matching segments
            max_match = max(matching_indices, key=lambda x:x['Length'])
            #Remove Length from dictionary to be able to compare only the 2 matching sequence lengths
            max_match.pop("Length", None)

            #Find the sequence of the shorter of the segments which aligned
            index_max_match = min(max_match, key=lambda k: max_match[k])
            sequence = extension_segments[max_match[index_max_match]]

    return sequence

#Generate a histogram of lengths of extensions
def generate_histogram(extension_list,name,output_dir):

    plt.hist(extension_list, bins = max(extension_list))
    plt.title("Length of "+name+" Extensions")
    plt.ylabel('# Occurences')
    plt.xlabel('Length of Extension')

    #Save histogram to user specified output folder.
    plt.savefig(output_dir+"/"+name+"Extensions.png")
    plt.clf()

def generate_alleles_file(final_query,output_dir):

    with open(output_dir+"/ALLELES.fasta", "w") as fasta_file:
        fasta_file.write(">EXTENDED_QUERY\n" + final_query)

def generate_table(input_dict,output_dir):

    subset_input_dict = []

    #Only include columns of interest
    columns = ['sseqid','qseqid','sstart','send','qstart','qend']
    for i in range(len(input_dict)):
        subset_input_dict.append(dict((k, input_dict[i][k]) for k in columns if k in input_dict[i]))

    with open(output_dir+"/ALLELES.aln", "w",newline='') as csvfile:
        dict_writer = csv.DictWriter(csvfile, subset_input_dict[0].keys(), delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(subset_input_dict)

def main():

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    SequencerFile = args.SequencerFile
    QueryFile = args.QueryFile
    MinMatchLength = args.MinMatchLength
    IntermediateFileOutput = args.IntermediateFileOutput
    OutputDirectory = args.OutputDirectory
    MatchAlgorithm = args.MatchAlgorithm

    #Command expected to execute the script
    cmd = "--sequencer-file %s --query-file %s --min-match-length %i --output-directory %s" % (SequencerFile,QueryFile,MinMatchLength,OutputDirectory)
    
    starttime = datetime.now()
    print("Minimum Match Length: ",MinMatchLength,"Intermediate File Output: ",IntermediateFileOutput,"Match Algorithm: ",MatchAlgorithm)

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

    print("Length of matching segments forward: ",len(seq_hash_forward))
    print("Length of matching segments reverse: ",len(seq_hash_reverse))

    #Validate that entirety of query match is found.
    query_presence = verify_query_presence(seq_hash_all,query)

    #Only perform extension if match to entire query is found, exit otherwise.
    if (query_presence is True):

        extend_start,extend_end = extend_query(seq_hash_all,query,OutputDirectory)

        #Use method according to user specification
        if MatchAlgorithm == "longest":
            query_extension = extend_with_longest_segment(query,extend_start,extend_end)
        
        if MatchAlgorithm == "dense":
            query_extension = extend_with_dense_match(query,extend_start,extend_end)

        print("Length of extended query: ",len(query_extension))

        generate_alleles_file(query_extension,OutputDirectory)

        generate_table(seq_hash_all,OutputDirectory)
    else:
        #Only printed if entire query match is not found.
        print("No output files will be generated as query match not found.")

    endtime = datetime.now()
    tdelta = endtime-starttime
    print("Total computation time: ",tdelta)

if __name__ == '__main__':
    main()
