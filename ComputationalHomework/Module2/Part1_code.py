import gzip
import numpy as np
import copy
from datetime import datetime

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

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def reverse_comp(sequence):

    sequence = sequence.replace("A", "%temp%").replace("T", "A").replace("%temp%","T").replace("G", "%temp%").replace("C", "G").replace("%temp%","C")
    reverse = ''.join(reversed(sequence))

    return reverse 

#Break reads in given fasta file into kmer of length n.
def generate_kmer(reads,kmer_length):

    kmers = []

    for i in range(0,len(reads)):
        count = 0
        for j in chunker(reads[i]["Read"],kmer_length):
            kmer_dict = {}
            kmer_dict["Sequence"] = j
            kmer_dict["Start"] = count
            kmer_dict["End"] = count+len(j)
            kmer_dict["ID"] = reads[i]["Seq_ID"]
            count += kmer_length 
            kmers.append(kmer_dict)

    return kmers

#Break reads in given fasta file into kmer of length n.
def generate_reverse_kmer(reads,kmer_length):

    kmers = []

    for i in range(0,len(reads)):
        count = 0
        for j in chunker(reads[i]["Read"],kmer_length):
            kmer_dict = {}
            kmer_dict["Sequence"] = reverse_comp(j)
            kmer_dict["Start"] = count+len(j)
            kmer_dict["End"] = count
            kmer_dict["ID"] = reads[i]["Seq_ID"]
            count += kmer_length 
            kmers.append(kmer_dict)

    return kmers

#Create a matrix with the first row/column initialized to all gap_scores
def create_matrix(query_kmer_length,seq_kmer_length,gap_score):

    #Initialize array of kmer lengths + 1 to account for the first row/column of non-matches, matrix will be x=query by y=sequence
    matrix = np.empty((query_kmer_length+1,seq_kmer_length+1))
    matrix.fill(np.nan)
    matrix[[0], [0]] = 0
    matrix_reverse = np.empty((query_kmer_length+1,seq_kmer_length+1))
    matrix_reverse.fill(np.nan)
    matrix_reverse[[0], [0]] = 0

    for i in range(1,query_kmer_length+1):
        matrix[0,i] = matrix[0,i-1] + gap_score
        matrix[i,0] = matrix[i-1,0] + gap_score
        matrix_reverse[0,i] = matrix_reverse[0,i-1] + gap_score
        matrix_reverse[i,0] = matrix_reverse[i-1,0] + gap_score

    return matrix,matrix_reverse

#Loop through every kmer pair and evaluate for matches/mismatch to produce a hash of optimal alignments.
def eval_kmer_match(input_matrix,input_reverse_matrix,query_kmers,seq_kmers,mismatch_score,match_score,gap_score):

    seq_hash = []
  
    #for i in range(1,len(matrix.T)+1):
    for j in range(len(seq_kmers)):
        seq_char = len(seq_kmers[j]["Sequence"])
        for l in range(len(query_kmers)):

            seq_dict = {}

            query_char = len(query_kmers[l]["Sequence"])

            #Generate matrix for this kmer comparison
            matrix = copy.deepcopy(input_matrix)
            #print("PreDiag: ",j,l,matrix,seq_kmers[j]["Sequence"],query_kmers[l]["Sequence"])

            #Loop through both query and sequence characters to produce matrix of all maximum scores
            for k_query in range(query_char):
                for k_seq in range(seq_char):

                    matrix = compute_matrix(k_query,k_seq,matrix,seq_kmers[j]["Sequence"][k_seq],query_kmers[l]["Sequence"][k_query],mismatch_score,match_score,gap_score)
    
            #print(matrix,"Q: ",query_kmers[l]["Sequence"],"S: ",seq_kmers[j]["Sequence"])

            #Traverse resulting matrix for optimal alignment
            query_seq,read_seq = traverse_matrix(seq_kmers[j]["Sequence"],query_kmers[l]["Sequence"],matrix)

            #Record each value in the dictionary for this query/kmer pair
            seq_dict["rindex"] = j
            seq_dict["sequence_r"] = read_seq
            seq_dict["sstart"] = seq_kmers[j]["Start"]
            seq_dict["send"] = seq_kmers[j]["End"]
            seq_dict["sseqid"] = seq_kmers[j]["ID"]
            seq_dict["qindex"] = l
            seq_dict["sequence_q"] = query_seq
            seq_dict["qstart"] = query_kmers[j]["Start"]
            seq_dict["qend"] = query_kmers[j]["End"]
            seq_dict["queryid"] = query_kmers[j]["ID"]


            seq_hash.append(seq_dict)

    return seq_hash

#Based on paired values from given kmer pair, calculate maximum score by moving in x,y, or diagonal direction
def compute_matrix(base_k,k,matrix,seq_kmers_val,query_kmers_val,mismatch_score,match_score,gap_score):

    #print("base_k3: ",base_k,"k: ",k)
    #print(seq_kmers_val,query_kmers_val)
    #Check for maximum diagonal move
    if (seq_kmers_val == query_kmers_val):
        matrix_diag = matrix[k,base_k] + match_score
        #print("Match: ",seq_kmers_val,query_kmers_val)
    elif (seq_kmers_val != query_kmers_val and seq_kmers_val != "N"):
        matrix_diag = matrix[k,base_k] + mismatch_score
        #print("Mismatch: ",seq_kmers_val,query_kmers_val)
    #If seq val is N, add 0
    else:
        matrix_diag = matrix[k,base_k] + 0

    #Check for gap move X
    matrix_gap_x = matrix[k+1,base_k] + gap_score
    matrix_gap_y = matrix[k,base_k+1] + gap_score

    #print("Comp_values: ",matrix_diag,matrix_gap_x,matrix_gap_y)
    matrix[k+1,base_k+1] = max(matrix_diag,matrix_gap_x,matrix_gap_y)

    return matrix


#Evaluate best path within matrix and generate a sequence comparison between sequence and query kmer. Query sequence should always be on x, sequence on y.
def traverse_matrix(seq_kmer,query_kmer,matrix):

    #Represents the outcome of kmer comparison in reverse as <query base pair on x, seq base pair on y>
    rev_sequence = []

    #Check if kmer comparison was not the same length first, and add gaps if not

    #Sequence kmer is greater than query kmer so add gap in query matched to character in sequence
    if (len(seq_kmer) - len(query_kmer) > 0):
        diff_x = len(seq_kmer) - len(query_kmer)
        diff_y = 0
        for i in range(0,diff_x):
            rev_sequence.append(["-",seq_kmer[len(seq_kmer)-1-i]])
    #Query kmer is greater than sequence kmer so add gap in sequence matched to character in query
    elif (len(seq_kmer) - len(query_kmer) < 0):
        diff_y = len(query_kmer) - len(seq_kmer)
        diff_x = 0
        for i in range(0,diff_y):
            rev_sequence.append([query_kmer[len(query_kmer)-1-i],"-"])
    #Query and sequence kmers are same length, so pair up the last 2 characters of each kmer
    elif (len(seq_kmer) == len(query_kmer)):
        diff_x = 0
        diff_y = 0
        rev_sequence.append([query_kmer[len(query_kmer)-1],seq_kmer[len(seq_kmer)-1]])

    #print("Rev_seq start: ",rev_sequence)

    #Evaluate max value from diagonal, gap_x, or gap_y move
    #Represent seq
    x_coord = matrix.shape[0]-1-diff_x
    #Represents query
    y_coord = matrix.shape[1]-1-diff_y

    #print("Start matrix val: ",matrix[x_coord,y_coord],x_coord,y_coord,len(seq_kmer),len(query_kmer))
    #print("Q: ",query_kmer[x_coord-1],"S: ",seq_kmer[y_coord-1])

    #Evaluate whether moving back in x by 1, back in y by 1, or back diagonally by 1 optimizes score. Stop when coordinates are > 1 such that final calculation will result in first base pair of each sequence.
    while x_coord > 1 and y_coord > 1:
        y_move = matrix[x_coord,y_coord-1]
        x_move = matrix[x_coord-1,y_coord]
        diag_move = matrix[x_coord-1,y_coord-1]
        moves = (y_move,x_move,diag_move)
        #print("Moves: ",moves)
        best_move = moves.index(max(moves))

        #Query and sequence kmer are always 1 less length than coordinate in matrix
        if (best_move == 0): 
            y_coord -= 1
            rev_sequence.append(["-",seq_kmer[y_coord-1]]) 
        elif (best_move == 1): 
            x_coord -= 1
            rev_sequence.append([query_kmer[x_coord-1],"-"]) 
        elif (best_move == 2): 
            x_coord -= 1
            y_coord -= 1
            rev_sequence.append([query_kmer[x_coord-1],seq_kmer[y_coord-1]])

    #Evaluate sequence in forward orientation
    forw_sequence = list(reversed(rev_sequence))
    query_sequence = []
    seq_sequence = []

    for i in range(len(forw_sequence)):
        query_sequence.append(forw_sequence[i][0])
        seq_sequence.append(forw_sequence[i][1])

    #Concatenate each character to generate a string
    query_sequence = ''.join(query_sequence)
    seq_sequence = ''.join(seq_sequence)

    #print("FinalX: ",x_coord,"FinalY: ",y_coord)

    print("Q: ",query_sequence)
    print("S: ",seq_sequence)

    return query_sequence,seq_sequence  

def main():

    ##Parse query and sequencing reads file and generate kmers
    seq_reads = parse_file('/Users/brooksantangelo/Documents/Methods7712/Module2/SANTANGELO_Module2Day3LikeHW_Part1_Code/READS.fasta.gz')

    query = parse_file('/Users/brooksantangelo/Documents/Methods7712/Module2/SANTANGELO_Module2Day3LikeHW_Part1_Code/QUERY.fasta')

    start = datetime.now()
    print("Start Time: ",start)

    query_kmer = generate_kmer(query,10)
    reads_kmer = generate_kmer(seq_reads,10)
    reads_kmer_reverse = generate_reverse_kmer(seq_reads,10)

    print(len(query_kmer))
    print(len(reads_kmer))

    ###This is for testing
    query_kmer_test = query_kmer[1:5]
    reads_kmer_test = reads_kmer[1:5]
    #####

    ##Create a scoring matrix and evaluate kmer sequences

    matrix,reverse_matrix = create_matrix(10,10,-2)

    seq_hash = eval_kmer_match(matrix,reverse_matrix,query_kmer_test,reads_kmer_test,-1,1,-2)

    print(len(seq_hash))
    print(seq_hash[0])
    print(seq_hash[len(seq_hash)-1])

    end = datetime.now()
    print("End Time: ",end)

    tdelta = end-start
    print(tdelta)

if __name__ == '__main__':
    main()

