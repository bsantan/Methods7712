# Extending a Query Sequence

The Query_Extend.py is an algorithm that takes a given set of Next Generation sequencing reads and aligns them to a given query sequence. The algorithm will output a fasta file with the extended query, and a table of all sequencing reads which aligned to the query. It can be specified whether to extend the query based on the longest read possible or the longest read which aligned with other sequencing segments. The algorithm is only able to find exact matches to the given query, and does not evaluate N in the sequencing read as a potential match. The algorithm will also exit if no matches to the query are found, and the final fasta and table will not be output. 

## Getting Started

These instructions will provide the necessary environments, programs with installation instructions, and input files in order to run the Query_Extend.py script.

### Prerequisites
The following software packages must be installed:
```
- Python (any version later than 2.7)
```

## Running the Script

To run the Query_Extend.py script, execute the following command, where SEQUENCERFILE is the sequencing fasta file, QUERYFILE is the query fasta file as defined below, and OUTPUTDIRECTORY is the desired output directory:

```
Query_Extend.py --sequencer-file SEQUENCERFILE --query-file QUERYFILE --output-directory OUTPUTDIRECTORY
```

*Note that the specified output directory must already exist

### Optional Arguments

By default, the program will output a match between the sequence and query with a length as small as 1. However you can change the minimum match length by including the following optional argument:
```
--minimum-match-length MINMATCHLENGTH (default 0)
```

### Required Input Files

To run the Query_Extend script, a FASTA file from a next generation sequencing run and a FASTA file of the query of interest are required.

The format of a FASTA file is:
```
> Read Identifier
Read Sequence
```

For more information on FASTA file format, see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp. A fixed-width fasta file is also supported.


An example of the sequencer FASTA file with many reads is:
```
>2S43D:03629:08794
TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT
>2S43D:08938:01257
GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATG
```

The query FASTA file is expected to only contain one read. 

The program will not output an intermediate file by default. This file can still be generated when a match is not found for the entire query for troubleshooting purposes. To output an intermediate json file of all matches identified, include the following optional arguments:
```
--intermediate-file-output
```

The algorithm can extend the query by finding the longest possible segment (1) or by finding the longest segment which also matches other segments in the given matches (2), which is the default. The following commands can be used to specify the algorithm as (1) or (2) respectively:
```
--match-algorithm longest
--match-algorithm dense (default)
```

## Expected Outputs

The Query_Extend.py script will generate the following files (only if the entire query is identified):

### ALLELES.aln
A table containing the set of sequence reads which matched the query read in the following format:
sseqid  qseqid	sstart	send	qstart	qend
2S43D:08461:04180	contig1	13	40	1	64
2S43D:07701:07310	contig1	20	112	240	332
2S43D:07489:10315	contig1	123	90	20	53
2S43D:04035:14719	contig1	105	41	10	74

*Note that when the reverse complement of a seqencing read matches the query, the sstart will be greater than the send

### ALLELES.fasta
A fasta file containing the extended query in the following format (only if the entire query is identified):
```
>EXTENDED_QUERY
GAACAAGATGGATTGCAGGGATCGGCCATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGCCAGCTTGGGTGGAGA
```

### StartExtensions and EndExtensions Histograms
Histograms representing the length of partial matches to the beginning or end of the query sequence will be output as png files (only if the partial matches to beginning or end of query are identified, respectively).

The Query_Extend.py script will generate the following file when the --intermediate-file-output is enabled: 

### sequence_match.json
A file containing the set of sequence reads which matched the query read, including the following information about that match:

- Sequence: entire sequence of the sequencing read which matched
- ID: identifier of the sequence read which matched
- sstart: starting coordinate of the sequence read which matched
- send: ending coordinate of the sequence read which matched
- qstart: starting coordinate of the query read which matched
- qend: ending coordinate of the query read which matched

*Note that when the reverse complement of a seqencing read matches the query, the sstart will be greater than the send




