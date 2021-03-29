# Extending a Query Sequence, Part 1

The Query_Extend_Part1.py is the first part of an algorithm that takes a given set of Next Generation sequencing reads and aligns them to a given query sequence. The outcome of the algorithm in its entirety is the longest sequence contig that can be constructed from the sequencing reads along with metadata indicating the coordinates of these alignments. Both forward and reverse complement sequencing reads will be evaluated in alignment to the query read. Part1 will produce an intermediate output indicating all sequence matches identified between the sequencing reads and the query read as a json file. Part2 will extend this result and produce the contig and metadata in a human readable format. 

## Getting Started

These instructions will provide the necessary environments, programs with installation instructions, and input files in order to run the Query_Extend_Part1.py script.

### Prerequisites
The following software packages must be installed:
```
- Python (any version later than 2.7)
```

### Required Input Files

To run the Query_Extend_Part1 script, a FASTA file from a next generation sequencing run and a FASTA file of the query of interest are required.

The format of a FASTA file is:
```
> Read Identifier
Read Sequence
```

For more information on FASTA file format, see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp

An example of the sequencer FASTA file with many reads is:
```
>2S43D:03629:08794
TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT
>2S43D:08938:01257
GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATG
```

The query FASTA file is expected to only contain one read. An example of the query FASTA file is:
```
>INITIAL_QUERY
GGGATCGGCCATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACT
```
## Running the Script

To run the Query_Extend_Part1.py script, execute the following command, where SEQUENCERFILE is the sequencing fasta file and QUERYFILE is the query fasta file as defined above:

```
Query_Extend_Part1.py --sequencer-file SEQUENCERFILE --query-file QUERYFILE
```

### Optional Arguments

By default, the program will output a match between the sequence and query with a length as small as 1. However you can change the minimum match length by including the following optional argument:
```
--minimum-match-length MINMATCHLENGTH (default 0)
```

The program will not output an intermediate file by default, however you can specify that an intermediate json file of all matches identified be output by including the following optional arguments:
```
--intermediate-file-output
--output-directory OUTPUTDIRECTORY
```

*Note that when --intermediate-file-output is provided, --output-directory is required

## Expected Outputs

The Query_Extend_Part1.py script will generate the following file when the --intermediate-file-output is enabled: 

### sequence_match.json
A file containing the set of sequence reads which matched the query read, including the following information about that match:

- Sequence: entire sequence of the sequencing read which matched
- ID: identifier of the sequence read which matched
- sstart: starting coordinate of the sequence read which matched
- send: ending coordinate of the sequence read which matched
- qstart: starting coordinate of the query read which matched
- qend: ending coordinate of the query read which matched

*Note that when the reverse complement of a seqencing read matches the query, the sstart will be greater than the send



