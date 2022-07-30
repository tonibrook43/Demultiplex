# Define the problem
-The problem we are dealing with is determining where there is index hopping and unknown indexes in our R1 and R2 files. Index hopping typically occurs when you have a molecule cluster that is close enough to bridge over to the adjacent cluster and they swap indexes. In your sequenced read files, you might end up seeing a location where I1=A8 but I2=A5 when what you should really see is, I1 = A8 and I2=A8. We can solve this problem by demultiplexing which will separate out index swapped and unknown indexes from the known indexes. If this problem is not resolved, it could cause a lot of issues downstream.

# Determine/describe what output would be informative
-We will have the following output with a total of 52 FASTQ files:
1. 48 FASTQ files; 24 files will be for the 24 known index pairs for R1 and 24 files will be for the 24 known index pairs for R2.
	*These files will look something like A8.R1.fq.gz and A8.R2.fq.gz
2.  2 FASTQ files; 1 file for unknown indexes for R1 and 1 file for unknown indexes for R2
	*These files will look something like unknown.R1.fq.gz and unknown.R2.fq.gz
3. 2 FASTQ files; 1 file for known index swaps for R1 and 1 file for known index swaps for R2.
	*These files will look something like swapped.I1.fq.gz and swapped.I2.fq.gz
-Histograms of output information we might want to know:
1. Distribution of swapped indexes for R1 relative to swapped indexes for R2
2. Distribution of unknown indexes across R1 and R2

# Write examples (unit tests!)
Include four properly formatted input FASTQ files with read pairs that cover all three categories (dual matched, index-hopped, unknown index)
Include the appropriate number of properly formatted output FASTQ files given your input files
- 4 FASTQ input files: 
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz

-6 FASTQ output files examples:
A8.R1.fq.gz
A8.R2.fq.gz
unknown.R1.fq.gz
unknown.R2.fq.gz
swapped.I1.fq.gz
swapped.I2.fq.gz


# Develop your algorithm using pseudocode
	*import gzip, matplotlib as plt, bioinfo
1.	Open R1, R2, R3, R4 (provide absolute path to files) with gzip.open
2.	Create output file paths (needs variable) for known, unknown, and swapped
3.	Open Known Index file (provide absolute path from talapas)
    a.	Create set of known indexes
4.	‘for’ loop
    a.	Iterate over R1 FASTQ file
        i.	If index is known, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record
        ii.	If index unknown, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record
        iii.	If index is known by swapped, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record
    b.	Iterate over R2 FASTQ file
        i.	If index is known, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record
        ii.	If index unknown, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record
        iii.	If index is known by swapped, add index to header and output to appropriate file
            1.	Make sure to extract all four lines of record


# Determine High Level Functions
1.	Reverse Complement function:
    a.	Docstring: This function will first create the complement of my DNA sequence then create the reverse complement of my DNA sequence.
    b.	Pseudocode for reverse complement function:
        i.	def the function
        ii.	create a dictionary of keys and values (base and base complement)
        iii.	‘for’ loop
            1.	For each base in our sequence:
                a.	Iterate through dictionary keys and make complement
                    i.	Make variable for this output
                b.	Make an empty list for complement base
                c.	Append the complement base to empty list
            2.	return “empty” list (now has complement bases in it) 
            3.	add [::-1] to return statement to reverse the complement sequence
c.	Unit Test Example:
i.	Input: AATTCG
ii.	Output: CGAATT

