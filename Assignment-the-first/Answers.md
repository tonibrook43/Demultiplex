# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 bp | Phred+ 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 bp | Phred+ 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 bp | Phred+ 33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 | 101 bp | Phred+ 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE** Histograms in repository
    3. **YOUR ANSWER HERE**

number of reads in each file:

    1294_S1_L008_R1_001.fastq.gz
    number of reads: 1452986940/4= 363246735

    1294_S1_L008_R2_001.fastq.gz
    number of reads: 1452986940/4= 363,246,735

    1294_S1_L008_R3_001.fastq.gz
    number of reads: 1452986940/4= 363,246,735

    1294_S1_L008_R4_001.fastq.gz
    number of reads: 1452986940/4= 363,246,735

Part One Questions and Answers:

Q: What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.

A: A good quality score cutoff for biological read pairs and our index reads would be between 30 since that is the minimum mean quality score for both R1 and R2 and our index reads. This can be seen in all four of the histograms, the smallest quality score seen is 30. If we did a cutoff at 30, we would still capture all of our read pairs with a 99% accuracy.

Q: How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)

For R2:
command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | grep -v "@\|#\|+" | grep -o "N" | wc -l | head
output: 3976613 (this is about 1.0% of our total base calls)

For R3:
command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | grep -v "@\|#\|+" | grep -o "N" | wc -l | head
command: 3329901 (this is about 0.9% of our total base calls)

Total undertermined base calls: 7306514
    
## Part 2 (found in pseudocode.md)
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
    a. Known Indexes for test files:
        TT
        CC
        AC
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
