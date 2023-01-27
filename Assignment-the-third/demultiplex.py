#!/usr/bin/env python3.10
import argparse
import numpy
import gzip
import itertools

#use itertools for permutations

def get_args():
    parser = argparse.ArgumentParser(description="A program to read and output our input files for demultiplexing")
    #parser.add_argument("-qs", "--quality-score", help="quality score cutoff value", type=int)
    parser.add_argument("-f1", "--file1", help="output file for read 1")
    parser.add_argument("-f2", "--file2", help="output file for read 2")
    parser.add_argument("-f3", "--file3", help="output file for read 3")
    parser.add_argument("-f4", "--file4", help="output file for read 4")
    return parser.parse_args()
args=get_args()



#Read files I will be parsing through: we have to take these files from the shared directory
file1= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz", "rt")
file2= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz", "rt")
file3= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz", "rt")
file4= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz", "rt")



index_txt= open("/projects/bgmp/shared/2017_sequencing/indexes.txt", "r") 
index_set= set() 
rec_dict={} 
base_dict={"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
count_matched_dict={"matched": 0}
count_unknown_dict_r2={"unknown_r2": 0} #creating a dictionary for unknown indexes in r2 -- used to make a count for unknown + low qs
count_unknown_dict_r3={"unknown_r3": 0} #creating a dictionary for unknown indexes in r3 -- used to make a count for unknown + low qs
count_hopped_dict={"hopped": 0}


# **FUNCTIONS NEEDED**
def rev_comp(DNA): #we are looping over the dna string so we are looking at each position in our string and saving it to base
    cDNA=""
    for base in DNA:
        cDNA += base_dict[base] #need to make a list to store my complemented bases
    return cDNA[::-1]

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    convert_phred_int = (ord(letter) - 33)
    return convert_phred_int


#this loop will go through the index file, strip the new character line, split the columns, take column 4, and return my index set. 
for x in index_txt:
    x= x.strip().split()[4]
    if x == 'index':
        continue #if x is equal to index, go on to the next item/line but do not print
    index_set.add(x)
#print(index_set)

## itertools.product -- needed so that we capture both AB and BA as 2 occurances and not 1. In this case, we want to take matches from index_set##
match_products={}
p=itertools.product(index_set, index_set)  
for match in p:
    match_products[match]= 0
    #print(match)

# Creating my index files in a for loop so I can capture each index and add it to my file name -- used later to store my condition output to my files
indexes={}
for index in index_set:
    known_r1=open(index + '_R1.fq', 'wt' )
    known_r2=open(index + '_R2.fq', 'wt')
    indexes[index]=[known_r1, known_r2]

# Creating my unknown index file variables -- used later to store my condition output to my files
unknown_r1=open('unknown_r1.fq', 'wt')
unknown_r2=open('unknown_r2.fq', 'wt')
hopped_r1=open('hopped_r1.fq', 'wt')
hopped_r2=open('hopped_r2.fq', 'wt')

# Opening and reading my R1, R2, R3, and R4 files
while True:
    header_r1=file1.readline().strip()
    header_r2=file2.readline().strip()
    header_r3=file3.readline().strip()
    header_r4=file4.readline().strip()
    if header_r1 == "":
        break
    seq_r1=file1.readline().strip()
    seq_r2=file2.readline().strip()
    seq_r3=file3.readline().strip()
    seq_r4=file4.readline().strip()

    plus_r1=file1.readline().strip()
    plus_r2=file2.readline().strip()
    plus_r3=file3.readline().strip()
    plus_r4=file4.readline().strip()

    qs_r1=file1.readline().strip()
    qs_r2=file2.readline().strip()
    qs_r3=file3.readline().strip()
    qs_r4=file4.readline().strip()

    index_header_r1= header_r1+ ' ' + seq_r2 + '+' + rev_comp(seq_r3) #use both indexes seq files to check for hopped indexes
    index_header_r4= header_r4 + ' ' + seq_r2 + '+' + rev_comp(seq_r3)
   # print(index_header_r1, index_header_r4)

    written= False
    qs=30
    if seq_r2 not in index_set:
        unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
        unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
        #print("unknown_index2:", seq_r2)
        count_unknown_dict_r2["unknown_r2"] += 1
    elif rev_comp(seq_r3) not in index_set:
        unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
        unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
        #print("unknown_index3:", rev_comp(seq_r3))
        count_unknown_dict_r3["unknown_r3"] += 1
    elif seq_r2 == rev_comp(seq_r3):
        for phred_score in qs_r2:
            if convert_phred(phred_score) < qs:
                unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                count_unknown_dict_r2["unknown_r2"] += 1
                #print("matched with <30 qs:", seq_r2, rev_comp(seq_r3))
                written=True
                break
        if written == False:
            for phred_score in qs_r3:
                if convert_phred(phred_score) < qs:
                    unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                    unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                    count_unknown_dict_r3["unknown_r3"] += 1
                    written=True
                    break
            if written == False:
                indexes[seq_r2][0].write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                indexes[seq_r2][1].write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                #print("matched with >= 30 qs:", seq_r2, rev_comp(seq_r3))
                count_matched_dict["matched"] += 1
                match_products[(seq_r2, rev_comp(seq_r3))] += 1  #use another set of parenthesis to make a tuple
    else:
        for phred_score in qs_r2:
            if convert_phred(phred_score) < qs:
                unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                count_unknown_dict_r2["unknown_r2"] += 1
                written=True
                break
        if written == False:
            for phred_score in qs_r3:
                if convert_phred(phred_score) < qs:
                    unknown_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                    unknown_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                    count_unknown_dict_r2["unknown_r2"] += 1
                    written=True
                    break
            if written == False:
                hopped_r1.write(index_header_r1 + "\n" + seq_r1 + "\n" + plus_r1 + "\n" + qs_r1 + "\n")
                hopped_r2.write(index_header_r4 + "\n" + seq_r4 + "\n" + plus_r4 + "\n" + qs_r4 + "\n")
                #print("Hopped match:", seq_r2, rev_comp(seq_r3))
                count_hopped_dict["hopped"] += 1
                match_products[(seq_r2, rev_comp(seq_r3))] += 1
                written=True
#print(count_unknown_dict_r2, count_unknown_dict_r3, count_matched_dict, count_hopped_dict)


#To close our files
unknown_r1.close()
unknown_r2.close()
hopped_r1.close()
hopped_r2.close()
known_r1.close()
known_r2.close()


## Statistics ##

hopped_sum= sum(count_hopped_dict.values())
print("total hopped:", hopped_sum)

matched_sum= sum(count_matched_dict.values())
print("total matched:", matched_sum)

unknown_r1_sum= sum(count_unknown_dict_r2.values())
print("total read 1 unknown:", unknown_r1_sum)

unknown_r2_sum= sum(count_unknown_dict_r3.values())
print("total read 2 unknown:", unknown_r2_sum)

Total_Indexes= hopped_sum + matched_sum + unknown_r1_sum + unknown_r2_sum
print(Total_Indexes)

percent_matched= (matched_sum/Total_Indexes)*100
print(percent_matched)

percent_hopped= (hopped_sum/Total_Indexes)*100
print(percent_hopped)

percent_unknown= ((unknown_r1_sum + unknown_r2_sum)/Total_Indexes)*100
print(percent_unknown)

#this was from itertools and made a loop to access my pairs of matched products
for pairs in match_products:
    print(match_products[pairs], pairs)

with open("Index_Statistics.md", "wt") as statistics:
    statistics.write("Hopped Indexes= " + str(hopped_sum)+ "\n")
    statistics.write("Matched Indexes= " + str(matched_sum)+ "\n")
    statistics.write("Read 1 Unknown Indexes= " + str(unknown_r1_sum)+ "\n")
    statistics.write("Read 2 Unknown Indexes= " + str(unknown_r2_sum)+ "\n")
    statistics.write("Total Indexes in Read Files= " + str(Total_Indexes)+ "\n")
    statistics.write("Percent of Matched Reads= " + str(percent_matched) + "\n")
    statistics.write("Percent of Hopped Reads= " + str(percent_hopped) + "\n")
    statistics.write("Percent of Unknown Reads= " + str(percent_unknown) + "\n")


#For further output information:
    #Create a dictionary with actual index names so we can write the hopped and matched to our statistics file
    #file names and where they went (hopped, matched, unknown)
    #graph for phred scores?
    #graph for number of hopped, matched, and unknown indexes for each file

                    
