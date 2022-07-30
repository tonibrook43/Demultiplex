#!/usr/bin/env python3.10
import argparse
import numpy as np
import matplotlib.pyplot as plt
import bioinfo
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to read and output our input files for demultiplexing")
    parser.add_argument("-r", "--read-length", help="read length", type=int)
    parser.add_argument("-n", "--number-reads", help="number of read in file", type=int)
    parser.add_argument("-f", "--file", help="Name of file", required=True)
    parser.add_argument("-o", "--output", help="Name of output file", required=True)
    return parser.parse_args()
args=get_args()

#Read files I will be parsing through: we have to take these files from the shared directory
R1= "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2= "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3= "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4= "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"


with gzip.open(args.file, 'rt') as myfile: #add gzip.open whne running full file. use rt so that gzip.open will read it as text and not binary.
    list1 = np.zeros([args.read_length]) #np is actually creating a placeholder in the array with floats and is not actually empty. replaced with zeros then we know its a zero. The error I was receiving is that my list of sums were coming up as floats and completely wrong. once replaced with np.zero, the list came out correctly!
    #mean=np.empty([args.read_length])
    count_line = 0 
    for line in myfile:
    #myline = myfile.readline()
        #count_line += 1
        #line=line.decode("ascii") # when you unzip the file, it will give an encoder error because of how its intrepreting characters. This will turnr it into something that it can intrepret. 
        line = line.strip()
        count_line += 1
        if count_line %4 == 0:
            character_count=0
            for character in line:
                qs = bioinfo.convert_phred(character) #we have to assign this command a variable in order for python to remember it
                list1[character_count] += qs
                character_count +=1
                #print(character_count)
#print(list1) #used to check that my list is populating 
#print(len(list1)) #used to check that my list is the length it should be

mean=[]
for sum in list1:
   mean.append(sum/(count_line/4))
print(mean)

#array
x= range(args.read_length) #this will make x a list of the range 101. (x is technically a range object)
y= mean
# print(x)
# print(y)
# Function to plot
plt.bar(x, y)

#Axis Labels and Title name
plt.ylabel("Mean Quality score", fontsize=10)
plt.xlabel("Position in Read", fontsize=10)

plt.title(f"Distribution of Mean Quality Scores vs Position", fontsize=10)


plt.savefig(args.output+".png") #this will pull the first part of the figure name from my command line (-o)

   
# mean=np.mean(list1)
# print(mean)      

                    
