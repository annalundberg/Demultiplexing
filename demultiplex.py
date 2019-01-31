#!/usr/bin/env python

#SBATCH --partition=short
#SBATCH --job-name=demult
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=annalundberg92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

'''Bi624: Demultiplexing assignment
This is an algorithm for demultiplexing sequence reads/indexes.
We have 4 fastq files, 2 indexes and 2 sequences. We need to determine if the
indexes match per read. One sequence file is forward reads, the other is reverse.
The reads appear in a consistent order between files. We need to filter the
reads so that each index gets a file of forward reads and one of reverse reads.
Since we don’t know which file has forward and which has reverse, we will
arbitrarily assign. Reads with unclear indexes will be thrown into the pair of
files with unknown indexes.'''
'''Output: File names of all files created (2 files per index (forward and
reverse) and 2 files of reads with unclear indexes), integer of properly matched
indexes and float percentage of index hopping observed.'''

import argparse
import gzip

def get_arguments():
    parser = argparse.ArgumentParser(description="input paired-end reads and index files for demultiplexing")
    parser.add_argument("-r1", "--read1", help="name of file for read 1",\
                        required=True, type=str)
    parser.add_argument("-r2", "--read2", help="name of file for read 2",\
                        required=True, type=str)
    parser.add_argument("-i1", "--index1", help="name of file for index 1",\
                        required=True, type=str)
    parser.add_argument("-i2", "--index2", help="name of file for index 2",\
                        required=True, type=str)
    return parser.parse_args()

def complement_seq(seq):
    '''(string)-> string
    This fxn returns the reverse complement of a sequence'''
    comp_base = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    newseq=''.join([comp_base[base] for base in seq[::-1]])
    return newseq

def eval_quality(qual_scores):
    '''(string) -> float
    This fxn takes in a string of phred quality scores. It uses convert_phred()
    to convert ASCII characters into numberic quality. It takes the average of
    the quality scores and returns that float'''
    qual=0
    nums=0
    for ch in qual_scores:
        qual+=(ord(ch)-33)
        nums+=1
    avg_qual = qual/nums
    return avg_qual

def fix_n(in1,in2):
    '''(string,string)-> string, string
    This fxn takes in 2 indexes and uses each to correct 1st base N, or last
    in index 2 (b/c reverse complement)'''
    if in1[0] == 'N':
        in1=''.join((in2[0],in1[1:len(in1)]))
    if in2[len(in2)-1] == 'N':
        in2=''.join((in2[0:(len(in1)-1)],in1[len(in1)-1]))
    return in1, in2

def compare_fastq(r1, i1, i2, r2):
    '''(str, str, str, str)->dict, dict
    This fxn opens fwd and rev index and sequence files, compares the indexes,
    and uses write_files() to write the fwd and rev files into
    new files by indexes'''
    lnct = match_ct = hop_ct = 0 #init counts
    ns='Nn' #establish what an n is for index checking
    indexes = {'GTAGCGTA':[open('b1_fwd.fastq','w+'),open('b1_rev.fastq','w+')], 'CGATCGAT':[open('a5_fwd.fastq','w+'),open('a5_rev.fastq','w+')],
            'GATCAAGG':[open('c1_fwd.fastq','w+'),open('c1_rev.fastq','w+')], 'AACAGCGA':[open('b9_fwd.fastq','w+'),open('b9_rev.fastq','w+')],
            'TAGCCATG':[open('c9_fwd.fastq','w+'),open('c9_rev.fastq','w+')], 'CGGTAATC':[open('c3_fwd.fastq','w+'),open('c3_rev.fastq','w+')],
            'CTCTGGAT':[open('b3_fwd.fastq','w+'),open('b3_rev.fastq','w+')], 'TACCGGAT':[open('c4_fwd.fastq','w+'),open('c4_rev.fastq','w+')],
            'CTAGCTCA':[open('a11_fwd.fastq','w+'),open('a11_rev.fastq','w+')], 'CACTTCAC':[open('c7_fwd.fastq','w+'),open('c7_rev.fastq','w+')],
            'GCTACTCT':[open('b2_fwd.fastq','w+'),open('b2_rev.fastq','w+')], 'ACGATCAG':[open('a1_fwd.fastq','w+'),open('a1_rev.fastq','w+')],
            'TATGGCAC':[open('b7_fwd.fastq','w+'),open('b7_rev.fastq','w+')], 'TGTTCCGT':[open('a3_fwd.fastq','w+'),open('a3_rev.fastq','w+')],
            'GTCCTAAG':[open('b4_fwd.fastq','w+'),open('b4_rev.fastq','w+')], 'TCGACAAG':[open('a12_fwd.fastq','w+'),open('a12_rev.fastq','w+')],
            'TCTTCGAC':[open('c10_fwd.fastq','w+'),open('c10_rev.fastq','w+')], 'ATCATGCG':[open('a2_fwd.fastq','w+'),open('a2_rev.fastq','w+')],
            'ATCGTGGT':[open('c2_fwd.fastq','w+'),open('c2_rev.fastq','w+')], 'TCGAGAGT':[open('a10_fwd.fastq','w+'),open('a10_rev.fastq','w+')],
            'TCGGATTC':[open('b8_fwd.fastq','w+'),open('b8_rev.fastq','w+')], 'GATCTTGC':[open('a7_fwd.fastq','w+'),open('a7_rev.fastq','w+')],
            'AGAGTCCA':[open('b10_fwd.fastq','w+'),open('b10_rev.fastq','w+')], 'AGGATAGC':[open('a8_fwd.fastq','w+'),open('a8_rev.fastq','w+')],
            'fail':[open('noindex_fwd.fastq','w+'),open('noindex_rev.fastq','w+')]} #make index dictionary to hold file names
    with gzip.open(r1,'rt') as r1, gzip.open(r2,'rt') as r2, gzip.open(i1,'rt') as i1, gzip.open(i2,'rt') as i2: #open all index & read files for reading
        while r1 and r2 and i1 and i2: #iterate through both index & both read files
            h1=r1.readline().strip() #read 1 block
            if not h1: #exit when no more lines
                break
            seq1,opt1,qual1 = r1.readline().strip(),r1.readline().strip(),r1.readline().strip()
            lnct+=1 #read counter
            h2,seq2,opt2,qual2 = r2.readline().strip(),r2.readline().strip(),r2.readline().strip(),r2.readline().strip()
            h3,index1,opt3,quali1 = i1.readline().strip(),i1.readline().strip(),i1.readline().strip(),i1.readline().strip()
            h4,index2,opt4,quali2 = i2.readline().strip(),i2.readline().strip(),i2.readline().strip(),i2.readline().strip()
            index2=complement_seq(index2) #reverse complement and store index 2
            iqual1,iqual2 = eval_quality(quali1),eval_quality(quali2) #get avg index quality
            if iqual1 >= 30 and iqual2 >= 30: #set standard for index read quality
                index1,index2 = fix_n(index1,index2) #use other index to fix starting 'n' error
                if index1 in indexes:
                    if index1 == index2: #decided not to allow 1 off, alt code allows hamming distance of 1
                        match_ct+=1 #add to match count
                        indexes.get(index1)[0].write(h1+':'+index1+'\n'+seq1+'\n'+opt1+'\n'+qual1+'\n') #append entry to fwd index file
                        indexes.get(index1)[1].write(h2+':'+index1+'\n'+seq2+'\n'+opt2+'\n'+qual2+'\n') #append entry to rev index file
                    else:
                        hop_ct+=1 #register hopped index to count
                        index2='fail' #destination fail file
                else:
                    index2='fail'#sequencing error
            else: #quality check fails
                index2='fail' #destination fail file
            if index2=='fail': #reads marked for fail file
                indexes.get(index2)[0].write(h1+':'+index1+'\n'+seq1+'\n'+opt1+'\n'+qual1+'\n') #append entry to fwd fail file
                indexes.get(index2)[1].write(h2+':'+index1+'\n'+seq2+'\n'+opt2+'\n'+qual2+'\n') #append entry to rev fail file
    for key in indexes: #use index dict to close all open files
        indexes[key][0].close() #close fwd file
        indexes[key][1].close() #close rev file
    hop_rate=hop_ct/lnct #get index hopping rate
    return match_ct, hop_rate, len(indexes)
def main():
    '''run all fxns'''
    args = get_arguments() #uses argparse to get files
    matches,ihop_rate,index_dict = compare_fastq(args.read1,args.index1,args.index2,args.read2)
    print('Number of matches:',matches,'Index hopping rate:',ihop_rate,'\nIndex dictionary size:',index_dict)
    return None

main()

#files on talapas for use
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
