#!/usr/bin/env python

#SBATCH --partition=short
#SBATCH --job-name=demult_part1
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-user=annalundberg92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt


def get_arguments():
    parser = argparse.ArgumentParser(description="reads fasta file, stores gene id in dict")
    parser.add_argument("-f", "--filename", help="name of file",\
                        required=True, type=str)
    parser.add_argument("-t", "--type", help="index or seq",\
                        required=True, type=str)
    return parser.parse_args()

def convert_phred(letter):
    """(string)->int
    Converts a single character into a phred score"""
    phred = ord(letter)-33
    return phred

def populate_array(file, typ):
    '''(string, string)-> array
    Function opens fastq file, converts phred quality score, adds up q scores of
    each position, stores in array mean_scores, divides mean_scores by number
    of reads to return an array of mean phred scores per position'''
    p_scores = []
    ct = 0
    if typ == 'index':
        num = 8
    elif typ == 'seq':
        num = 101
    else:
        print('Not a valid type')
    while ct < num:
        p_scores.append(0.0)
        ct += 1
    with open(file) as f:
        LN = 0
        for line in f:
            LN+=1
            line = line.strip('\n')
            if LN%4 == 0:
                pos = 0
                for ch in line:
                    qscore = convert_phred(ch)
                    p_scores[pos]+= qscore
                    pos += 1
        mean_scores=np.array(p_scores)
        reads = LN//4
        mean_scores = mean_scores/reads
    return mean_scores, num

def plot_phred(phred_array, nscores, fname):
    '''(numpy array, int) -> filename
    This fxn takes an array representing mean phred scores per position, plots
    mean phred scores, saves plot as a file, which fxn returns name of'''
    plt.plot(range(nscores), phred_array, marker = 'o', linestyle= 'None')
    plt.title('Phred Quality Scores')
    plt.xlabel('Sequence position')
    plt.ylabel('Mean quality score')
    fig_name = fname+'plot.png'
    plt.savefig(fig_name)
    return fig_name

def main():
    '''run all the fxns'''
    args = get_arguments()
    mean_phred, num_scores = populate_array(args.filename, args.type)
    print(mean_phred)
    plot_phred(mean_phred, num_scores, args.filename)
    return None

main()
