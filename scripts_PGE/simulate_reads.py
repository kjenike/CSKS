import argparse
from random import seed
from random import randint
from statistics import mean
from pyfaidx import Fasta
from sys import argv
from os import system

if __name__ == "__main__":
    args = None
    if len(argv) == 1:
        args = ["--help"]

    parser = argparse.ArgumentParser(description="Create sequencing reads using wgsim simulator.")
    parser.add_argument('-g', '-genome', help='the reference genome that will be in silico "sequenced" (.fasta)')
    parser.add_argument('-c', '-coverage', help='expected sequencing coverage (default: 20x)', default = 20, type = float)
    parser.add_argument('-r', '-read_length', help='read length (default: 150bp)', default = 150, type = int)
    parser.add_argument('-e', '-error_rate', help='the probability of a sequencing error', default = 0.01, type = float)
    parser.add_argument('-o', '-output', help='output pattern. "_R1.fq" and "_R2.fq" files will be generated (default:basename of the provided reference)', default = None)
    parser.add_argument('-s', '-seed', help='seed for generating random numbers (default: defined by time)', default = None, type = int)
    args = parser.parse_args(args)

    fasta_records = Fasta(args.g)
    number_of_loci = 0
    for single_record in fasta_records:
        number_of_loci += single_record.unpadded_len

    ##### Setting up the random number generator
    if args.s == None:
        used_seed = randint(100000, 99999999)
    else:
        used_seed = args.s
    # used_seed = 15641

    simulated_reads = round((number_of_loci * args.c) / (2 * args.r))

    R1_name = args.o + "_R1.fq"
    R2_name = args.o + "_R2.fq"

    system_call = "wgsim -N {} -1 {} -2 {} -r 0 -S {} -e {} {} {} {}".format(simulated_reads, args.r, args.r, used_seed, args.e, args.g, R1_name, R2_name)
    system(system_call)
