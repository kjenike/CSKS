#!/usr/bin/env python3

import argparse
import sys
from pyfaidx import Fasta
from random import choice
from random import seed
from random import randint

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]

    parser = argparse.ArgumentParser(description="Create a randomised subset of a genome.")
    parser.add_argument('-g', '-genome', help='genome file (.fasta, can be gzipped, but must be indexed)')
    parser.add_argument('-n', '-number_of_chromosomes', help="Number of chromosomes (of length -l) to be generated (defalt: 1)", default = 1, type = int)
    parser.add_argument('-l', '-length', help='The total length of newly sampled reference (in Mbp, default: 20)', default = 1, type = int)
    parser.add_argument('-o', '-output', help='output pattern (default: sampled_genome)', default = 'sampled_genome.fasta')
    parser.add_argument('-s', '-seed', help='seed for generating random numbers (default: defined by time)', default = None, type = int)
    args = parser.parse_args(args)

    ##### Setting up the random number generator
    if args.s == None:
        used_seed = randint(100000, 99999999)
    else:
        used_seed = args.s

    seed(used_seed)
    sys.stderr.write('This sampling is can be regenerated using seed {} (parameter -s)\n'.format(used_seed))

    ##### Now we load the reference genome index and generate the sub sampled genome
    genome = Fasta(args.g)
    list_of_scaffolds = list(genome.keys())
    sys.stderr.write('loaded {} genome file\n'.format(args.g))

    # load the index?

    chromosome = 0
    total_desired = int(args.l * 1e6)

    with open(args.o, 'w') as output_fasta:
        while chromosome < args.n:
            chromosome += 1
            output_fasta.write('>simulated_chromosome_' + str(chromosome) + '_with_' + str(used_seed) + '_seed\n')

            sampled_length = 0
            while sampled_length < total_desired:
                picked_scf = choice(list_of_scaffolds)
                scaffold_record = genome[picked_scf]
                # use the INDEX instead!!!
                scaffold_length = len(scaffold_record)

                if scaffold_length < 10000:
                    continue
                else:
                    picked_pos = randint(0, scaffold_length - 10000)

                if 10000 + sampled_length > total_desired:
                    seq_to_print = scaffold_record[picked_pos:(picked_pos + (total_desired - sampled_length))]
                else:
                    seq_to_print = scaffold_record[picked_pos:(picked_pos + 10000)]

                if str(seq_to_print).count('N') / len(seq_to_print) > 0.05:
                    continue

                sampled_length += len(scaffold_record)
                output_fasta.write(str(seq_to_print))
            output_fasta.write('\n')

    sys.stderr.write('Done\n')
