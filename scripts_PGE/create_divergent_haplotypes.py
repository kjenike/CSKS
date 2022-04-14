import msprime
import argparse
from pyfaidx import Fasta
from random import sample
from random import seed
from random import randint
from statistics import mean
from sys import argv

if __name__ == "__main__":
    args = None
    if len(argv) == 1:
        args = ["--help"]

    parser = argparse.ArgumentParser(description="Create reasobnably distributed heterozygosity along dataset on an genome input.")
    parser.add_argument('-g', '-genome', help='indexed genome used as the reference template (.fasta)')
    parser.add_argument('-het', '-heterozygosity', help='total simulated heterozygosity (theta, default: 0.003)', default = 0.003, type = float)
    parser.add_argument('-Ne', '-population_size', help='Ne (default: 1e4)', default = 1e4, type = int)
    parser.add_argument('-o', '-output', help='output pattern. "_mat.fasta" and "_pat.fasta" files will be generated (default:basename of the provided reference)', default = None)
    parser.add_argument('-s', '-seed', help='seed for generating random numbers (default: defined by time)', default = None, type = int)
    parser.add_argument('--monosomic', help='will generate only one (_mat) haplotype (rather than _mat and _pat)', dest='monosomic', action='store_const', const = True, default = False,)
    args = parser.parse_args(args)

    ##### Setting up the random number generator
    if args.s == None:
        used_seed = randint(100000, 99999999)
    else:
        used_seed = args.s
    # used_seed = 15641

    # calculating mutation rate (mu) from heterozygosity (theta) and population size (ne)
    # heterozygosity (theta) = 4 * mu * Ne
    # therefore mu = theta / (4 * Ne)
    # args.het = 0.003
    mu = args.het / (4 * args.Ne)
    sequence_headers = []
    maternal_sequences = []
    paternal_sequences = []

    # args.g = 'data/generated/sim_reference_X.fasta'
    fasta_records = Fasta(args.g)
    for single_record in fasta_records:
        seq_header = single_record.name
        sequence_headers.append(seq_header)
        number_of_loci = single_record.unpadded_len

        # sys.stderr.write('loaded {} file\n'.format(args.g))

        # this simulates ancestry along all the loci on the chromosome given a reasonable recombination rates (50cM over the whole chromosome)
        ts = msprime.sim_ancestry(
                samples=1,
                ploidy=2,
                population_size=args.Ne,
                discrete_genome=True,# default setting now?
                recombination_rate=(0.5 / number_of_loci), # 0.5 / 1e6
                sequence_length=number_of_loci, random_seed=used_seed)
        mts = msprime.sim_mutations(ts, rate=mu, random_seed=(used_seed + 1))


        # setting up template seuquence
        chromosome_seq = fasta_records[seq_header][0:]
        ref = list(chromosome_seq.seq)
        alt = list(chromosome_seq.seq)

        for variant_locus in mts.variants():
            pos = int(variant_locus.site.position)
            pool = set(['A', 'T', 'C', 'G'])
            ancestral = str(chromosome_seq[pos])
            try: # this is catching only ancestral = N (ancestral sequence is generated out of real sequenced nucleotides)
                pool.remove(ancestral) # remove the ancestral allele out of the pool
            except KeyError:
                ancestral = 'N'
            mutated_variants = [ancestral] + sample(list(pool), 3) # randomly generate the two derived states
            ref[pos] = mutated_variants[variant_locus.genotypes[0]] # save the states in the ref and alt lists of nucleotides
            alt[pos] = mutated_variants[variant_locus.genotypes[1]]

        maternal_sequences.append(''.join(ref))
        paternal_sequences.append(''.join(alt))

#  the next section needs to be updated for all the fasta
    if args.o == None:
        output_pattern = ''.join(args.g.split('.')[:-1])
    else:
        output_pattern = args.o

    with(open(output_pattern + '_mat.fasta', 'w')) as mat_file:
        for i,seq in enumerate(maternal_sequences):
            mat_file.write('>' + sequence_headers[i] + '_mutated_with_' + str(used_seed) + '_seed_mat\n')
            mat_file.write(seq + '\n')

    if not args.monosomic:
        with(open(output_pattern + '_pat.fasta', 'w')) as pat_file:
            for i,seq in enumerate(maternal_sequences):
                pat_file.write('>' + sequence_headers[i] + '_mutated_with_' + str(used_seed) + '_seed_pat\n')
                pat_file.write(seq + '\n')

# ref = maternal_sequences[0]
# alt = paternal_sequences[0]
# mean([ref[i] == alt[i] for i in range(number_of_loci)])
