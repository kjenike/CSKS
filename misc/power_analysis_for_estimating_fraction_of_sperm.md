### The actual power analysis

Here I am designing individual components of the pipeline. In the end I will wrap everything in a single script so each job will do all the tasks on local scatch disc,

#### Simulating the reference

Script `scripts/generate_reference_subset.py` can generate a subset of the genome reference with desired total length. They will be randomly sampled scaffolds that are already in the genome. And it was tested to have deterministic behavior when seed is specified. Example use would be

```
python scripts_PGE/generate_reference_subset.py -g data/genome_short_headers.fa -l 1 -n 19 -o simulation/Afus_sim_reference_A.fasta
python scripts_PGE/generate_reference_subset.py -g data/genome_short_headers.fa -l 1 -o simulation/Afus_sim_reference_X.fasta

cat simulation/sim_reference_A.fasta simulation/sim_reference_X.fasta > data/sim_reference_complete.fasta
```

which generates 1Mbp of X-linked reference sequences.

#### Simulating variants on top of the refenrece

```
python3 scripts_PGE/create_divergent_haplotypes.py -g simulation/sim_reference_A.fasta -het 0.003 -o simulation/sim_genome_A
python3 scripts_PGE/create_divergent_haplotypes.py -g simulation/sim_reference_X.fasta -het 0.003 -o simulation/sim_genome_X
```

#### Simulating the reads

```
# maternal A
python3 scripts_PGE/simulate_reads.py -g simulation/sim_genome_A_mat.fasta  -c 50 -r 150 -o simulation/maternal_A

# maternal X
python3 scripts_PGE/simulate_reads.py -g simulation/sim_genome_X_mat.fasta -c 50 -r 150 -o simulation/maternal_X

# paternal A
python3 scripts_PGE/simulate_reads.py -g simulation/sim_genome_A_pat.fasta  -c 50 -r 150 -o simulation/paternal_A

cat simulation/maternal_A_R1.fq simulation/maternal_X_R1.fq simulation/paternal_A_R1.fq > simulation/sim_reads_R1.fq
cat simulation/maternal_A_R2.fq simulation/maternal_X_R2.fq simulation/paternal_A_R2.fq > simulation/sim_reads_R2.fq
```

#### K-mer based estimate of 1n and 2n

```
mkdir -p tmp

# inidividual subpart libraries
ls simulation/maternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_maternal.hist -cx1000

ls simulation/maternal_X*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_maternal_X.hist -cx1000

ls simulation/paternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_paternal.hist -cx1000

# joint histogram
ls simulation/maternal_A*.fq simulation/maternal_X*.fq simulation/paternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21.hist -cx1000
```

Alright, now we generated from the simualted data a k-mer histogram!
