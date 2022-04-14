### The actual power analysis

Here I am designing individual components of the pipeline. In the end I will wrap everything in a single script so each job will do all the tasks on local scatch disc,

#### Simulating the reference

Script `scripts/generate_reference_subset.py` can generate a subset of the genome reference with desired total length. They will be randomly sampled scaffolds that are already in the genome. And it was tested to have deterministic behavior when seed is specified. Example use would be

```
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "A" -l 1 -n 19 -o simulation/Afus_sim_reference_A.fasta
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "X" -l 1 -o simulation/Afus_sim_reference_X.fasta

cat simulation/sim_reference_A.fasta simulation/sim_reference_X.fasta > data/sim_reference_complete.fasta
```

which generates 1Mbp of X-linked reference sequences.

#### Simulating variants on top of the refenrece

```
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_A.fasta -het 0.003 -o simulation/sim_genome_A
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_X.fasta -het 0.003 -o simulation/sim_genome_X
```

#### Simulating the reads

```
# maternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_mat.fasta  -c 50 -r 150 -o simulation/maternal_A

# maternal X
python3 scripts/simulate_reads.py -g simulation/sim_genome_X_mat.fasta -c 50 -r 150 -o simulation/maternal_X

# paternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_pat.fasta  -c 50 -r 150 -o simulation/paternal_A

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

```{R}
tab_mat_A <- read.table('simulation/kmcdb_k21_maternal.hist')
tab_pat_A <- read.table('simulation/kmcdb_k21_paternal.hist')
tab_mat_X <- read.table('simulation/kmcdb_k21_maternal_X.hist')

tab <- read.table('simulation/kmcdb_k21.hist')

par(mfrow = c(2, 2))
plot(tab$V1[4:100], tab$V2[4:100])
plot(tab_mat_A$V1[4:100], tab_mat_A$V2[4:100])
plot(tab_pat_A$V1[4:100], tab_pat_A$V2[4:100])
plot(tab_mat_X$V1[4:100], tab_mat_X$V2[4:100])
# that all looks alright!

x <- tab_mat_A$V1[4:100]
y <- tab_mat_A$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 19e6, bias = 1))


x <- tab_pat_A$V1[4:100]
y <- tab_pat_A$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 19e6, bias = 1))


x <- tab_mat_X$V1[4:100]
y <- tab_mat_X$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 1e6, bias = 1))
```

but let's see about two tissue model

```{R}
Rscript scripts/two_tissue_model.R -i simulation/kmcdb_k21.hist
```

#### Mapping reads and estimating the coverage

```
REFERENCE=simulation/generated_complete_sim_reference
bowtie2-build simulation_complete_sim_reference.fasta $REFERENCE

R1=simulation/sim_reads_R1.fq
R2=simulation/sim_reads_R2.fq

bowtie2 --very-sensitive-local -p 4 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id 1 --rg SM:rand --rg PL:ILLUMINA --rg LB:LIB_rand \
        | samtools view -h -q -20 \
        | samtools sort -O bam - > simulation/rand.rg.sorted.bam
```

Mapped reads to a table of coverages (coverage per 10,000 nt).

```
samtools depth simulation/rand.rg.sorted.bam | perl scripts/depth2windows.pl 10000 > simulation/per_window_coverage.tsv
```

```
Rscript scripts/mapping_coverage_tab2cov_estimates.R -i simulation/per_window_coverage.tsv -o simulation/test
```
### Explored parametric space

These following values are what I consider a reasonable parametric space to explore. It's totalling ~300 simulations, each takes about 30'. The scripts require specified maternal and paternal coverage (rather than maternal coverage and fraction of sperm), so that will be handles internally in the R script bellow.


```
heterozygosity = 0.01, 0.1, 0.25, 0.5 (autosomal average heterozygosity in %)
X_chromosomes = 1 2 5 10 (Mbp out of 20 Mbp reference)
fraction_of_sperm = 0, 0.01, 0.05, 0.10, 0.25, 0.50 (0 is practically non-PGE system)
sequencing_depth = 10x, 15x, 25x (1n coverage)
```

The empirical data suggest that the springtail will be the closest to the simulation "het = 0.25, X = 10, fraction_of_sperm = 0.25 and sequencing_depth = 15x". Note that while we were able to manually curate the data for the three focal species, these simulations have the models fitted with naively estimated priors from the data. We have not dedicated extensive effort in making this method complete solution, that will be a task for a future study.

### Running it on cluster

The easiest way to reproduce the analysis is using conda, lets create an enviorment `CSKS` (evironment tested on Ubuntu 18.04 and Debian GNU/Linux 9).

```
conda create -n CSKS -c bioconda -c conda-forge msprime kmc python=3 numpy matplotlib wgsim pyfaidx samtools bowtie2 freebayes
# I had to add a few package later on:
# conda install -c bioconda wgsim
# conda install -c bioconda pyfaidx
# conda install -c bioconda samtools
# conda install -c bioconda bowtie2
# conda install -c bioconda freebayes
conda install -c r r
```

Now open `R` and install the required R packages

```
# run within R
install.packages('argparse')
install.packages('minpack.lm')
install.packages('nlstools')
```

now quit R again. With this newly created conda enviorment we should be able to test one replicate

```
qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 2 -b yes "scripts/run_power_analysis_replicate.sh 0.001 10 25 15"

qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 4 -b yes "scripts/run_power_analysis_replicate.sh 0.001 1 25 20"
```

That works. We wronte a script [scripts/create_list_of_commands.R](scripts/create_list_of_commands.R) that generates a long list of commands for all the parameter combinations we are interested in.

```bash
Rscript scripts/create_list_of_commands.R
```

and now we can run all the commands. Note this is 288 submited jobs to the cluster, so check with your sysadmins that this is something ok to do. You will also need to adjust the cluster submission command (`qsub ...`). I will do that using GNU parallel.

```
# conda activate CSKS
conda install -c conda-forge parallel

parallel -j 1 'qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 4 -b yes -l h="bigbang" {}' :::: scripts/power_analysis_commands.sh
```

### post-simulation adjustments to the model

I originally kept one point the error peak to give have `./\/\` kind of distributions (thinking that `.` will never mess it up), but it does, namely when the coverage is low and the disr looks more like this `./-\` (two peaks are blended). So I will remove this one small detail (TODO: add commit url) and rerun all the morels - I inteligently saved all the kmer spectra, so I don't need to rerun the whole analysis.

```bash
for sim in data/simulations/*; do
    cp scripts/two_tissue_model* $sim/scripts
    cd $sim
    Rscript scripts/two_tissue_model.R -i output/kmcdb_k21.hist -o output/two_tissue_model
    cd -
done
```


And I wanted to pull all the figures to one place.

```bash
mkdir -p figures/simulated_two_tissue_models
for sim in data/simulations/*; do sim=$(basename $sim); ln -s `pwd`/data/simulations/$sim/output/two_tissue_model_plot.png figures/simulated_two_tissue_models/"$sim".png ; done
```

### Explororing exceptions

I found that in many cases of very low coverage, my automatic way of guessing priors made them too high and the model converged on very wild values. Here is one example

```
data/simulations/sim_95407_h0.0030_X5_m15_p7.5
```

So I adjusted those models (~20) by adding

```
if fraction_of_sperm < -0.2 or > 0.7:
  rerun the model with 0.7 of the original prior
if still within those crazy boundaries:
  rerun the model with 0.5 of the original prior
```

which resolved all the cases (and would not affect any of the other estimates).
