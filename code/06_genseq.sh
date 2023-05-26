### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 06. Genome sequencing and assembly                                       ###
### ======================================================================== ###


### Connect to HPC cluster
ssh <account>@carl.hpc.uni-oldenburg.de

# Account ids and passwords can be found on StudIP in Files | course_accounts.csv
# UOL HPC website: https://uol.de/en/school5/sc/high-perfomance-computing
# UOL HPC Wiki: https://wiki.hpcuser.uni-oldenburg.de/index.php?title=HPC_User_Wiki_2016


### Clone course repository
git clone https://github.com/mhelmkampf/meg23_repo.git


### Make and move to local directory


### Explore the cluster
sinfo -s   # view summary of cluster partitions / compute nodes

htop   # live view of processes on current node

module avail   # list available software modules



### ============================================================================
### Exercise 1: Compare Illumina and Pacbio reads

### Copy files to local directory
cp ../meg23_repo/data/genome/*.gz .


### Decompress files
gzip -d *.gz


### View file in Fastq format
cat HypPue1_illumina_raw_F.fastq

head -n 4 HypPue1_illumina_raw_F.fastq   # display first 4 lines of file 


### Count the number of sequences in file
# use grep -c 'pattern' file


### View and count number of sequences in PacBio reads


### How do the two types of read compare?



### ============================================================================
### Exercise 2: Calculate assembly metrics

### Copy assembly files to local
cp /nfs/data/haex1482/shared/teaching/HypPue1.1_illumina_scf.fas.gz .
cp /nfs/data/haex1482/shared/teaching/HypPue2.1_pacbio_pctg.fas.gz .


### Decompress and view assembly files


### Install assembly-stats v0.1.4 and calculate assembly stats
ml hpc-env/8.3 Python/3.7.4-GCCcore-8.3.0   # load Python module
pip install assembly_stats                  # use Python's pip to install

assembly_stats HypPue1.1_illumina_scf.fas
assembly_stats HypPue2.1_pacbio_pctg.fas


### Discuss the differences between the assemblies



### ============================================================================
### Exercise 3: Trim reads and assess read quality before and after

### Load modules
ml hpc-uniol-env
ml cutadapt/1.9.1
ml FastQC/0.11.5

# FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


### Remove adapters and low quality bases
cutadapt -h

cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \   # adapter sequence to be removed from 3' end
  -q 20 \                                  # minimum quality score at 3' end
  -m 120 \                                 # minimum length after trimming
  -o HypPue1_illumina_trimmed_F.fastq \
  HypPue1_illumina_raw_F.fastq


### Run FastQC, a quality control tool for high-throughput sequence reads
fastqc HypPue1_illumina_raw_F.fastq HypPue1_illumina_raw_R.fastq


### Re-run FastQC on trimmed reads


### Download and view QC reports in HTML format
scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/*.html .

# Note: Run from local computer, e.g. meg23_exercises/local



### ============================================================================
### Bonus: assembly code example (PacBio HiFi reads, hifiasm assembler)

# hifiasm \
#   -o GlaCem_hifiasm_v1 \
#   -t 12 \
#   $base/2_ccs/GlaCem?_ccs.fastq.gz

## Convert assembly graph to fasta
# awk '/^S/ { print ">"$2 ; print $3 }' GlaCem_hifiasm_v1.bp.p_ctg.gfa \
#   > GlaCem_hifiasm_v1.bp.p_ctg.fas


### Compare multiple assemblies with Quast: http://cab.cc.spbu.ru/quast/



### ============================================================================
### Solutions

### Make local directory
mkdir local
cd local


### Count the number of sequences in file
grep -c '@HWI' HypPue1_illumina_raw_F.fastq   # note: '@' is also found in quality score


### View and count number of sequences in PacBio reads
head -n 4 HypPue2_pacbio_ccs.fastq
grep -c '@m54' HypPue2_pacbio_ccs.fastq


### Decompress and view assembly files
gzip -d *.gz
head HypPue1.1_illumina_scf.fas
head HypPue2.1_pacbio_pctg.fas


### Re-run FastQC on trimmed reads
fastqc HypPue1_illumina_trimmed_F.fastq