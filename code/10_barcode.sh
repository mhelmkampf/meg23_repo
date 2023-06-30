### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 10. DNA barcoding                                                        ###
### ======================================================================== ###


### Preparations

# On Windows, install and launch BioEdit
# Zipped installation file provided in meg23_repo/other
# Or download from https://bioedit.software.informer.com/

# On Mac, launch 4Peaks
# Executable provided in meg23_repo/other

# Copy trace files to local working directory
cp meg23_repo/data/barcode/*.ab1 .



### ============================================================================
### Exercise 1: Extract COI barcodes from Sanger trace files

### Background – PCR / sequencing primer cocktail used for COI:

# VF2_t1    tgtaaaacgacggccagtCAACCAACCACAAAGACATTGGCAC
# FishF2_t1 tgtaaaacgacggccagtCGACTAATCATAAAGATATCGGCAC

# FishR2_t1 caggaaacagctatgacACTTCAGGGTGACCGAAGAATCAGAA
# FR1d_t1   caggaaacagctatgacACCTCAGGGTGTCCGAARAAYCARAA

# Amplifies approx. 650 bp fragment in 5' region of COI
# Primers are M13-tailed to allow DNA sequencing
# (Ivanova et al. 2007, Molecluar Ecology Notes, lower case = tail)


### View trace files

# Open trace files in BioEdit (Windows) or 4Peaks (Mac)
# Note that each sample was sequenced in forward and reverse direction: seq_*F.ab1, seq_*R.ab1
# How many positions at the beginning (5') and end (3') look unreliable?


### Trim low quality 5' and 3' positions

# BioEdit: select range of positions to keep, copy selection (Ctrl + C) to plain text file
# 4Peaks: select range of positions to keep, crop (Cmd + K) and export as plain text file
# Save as seq_*F.tr, seq_*R.tr


### Reverse complement trimmed R sequence

# http://www.reverse-complement.com/
# Reverse complement seq_*R.tr file
# Save as seq_*R.rc


### Align sequences using MAFFT

# https://mafft.cbrc.jp/alignment/server/
# Paste seq_*F.tr and seq_*R.rc into input field and submit
# Go to "Fasta format" save downloaded file as seq_*.aln


### Create consensus sequence with Cons

# https://www.ebi.ac.uk/Tools/msa/emboss_cons/
# Select DNA and submit seq_*.aln alignment
# Save to file as seq_*.con
# Ideally, change sequence name in header to seq_*
# This is your barcode



### ============================================================================
### Exercise 2: Identify species based on COI barcode with BOLD

# Barcode of Life Data System (BOLD): http://www.boldsystems.org/index.php/
# Select Identification, Animal Identification (COI), Species Level Barcode Records
# Enter barcode sequence and submit

# How confident are you in the identification based on similarity score, 
# and the gap between within-BIN and NN distances (see BIN page)?

# Optinally, explore species and tree pages



### ============================================================================
### Bonus material: Process sequence files on cluster using bash loops

### Connect to HPC cluster
ssh <account>@carl.hpc.uni-oldenburg.de
# Account ids and passwords can be found on StudIP in Files | course_accounts.csv


### Update course repository
cd meg23_repo
git pull
# first time (from ~): git clone https://github.com/mhelmkampf/meg23_repo.git


### Copy sequence text files to local working directory
cd ../local
rm *
cp ../meg23_repo/data/barcode/*.txt .


### Trim low quality 5' and 3' positions using the Fastx toolkit
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
rm *.bz2

for i in *.txt
do
  bin/fastx_trimmer \
  -f 50 \
  -l 550 \
  -i ${i} \
  -o ${i%.txt}.tr
done


### Reverse complement trimmed R sequence with EMBOSS revseq
ml EMBOSS/6.6.0

for i in *R.tr
do
  revseq \
  ${i} \
  ${i%.tr}.rc
done


### Align sequences using MAFFT
ml hpc-env/8.3
ml hpc-env/8.3 MAFFT/7.475-GCC-8.3.0-with-extensions

for i in *F.txt
do
  s=${i%F.txt}
  mafft <(cat ${s}F.tr ${s}R.rc) > ${s}.aln
done


### Create consensus sequence with EMBOSS Cons
ml hpc-uniol-env
ml EMBOSS/6.6.0

for i in *.aln
do
  cons \
  -sequence ${i} \
  -outseq ${i%.aln}.con \
  -name ${i%.aln}
done


### Compile consensus sequences in one file
cat *.con > co1.fas



### ----------------------------------------------------------------------------
### Bonus material: identify barcode / species with BLAST and the NCBI nt database

# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# Select Nucleotide BLAST, Nucleotide collection (nr/nt), megablast



### ============================================================================
### Solutions

### Species ID of samples

# 10: Oncorhynchus keta (Chum salmon / Keta-Lachs), 100%, large barcode gap
# 11: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, small barcode gap
# 14: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., 100%, no barcode gap
# 15: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, see above, small barcode gap
# 19: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, see above, small barcode gap
# 20: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., 100%, small barcode gap
# 22: Gadus morhua (Atlantic cod / Kabeljau), 100%, small barcode gap
# 26: Gadus morhua (Atlantic cod / Kabeljau), 100%, small barcode gap
# 27:
# 28:
# 35: