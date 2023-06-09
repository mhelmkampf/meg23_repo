### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 08. Recombination and linkage disequilibrium                             ###
### ======================================================================== ###

### Open a copy of this script in R

### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory



### ============================================================================
### Exercise 1: Plot SNPs as PCA in R

### Install / load packages
install.packages("vcfR")

library(vcfR)
library(adegenet)
library(tidyverse)


# We filtered and downloaded a VCF file with hamlet SNPs to "local" last time:

# vcftools \
# --gzvcf snps_hamlets_lg12.vcf.gz \
# --max-missing 1 \
# --mac 2 \
# --recode \
# --stdout | bgzip > snps_hamlets_filtered.vcf.gz

# scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/snps_hamlets_filtered.vcf.gz .

# The file can also found be in meg23_repo/data/genome


### Read VCF file into R
lg12_vcf <- read.vcfR("snps_hamlets_filtered.vcf.gz")
lg12_vcf


### Convert from vcfR to genlight object
lg12_gl <- vcfR2genlight(lg12_vcf)


### Principal Component Analysis (PCA)
pca_lg12 <- glPca(lg12_gl, nf = 2)
pca_lg12
pca_lg12$scores   # view Principal Components (n = 2)


### Convert to tibble, add species and location information
scores_lg12 <- as.data.frame(pca_lg12$scores) %>%
  rownames_to_column("Sample") %>%
  as_tibble() %>%
  mutate(Species = str_sub(Sample, -6, -4),
         Location = str_sub(Sample, -3, -1)
         )


### Plot PCA
p <- ggplot(data = scores_lg12, aes(x = PC1, y = PC2, color = Species, shape = Location)) +
  geom_point(size = 4, alpha = 0.75)

p +
  labs(title = "LG12") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5)
        )


### Plot eigenvalues (proportion of variance explained by each PC)
pca_lg12$eig

var <- pca_lg12$eig / sum(pca_lg12$eig)

barplot(var, 
        main = "Proportion of variance explained", 
        las = 2, 
        ylim = c(0, 0.1))



### ============================================================================
### Exercise 2: Plot per-population Fis from genome-wide Fis per individuum

# We calculated heterozygosity and Fis from the hamlet SNPs last time:

# vcftools \
# --gzvcf snps_hamlets_filtered.vcf.gz \
# --het \
# --stdout > Het_hamlets_snps.tsv

# scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/Het_hamlets_snps.tsv .

# The file can also found be in meg23_repo/data/genome


### Read in TSV file and add population information
het <- read_tsv("Het_hamlets_snps.tsv") %>%
  mutate(Species = str_sub(INDV, -6, -4),
         Location = str_sub(INDV, -3, -1),
         Population = str_sub(INDV, -6, -1))


### Summarize and visualize with boxplot
f <- ggplot(het, aes(x = Population, y = F, fill = Species)) +
  geom_boxplot(color = "grey20",
               alpha = 0.75,
               lwd = 0.3)

f +
    scale_fill_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
    labs(title = NULL,
         x = "Population",
         y = "Mean genome-wide Fis") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(angle = 35)
    )



### ============================================================================
### Exercise 3: Filter VCF by linkage

### Connect to HPC cluster from a terminal (e.g. git bash)
ssh <account>@carl.hpc.uni-oldenburg.de

# Account ids and passwords can be found on StudIP in Files | course_accounts.csv
# UOL HPC website: https://uol.de/en/school5/sc/high-perfomance-computing
# UOL HPC Wiki: https://wiki.hpcuser.uni-oldenburg.de/index.php?title=HPC_User_Wiki_2016


### Update course repository
cd meg23_repo
git pull

# first time (from ~): git clone https://github.com/mhelmkampf/meg23_repo.git


### BCFtools: command line tool to call SNPs and manipulate VCF/BCF files
ml hpc-env/8.3 BCFtools/1.15.1-GCC-8.3.0   # load BCFtools module

# Manual: https://samtools.github.io/bcftools/bcftools.html


### Check BCFtools prune options
bcftools +prune


### Prune by a physical distance that guarantees loci are independent (compare to decay plot)
bcftools +prune snps_hamlets_filtered.vcf.gz \
  -n 1 -w 2kb \
  -Oz -o snps_hamlets_2kb.vcf.gz


### Prune by LD (r2)
bcftools +prune snps_hamlets_filtered.vcf.gz \
  -m 0.4 \   # default window size is 100kb
  -Oz -o snps_hamlets_04r.vcf.gz


### Calculate r2 and D statistics with VCFtools
ml VCFtools/0.1.16-GCC-8.3.0

# Manual: https://vcftools.github.io/man_latest.html

vcftools --gzvcf snps_hamlets_filtered.vcf.gz \
  --hap-r2 \
  --ld-window-bp 10000 \
  --stdout | head -n 10000 > LD_snps_hamlets_filtered.tsv

vcftools --gzvcf snps_hamlets_2kb.vcf.gz \
  --hap-r2 \
  --ld-window-bp 10000 \
  --stdout | head -n 10000 > LD_snps_hamlets_2kb.tsv

vcftools --gzvcf snps_hamlets_04r.vcf.gz \
  --hap-r2 \
  --ld-window-bp 10000 \
  --stdout | head -n 10000 > LD_snps_hamlets_04r.tsv


# scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/LD* .


### Plot r2 and D statistics

## Read in stats from vcftools output
un <- read_tsv("LD_snps_hamlets_filtered.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Un") %>%
  select(r2, Dprime, Set)

kb <- read_tsv("LD_snps_hamlets_2kb.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Kb") %>%
  select(r2, Dprime, Set)

ld <- read_tsv("LD_snps_hamlets_04r.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Ld") %>%
  select(r2, Dprime, Set)


## Plot r2
boxplot(un$r2, kb$r2, ld$r2, names = c("Un", "Kb", "Ld"), ylab = "r2")   # quick boxplot
boxplot(un$r2, kb$r2, ld$r2, names = c("Un", "Kb", "Ld"), ylab = "r2", outline = FALSE)   # without outliers


## Exercise: Plot D'


## Optional: Plot r2 using ggplot




### ============================================================================
### Exercise 4: Estimate Ne using heterozygote excess

### Unzip NeEstimator.zip in meg23_repo/other and move new program directory to local
### There, execute Ne2-1.exe (or Ne2M on Mac) through terminal or double-click

### Provide path to input file puella_caribbean.gen when prompted
### For method choose linkage disequilibrium
### Provide output file name (e.g. NeLD_caribbean.txt)
### Results will be saved in local/NeEstimator


### Compare these results to the results based on heterozygote excess (see 03_drift.R)




### ============================================================================
### Bonus material

### Test for linkage disequilibrium between microsatellite loci (G-test)
library(genepop)

test_LD("../meg23_repo/data/msats/puella_caribbean.txt", outputFile = "LD_msats_caribbean.txt")
test_LD("../meg23_repo/data/msats/hamlets_caribbean.txt", outputFile = "LD_msats_hamlets.txt")

# Null hypothesis: loci are independent



### ============================================================================
### Solutions



### How many sites were retained after filtering?
tabix -p vcf snps_hamlets_filtered.vcf.gz
bcftools stats snps_hamlets_filtered.vcf.gz


### Summarize and visualize with boxplot
(g <- ggplot(het, aes(x = Population, y = F, fill = Species)) +
    geom_boxplot(color = "grey20",
                 alpha = 0.75,
                 lwd = 0.3) +
    scale_fill_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
    labs(title = NULL,
         x = "Population",
         y = "Mean genome-wide Fis") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(angle = 35)
    )
)
