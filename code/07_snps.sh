### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 07. Genotyping, SNPs and population genomics                             ###
### ======================================================================== ###


### Connect to HPC cluster
ssh <account>@carl.hpc.uni-oldenburg.de

# Account ids and passwords can be found on StudIP in Files | course_accounts.csv
# UOL HPC website: https://uol.de/en/school5/sc/high-perfomance-computing
# UOL HPC Wiki: https://wiki.hpcuser.uni-oldenburg.de/index.php?title=HPC_User_Wiki_2016


### Update course repository
cd meg23_repo
git pull

# first time (from ~): git clone https://github.com/mhelmkampf/meg23_repo.git



### ============================================================================
### Exercise 1: View and filter VCF file

### Navigate to "local" directory (create if needed)


### Copy file to "local" directory
cp ../meg23_repo/data/genome/snps_hamlets_lg12.vcf.gz .


### Decompress file


### View file
cat snps_hamlets_lg12.vcf


### View only first 10 lines


### Calculate stats
ml hpc-env/8.3 BCFtools/1.15.1-GCC-8.3.0   # load BCFtools module

bgzip snps_hamlets_lg12.vcf   # re-compress file with bgzip
tabix -p vcf snps_hamlets_lg12.vcf.gz   # index VCF file
bcftools stats snps_hamlets_lg12.vcf.gz


### VCFtools: command line tool to filter VCF files and run calculations on SNP data
ml VCFtools/0.1.16-GCC-8.3.0   # load VCFtools module

# Manual: https://vcftools.github.io/man_latest.html


### Filter sites (rows)
vcftools \
    --gzvcf snps_hamlets_lg12.vcf.gz \
    --max-missing 1 \
    --mac 2 \
    --recode \
    --stdout | bgzip > snps_hamlets_filtered.vcf.gz


### How many sites were retained after filtering?


### Original code to obtain example dataset (filters by sample and chromosome id)
# vcftools \
#   --gzvcf $DATA/shared/3_genotypes/chapter2_phased_mac2.vcf.gz \
#   --keep meg_36.ids \
#   --chr LG12 \
#   --recode \
#   --stdout | bgzip > snps_hamlets_lg12.vcf.gz



### ============================================================================
### Exercise 2: Calculate Fst along chromosome

### Calculate heterozygosity and Fis for each individual
vcftools \
  --gzvcf snps_hamlets_filtered.vcf.gz \
  --het \
  --stdout > Het_hamlets_snps.tsv


### Example parameters for other calculations
# --weir-fst-pop <file>   # Estimate Fst for populations defined in <file>
# --hardy                 # Reports a p-value for each site from a Hardy-Weinberg Equilibrium test


### Download files from cluster to local computer
# Open a new git bash windows on your local computer
# Navigate to "local" directory there

scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/snps_hamlets_filtered.vcf.gz .
scp <account>@carl.hpc.uni-oldenburg.de:/user/<account>/local/Het_hamlets_snps.tsv .



### ============================================================================
### Exercise 3: Plot PCA in R

### Open copy of this script in R


### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory


### Install / load packages
install.packages("vcfR")

library(vcfR)
library(adegenet)
library(tidyverse)


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
(l <- ggplot(data = scores_lg12, aes(x = PC1, y = PC2, color = Species, shape = Location)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "LG12") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5)
        )
)


### Plot eigenvalues (proportion of variance explained by each PC)
pca_lg12$eig

var <- pca_lg12$eig / sum(pca_lg12$eig)

barplot(var, main = "Proportion of variance explained", las = 2)




### ============================================================================
### Exercise 4: Plot per-site Fst

### Read in TSV file and add population information
het <- read_tsv("Het_hamlets_snps.tsv") %>%
  mutate(Species = str_sub(INDV, -6, -4),
         Location = str_sub(INDV, -3, -1),
         Population = str_sub(INDV, -6, -1))


### Summarize and visualize with boxplot
g <- ggplot() +
    geom_boxplot()



### ============================================================================
### Bonus material

### Distribution of allele frequencies
# freq_lg12 <- glMean(lg12_gl)   # compute mean of alternate alleles
# 
# hist(freq_lg12, breaks = 10, proba = TRUE, col = "grey",   # Shortcut for a quick histogram
#      xlab = "Allele frequencies",
#      main = "Distribution of alternate allele frequencies")



### ============================================================================
### Solutions

### Make local directory
mkdir local
cd local


### Decompress file
gzip -d *.gz


### View only first 10 lines
head snps_hamlets_lg12.vcf


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
