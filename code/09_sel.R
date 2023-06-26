### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 09. Selection and Mutation                                               ###
### ======================================================================== ###

### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory



### ============================================================================
### Exercise 1: Plot r2 and D statistics from LD-filtered VCF files

# Previously, SNPs were filtered by physical distance of 2 kb, 
# by a maximum correlation coefficient of linkage disequilibrium (r2) of 0.4,
# or left unfiltered

# LD statistics were computed with BCFtools and downloaded to the local directory

# These output files can also be found in meg23_repo/data/genome
# Copy to "local": cp ../meg23_repo/data/genome/LD* .


### Install / load packages
library(tidyverse)


### Read in stats from VCFtools output
unf <- read_tsv("../meg23_repo/data/genome/LD_snps_hamlets_filtered.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Unf")

kb <- read_tsv("../meg23_repo/data/genome/LD_snps_hamlets_2kb.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Kb")

ld <- read_tsv("../meg23_repo/data/genome/LD_snps_hamlets_04r.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Ld")


### Plot r2
boxplot(unf$r2, kb$r2, ld$r2, 
        names = c("Unf", "Kb", "Ld"), 
        ylab = "r2",
        ylim = c(0, 0.1),
        outline = FALSE)   # without outliers


### Plot D'
boxplot(unf$Dprime, kb$Dprime, ld$Dprime, 
        names = c("Unf", "Kb", "Ld"), 
        ylab = "Dprime")


### Add zero line for interpretation
abline(h = 0, col = "red")



### ============================================================================
### Exercise 2: Simulating the effect of selection and mutation on allele frequencies in finite populations

# Go to: http://evolutiongenetics.georgetown.edu/simulations/driftselection

# Set effective population size to 200, switch on Natural Selection

## Choose a combination of Natural Selection parameters for the following scenarios:

# - Selection for a dominant phenotype / against a recessive phenotype

# - Selection against a dominant phenotype

# - Heterozygote advantage

## How do population size and the strength of natural selection influence the process?

# Switch on Mutation (with a -> A e.g. 0.0001), and set initial allele frequency (fA) to 0

## Simulate the trajectory of a deleterious and beneficial mutation



### ============================================================================
### Exercise 3: Apply genome-wide scans and identify Fst outlier loci  

### Fst calculations were done with vcftools, output and log files are found in meg23_repo/data/genome)

# Calculate joint Fst per SNP
# vcftools --gzvcf snps_hamlets_filtered.vcf.gz \
# --weir-fst-pop pop.gumhon.txt \
# --weir-fst-pop pop.indbel.txt \
# --weir-fst-pop pop.nigbel.txt \
# --weir-fst-pop pop.nighon.txt \
# --weir-fst-pop pop.puebel.txt \
# --weir-fst-pop pop.unibel.txt \
# --stdout 1> Fst_snp_lg12.tsv 2> Fst_snp_lg12.log

# Calculate joint Fst in sliding windows of 50 kb
# vcftools --gzvcf snps_hamlets_filtered.vcf.gz \
# --weir-fst-pop pop.gumhon.txt \
# --weir-fst-pop pop.indbel.txt \
# --weir-fst-pop pop.nigbel.txt \
# --weir-fst-pop pop.nighon.txt \
# --weir-fst-pop pop.puebel.txt \
# --weir-fst-pop pop.unibel.txt \
# --fst-window-step 5000 \
# --fst-window-size 50000 \
# --stdout 1> Fst_50kb_lg12.tsv 2> Fst_50kb_lg12.log


### Plot joint Fst, per SNP
snp <- read_tsv("../meg23_repo/data/genome/Fst_snp_lg12.tsv") %>%
  rename("FST" = "WEIR_AND_COCKERHAM_FST")

(s <- ggplot(data = snp, aes(x = POS, y = FST)) +
  geom_point(size = 0.25, alpha = 0.5) +
  labs(x = "Position", y = "Fst") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
)


### Plot joint Fst, along 50 kb sliding windows (use start bin and weighted Fst as axes)
win <- read_tsv("../meg23_repo/data/genome/Fst_50kb_lg12.tsv")

(w <- ggplot()
)


### Add chromosome-wide mean Fst (value reported in log files)
chrwide_fst <- 

w + geom_hline(yintercept = chrwide_fst, color = "blue")


### Identify 99.5% quantile
quantile(win$WEIGHTED_FST, probs = 0.995)

out <- win %>%
  mutate(OUTLIER = case_when(
    WEIGHTED_FST > quantile(WEIGHTED_FST, probs = 0.995) ~ "yes",
    TRUE ~ "no")
  )


### Plot 99.5% quantile
(o <- ggplot(data = out, aes(x = BIN_START, y = WEIGHTED_FST, color = OUTLIER)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = chrwide_fst, color = "blue") +
    labs(x = "Window", y = "Fst") +
    scale_color_manual(values = c("black", "red")) +
    guides(color = "none") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


### Retrieve positions of windows containing Fst peaks
max(win$WEIGHTED_FST)
print(out %>% filter(OUTLIER == "yes"), n = nrow(out))


### Genomic regions of interest
# Region1: 20160001, 20285000 (125 kb)
# Region2: 22150001, 22245000  (95 kb)


### Extract genomic regions of interest from vcf
# vcftools --gzvcf snps_hamlets_filtered.vcf.gz \
# --chr LG12 \
# --from-bp 20160001 \
# --to-bp 20285000 \
# --recode \
# --stdout | gzip > region1.vcf.gz


### Calculate PCA of outlier region (region 1)
library(vcfR)
library(adegenet)

region1_vcf <- read.vcfR("../meg23_repo/data/genome/region1.vcf.gz")

region1_gl <- vcfR2genlight(region1_vcf)

pca_region1 <- glPca(region1_gl, nf = 2)

scores_region1 <- as.data.frame(pca_region1$scores) %>%
  rownames_to_column("Sample") %>%
  mutate(Species = str_sub(Sample, -6, -4)) %>%
  as_tibble()


### Plot PCA
(r1 <- ggplot(data = scores_region1, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "Region 1 (LG12, ~20.2 kb)") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5))
)


### Save plot as PDF file
ggsave(
  filename = "PCA_region1.pdf",
  plot = r1,
  width = 14,
  height = 12,
  units = "cm",
  device = cairo_pdf
)



### ============================================================================
### Solutions

### Exercise 2: Simulating the effect of selection and mutation on allele frequencies

# Selection for a dominant phenotype / against a recessive phenotype: e.g. wAA = 1.1, wAa = 1.1, waa = 0.9
# -> A rises to high frequency (near fixation), but is not fixed (positive selection)

# Selection against a recessive phenotype: e.g. wAA = 0.9, wAa = 0.9, waa = 1.1
# -> A is purged from the population (negative selection)

# Heterozygote advantage: e.g. wAA = 0.9, wAa = 1.1, waa = 0.9
# -> both A and a are maintained at equal frequency, may lead to heterozygote excess (balancing selection)

# The strength of selection affects time scale

# Population size affects the effectiveness of selection (as opposed to drift)

# Deleterious mutations are swiftly purged

# Beneficial mutations can rise to (near-)fixation


### Plot joint Fst, along 50 kb sliding windows
(w <- ggplot(data = win, aes(x = BIN_START, y = WEIGHTED_FST)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(x = "Window", y = "Fst") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


### Add chromosome-wide mean Fst (value reported in log files)
chrwide_fst <- 0.043



#### ============================================================================
### Bonus material
  
### Smooth (spline interpolation)
smooth <- as.data.frame(spline(out$BIN_START, out$WEIGHTED_FST))

(l <- ggplot() +
    geom_line(data = smooth, aes(x = x, y = y), color = "gray20") +
    labs(x = "Window position", y = "Mean Fst") +
    geom_hline(yintercept = meanfst, color = "blue") +
    scale_y_continuous(limits=c(0, 0.6)) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


### Identify gene of interest using BLAST (example: region 2)

# Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi
# Select Nucleotide BLAST
# Paste nucleotide sequence (LG12_22220-230.fas in meg23_repo/data/genome) into search field
# Run BLAST
