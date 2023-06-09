### ========================= ###
### 08 Linkage disequilibrium ###
### ========================= ###


## Preparations
library(tidyverse)
library(genepop)


## Set working directory to [...]/meg_ss22



### PCA

## Load data
# use data objects from previous script (07_snps.R), or:
# scores_lg12 <- read_csv("data/tables/pca_lg12_scores.csv")
# eig_lg12 <- as.numeric(read_lines("data/tables/pca_lg12_eig.csv"))


## Plot PC scores (minimal solution)
ggplot(data = scores_lg12, aes(x = PC1, y = PC2)) +
  geom_point()


## Add species and location information, improve plot
scores_lg12SL <- scores_lg12 %>%
  mutate(Species = str_sub(Sample, -6, -4),
         Location = str_sub(Sample, -3, -1))

l <- ggplot(data = scores_lg12SL, aes(x = PC1, y = PC2, color = Species, shape = Location)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "LG12") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5))
l


## Zooming into the plot
l + coord_fixed(xlim = c(-20, 20), ylim = c(-40, 0))


## Plot eigenvalues
eig_lg12

var_lg12 <- eig_lg12 / sum(eig_lg12)

barplot(var_lg12, main = "Proportion of variance explained", las = 2)   # quick barplot



### Test for linkage disequilibrium between microsatellite loci

## G-test (null hypothesis: loci are independent)
test_LD("data/msats/puella_caribbean.txt", outputFile = "local/puella_LD.txt")
test_LD("data/msats/hamlets_caribbean.txt", outputFile = "local/hamlets_LD.txt")


## Clean up temporary files (may not work on Windows)
system("rm fichier.in cmdline.txt")


### Manipulating VCF files (demo code for the UOL HPC cluster)

## Subset data, filter by MAC (already done)
# vcftools --gzvcf ../chapter2_phased_mac2.vcf.gz \
#   --keep meg_36.ids \
#   --chr LG12 \
#   --mac 4 --max-missing 0.33 --remove-indels \
#   --recode --stdout | gzip > snps_lg12_phased_36.vcf.gz


## Prune by physical distance
bcftools +prune snps_lg12_phased_36.vcf.gz \
  -n 1 -w 2kb \
  -Oz -o snps_lg12_phased_36_2kb.vcf.gz


## Prune by LD (r2)
bcftools +prune snps_lg12_phased_36.vcf.gz \
  -l 0.4 \
  -Oz -o snps_lg12_phased_36_ld4.vcf.gz


## Calculate r2 and D statistics
vcftools --gzvcf snps_lg12_phased_36.vcf.gz \
  --hap-r2 --ld-window-bp 10000 \
  --stdout | head -n 10000 > snps_lg12_phased_36.ld.tsv

vcftools --gzvcf snps_lg12_phased_36_2kb.vcf.gz \
  --hap-r2 --ld-window-bp 10000 \
  --stdout | head -n 10000 > snps_lg12_phased_36_2kb.ld.tsv

vcftools --gzvcf snps_lg12_phased_36_ld4.vcf.gz \
  --hap-r2 --ld-window-bp 10000 \
  --stdout | head -n 10000 > snps_lg12_phased_36_ld4.ld.tsv



### Plot r2 and D statistics

## Read in stats from vcftools output
un <- read_tsv("data/snps/snps_lg12_phased_36.ld.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Un") %>%
  select(r2, Dprime, Set)

kb <- read_tsv("data/snps/snps_lg12_phased_36_2kb.ld.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Kb") %>%
  select(r2, Dprime, Set)

ld <- read_tsv("data/snps/snps_lg12_phased_36_ld4.ld.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Ld") %>%
  select(r2, Dprime, Set)


## Plot r2
boxplot(un$r2, kb$r2, ld$r2, names = c("Un", "Kb", "Ld"), ylab = "r2")   # quick boxplot
boxplot(un$r2, kb$r2, ld$r2, names = c("Un", "Kb", "Ld"), ylab = "r2", outline = FALSE)   # without outliers


## Exercise: Plot D'


## Optional: Plot r2 using ggplot



### =============== Bonus material =============== ###
  
### Estimate Ne from LD with NeEstimator

## Software available in apps/NeEstimator
## To run, execute in terminal Ne2-1.exe on Windows, Ne2M on Mac, or Ne2L on Linux
## Provide path to input file when prompted, e.g. ../../data/msats/puella_caribbean.gen)
## For method choose linkage disequilibrium (1)
## Results will be saved in apps/NeEstimator



### ================== Solutions ================= ###

## Plot D'
boxplot(un$Dprime, kb$Dprime, ld$Dprime, names = c("Un", "Kb", "Ld"), ylab = "D'")


## Plot r2 using ggplot
all <- rbind(un, kb, ld)

ggplot(data = all, aes(x = fct_relevel(Set, c("Un", "Kb", "Ld")), y = Dprime)) +
  geom_boxplot() +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
