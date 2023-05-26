### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 02. Hardy-Weinberg equilibrium                                           ###
### ======================================================================== ###


### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory

# Solutions and hints at bottom of page



### ============================================================================
### Exercise 1: manual calculation (see 02_within.pdf)

### What are the allele frequencies?

### What are the expected genotype frequencies?

### Is the population in Hardy-Weinberg equilibrium?



### ============================================================================
### Exercise 2: Genepop format

### Install and load required R packages
install.packages("adegenet")
install.packages("pegas")
library(adegenet)
library(pegas)


### Read in data from Genepop format file
help(read.genepop)
yellowblue <- read.genepop("../meg23_repo/data/msats/yellowblue.gen", ncode = 1)


### Access data in the new genind object
yellowblue
yellowblue@tab
yellowblue@loc.n.all
yellowblue@all.names


### Test for HWE
hw.test(yellowblue)



### ============================================================================
### Exercise 3: microsatellite data of Caribbean reef fish populations

### Read in data from Genepop format file
barbados <- read.genepop("../meg23_repo/data/msats/puella_barbados.gen", ncode = 3)
barbados


### What is the most / least diverse locus in terms of number of alleles?


### Locus summary using poppr (no. alleles, heterozygosity, evenness)
install.packages("poppr")
install.packages("genepop")
library(poppr)
library(genepop)

locus_table(barbados)


### Is this population in HWE? What does this tell us?


### Alternative test for HWE, with result over all alleles
test_HW("../meg23_repo/data/msats/puella_barbados.txt", outputFile = "HW_barbados.txt")


### Compare these results to the Cayo de Media Luna population
medialuna <- read.genepop("../meg23_repo/data/msats/puella_medialuna.gen", ncode = 3)
test_HW("../meg23_repo/data/msats/puella_medialuna.txt", outputFile = "HW_medialuna.txt")


### Load and test all Caribbean populations. What patterns can you identify regarding loci and populations?
# Use data file: ../meg23_repo/data/msats/puella_caribbean.gen



### ============================================================================
### Exercise 1 solution

### Genotype-phenotype relationship: yy = yellow, yb = green, bb = blue phenotype


### Allele frequencies
n <- 10                        # number of individuals
y <- ((6 * 2) + 1) / (n * 2)   # yellow allele frequency
b <- ((3 * 2) + 1) / (n * 2)   # blue allele frequency


### Expected genotype frequencies
yy_e <- y ^ 2
yb_e <- 2 * y * b
bb_e <- b ^ 2


### Enter data into data frame for test
dat <- data.frame(row.names = c("yy", "yb", "bb"),
                  "observed" = c(6, 1, 3),
                  "expected" = c(yy_e, yb_e, bb_e)
                  )


### Perform Pearson's chi-squared test of goodness of fit
help(chisq.test)

test <- chisq.test(dat$observed, p = dat$expected)
test
pchisq(test$statistic, df = 1, lower.tail = FALSE)   # recalculate p with 1 degree of freedom

# Note: The Chi-square test is not actually recommended for such low sample sizes, 
# we use it here because it is appropriate for most real world data


### Heterozygosity
Ho <- 1 / n
He <- 2 * y * b


### Fixation index
Fis = (He - Ho) / He



### ============================================================================
### Exercise 3 solution

## What is the most / least diverse locus in terms of number of alleles?
barbados@loc.n.all   # Number of loci and alleles
barbados@all.names   # Names of alleles for each locus


### Is this population in HWE?
hw.test(barbados)


### Load and test all Caribbean populations
caribbean <- read.genepop("../meg23_repo/data/msats/puella_caribbean.gen", ncode = 3)
test_HW("../meg23_repo/data/msats/puella_caribbean.txt", outputFile = "HW_caribbean.txt")
hw.test(caribbean)
