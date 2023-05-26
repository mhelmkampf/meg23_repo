### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 03. Genetic drift and effective population size                          ###
### ======================================================================== ###


### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory

# Solutions and hints at bottom of page



### ============================================================================
### Homework / recap: Test for HW across multiple loci and populations

### Load packages
# install.packages("adegenet")
# install.packages("genepop")
library(adegenet)
library(genepop)


### Read in data
caribbean <- read.genepop("../meg23_repo/data/msats/puella_caribbean.gen", ncode = 3)
test_HW("../meg23_repo/data/msats/puella_caribbean.txt", outputFile = "HW_caribbean.txt")


### What conclusions can be drawn from these results?



### ============================================================================
### Exercise 1: Test for heterozygote deficiency

help(test_HW)


### Is there evidence for self-fertilization occurring in hamlets?



### ============================================================================
### Exercise 2: Simulate genetic drift

# Go to http://evolutiongenetics.georgetown.edu/simulations/
# Explore Simulations | Drift Selection Mutation
# Explore Simulations | Drift Decline in H


### What do you observe as you change initial allele frequency and effective population size?



### ============================================================================
### Exercise 3: Estimate Ne using heterozygote excess

### Unzip NeEstimator.zip in meg23_repo/other and move new program directory to local
### There, execute Ne2-1.exe (or Ne2M on Mac) through terminal or double-click

### Provide path to input file puella_caribbean.gen when prompted
### For method choose heterozygote excess
### Provide output file name (e.g. NeHet_caribbean.txt)
### Results will be saved in local/NeEstimator


### Discuss whether these estimates are reliable (see quotes in 03_drift.pdf)



### ============================================================================
### Exercise 4: Ne and fluctuating population size

### Simulate a population experiencing a severe bottleneck
N <- c(100, 120, 80, 110, 10)
N
Na <- mean(N)                 # arithmetic mean
Na                            # explicitly: Na <- (100 + 120 + 80 + 110 + 10) / 5


### Calculate the harmonic mean = effective population size
Ne <- 


### Add 45 generations with population size = 100 after the bottle neck
help(rep)   # alternatively, search for "rep" in Help tab in bottom right panel
N <- c(100, 120, 80, 110, 10, )


### Plot Ne over time using for loop
Ne <- 0   # initialize variable

for (i in 1:length(N)) {
  Ne[i] <- 1 / mean(1 / N[1:i])
}

df <- data.frame(Generation = 1:length(Ne), Ne)
df

plot(x = df$Generation, y = df$Ne)


### Optional: fancier plotting with ggplot
install.packages("tidyverse")
library(tidyverse)

ggplot(data = df, aes(x = Generation, y = Ne)) + 
  geom_point() + geom_line() +
  theme_classic(base_size = 16)


### How is the effective populations size affected by the bottleneck?



### ============================================================================
### Exercise 1 solution

### Test for heterozygote deficiency
test_HW("../meg23_repo/data/msats/puella_caribbean.txt", which = "global deficit", outputFile = "Hdef_caribbean.txt")



## ============================================================================
### Exercise 3 hint

### Input file path, relative from local: 
### ../../meg23_repo/data/msats/puella_caribbean.gen



## ============================================================================
### Exercise 4 solution

### Calculate the harmonic mean = effective population size
Ne <- 1 / mean(1 / N)   # explicitly: Ne <- 1 / ((1/100 + 1/120 + 1/80 + 1/110 + 1/10) / 5)
Ne                            


### Add 45 generations with population size = 100 after the bottleneck above
N <- c(100, 120, 80, 110, 10, rep(100, 45))
