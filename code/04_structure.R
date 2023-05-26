### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 04. Population structure and gene flow                                   ###
### ======================================================================== ###


### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory

# Solutions and hints at bottom of page



### ============================================================================
### Homework / recap: Ne and fluctuating population size

### Simulate a population experiencing a severe bottleneck
N <- c(100, 120, 80, 110, 10)
N
Na <- mean(N)                 # arithmetic mean
Na                            # explicitly: Na <- (100 + 120 + 80 + 110 + 10) / 5


### Calculate the harmonic mean = effective population size
Ne <- 1 / mean(1 / N)   # explicitly: Ne <- 1 / ((1/100 + 1/120 + 1/80 + 1/110 + 1/10) / 5)
Ne


### Add 45 generations with population size = 100 after the bottleneck above
N <- c(100, 120, 80, 110, 10, rep(100, 45))


### Plot Ne over time using for loop
Ne <- 0   # initialize variable

for (i in 1:length(N)) {
  Ne[i] <- 1 / mean(1 / N[1:i])
}

df <- data.frame(Generation = 1:length(Ne), Ne)
df

plot(x = df$Generation, y = df$Ne)



### ============================================================================
### Exercise 1: Compare the amount of genetic structure in three datasets

### Load packages
# install.packages("adegenet")
# install.packages("genepop")
# install.packages("pegas")
library(adegenet)
library(genepop)
library(pegas)


### Read in data
caribbean <- read.genepop("../meg23_repo/data/msats/puella_caribbean.gen", ncode = 3)
temporal <- read.genepop("../meg23_repo/data/msats/puella_temporal.gen", ncode = 3)
hamlets <- read.genepop("../meg23_repo/data/msats/hamlets_caribbean.gen", ncode = 3)


### Calculate F-statistics for each locus
pegas::Fst(as.loci(caribbean))
pegas::Fst(as.loci(temporal))
pegas::Fst(as.loci(hamlets))


### Calculate rho_st using Genepop for each dataset
help(Fst)

# hints: use genepop::Fst(), to differentiate this function from the one in the pegas package
# rho_st was developed for microsatellites, which are allele-size based markers


### Test for population structure / genetic differentiation (G-test)
help(test_diff)



### ============================================================================
### Exercise 2: Calculate population-specific Fst

### Load packages
install.packages("hierfstat")
install.packages("tidyverse")
library(hierfstat)
library(tidyverse)


### Population-specific Fst (beta)
betas(hamlets)


### Convert output to tidyverse data frame (tibble) -- execute step by step to follow the pipeline
b <- betas(hamlets)$betaiovl %>%                   # extract Fsts from betas object
  bind_rows() %>%                                  # convert to tibble
  pivot_longer(cols = everything()) %>%            # transform data from wide to long format
  rename("Population" = name, "Fst" = value) %>%   # rename columns
  arrange(desc(Fst))                               # sort data by descending Fst


### Basic plotting with ggplot
ggplot(data = b, aes(x = Population, y = Fst)) +
  geom_bar(stat = "identity")


### Make plot prettier
ggplot(data = b, 
       aes(x = reorder(Population, -Fst), y = Fst)) +   # determine basic data structure (mapping)
  geom_bar(stat = "identity",                           # use barplot geom, set outline and fill color
           color = "grey20", 
           fill = "skyblue3") +                                        
  labs(x = "Population") +                              # change x-axis label
  theme_minimal(base_size = 16) +                       # use theme "minimal" (e.g. no frame, white background)
  theme(
    panel.grid.minor = element_blank(),                 # adjust theme: remove minor grid lines
    panel.grid.major.x = element_blank(),               # adjust theme: remove major grid lines intersecting x-axis
    axis.text.x = element_text(angle = 45,              # adjust theme: change angle and position of x-axis labels
                               hjust = 1, 
                               vjust = 1.25)    
)


### Test differentiation between pairs of populations (hamlet dataset)
help(test_diff)



### ============================================================================
### Exercise 3: PCoA of genetic distances

### Calculate matrix of pairwise Fst
d <- genet.dist(hamlets, method = "Nei87")   # pairwise Fst following Nei (1978)
d


### Basic plotting
p <- pcoa(as.matrix(d))


### Create tidyverse data frame (tibble) with first two axes
t <- tibble(pco1 = p$vecp[, 1], 
            pco2 = p$vecp[, 2], 
            population = colnames(as.matrix(d)))

g <-
  ggplot(data = t,
         aes(x = pco1, y = pco2)) +
  geom_point(color = "skyblue3", size = 3) +
  labs(x = "PCoA 1",
       y = "PCoA 2") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent")
  )
g


### Add population labels
install.packages("ggrepel")
library(ggrepel)

g + geom_text_repel(aes(label = population), point.padding = 10)


### How does the Fst-based PCoA compare to the Fst barplot from exercise 2?



### ============================================================================
### Exercise 1 solution

### Calculate rho_st using Genepop for each dataset
genepop::Fst("../meg23_repo/data/msats/puella_caribbean.txt", sizes = TRUE, outputFile = "Rst_caribbean.txt")
genepop::Fst("../meg23_repo/data/msats/puella_temporal.txt", sizes = TRUE, outputFile = "Rst_temporal.txt")
genepop::Fst("../meg23_repo/data/msats/hamlets_caribbean.txt", sizes = TRUE, outputFile = "Rst_hamlets.txt")


### Test for population structure / genetic differentiation (G-test)
test_diff("../meg23_repo/data/msats/puella_caribbean.txt", outputFile = "Diff_caribbean.txt")
test_diff("../meg23_repo/data/msats/puella_temporal.txt", outputFile = "Diff_temporal.txt")
test_diff("../meg23_repo/data/msats/hamlets_caribbean.txt", outputFile = "Diff_hamlets.txt")



### ============================================================================
### Exercise 2 solution

test_diff("../meg23_repo/data/msats/hamlets_caribbean.txt", 
          pairs = TRUE,
          outputFile = "Pairs_hamlets.txt", 
          iterations = 500)



### ============================================================================
### Extra code: plot Ne with ggplot

ggplot(data = df, aes(x = Generation, y = Ne)) + 
  geom_point() + geom_line() +
  theme_classic(base_size = 16)
