### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 05. Isolation by distance                                                ###
### ======================================================================== ###


### Set working directory to "local" (use Files tab in bottom right panel)
getwd()    # check working directory

# Solutions and hints at bottom of page



### ============================================================================
### Exercise 1: Correlate geographic and genetic distance

### Load packages
# install.packages("adegenet")
# install.packages("hierfstat")
# install.packages("tidyverse")
library(adegenet)
library(hierfstat)
library(tidyverse)


### Read in microsatellite data (puella_caribbean.gen in data/msats)
caribbean <- read.genepop()


### Calculate pairwise Fst (method following Nei 1987)
d <- genet.dist()


### Convert to tibble
m <- as.matrix(d)
m[upper.tri(m, diag = TRUE)] <- NA   # discard upper triangle (including diagonal) by setting it to 'NA'

fst <- as_tibble(m[-1, -15], rownames = "popB")   # convert, discard first row and last column


### Read in pairwise geographic distances in km
geo <- read_csv("../meg23_repo/data/msats/puella_caribbean_geo.csv", col_names = TRUE, col_types = "cd")   # pairwise geographical distances in km


### Convert data from wide to long format
geo_long <- pivot_longer(geo,
                         cols = 2:15,
                         names_to = "popA",
                         values_to = "km",
                         values_drop_na = TRUE)

fst_long <- pivot_longer(fst,
                         cols = 2:15,
                         names_to = "popA",
                         values_to = "fst",
                         values_drop_na = TRUE)


### Create data frame with both geographical and genetic distances 
all_dist <- full_join(geo_long, fst_long, by = c("popB", "popA"))


### Plot geographical over genetic distance
a <- ggplot()


### Linear regression analysis
reg <- lm(formula = fst ~ km, data = all_dist)
summary(reg)

a + geom_smooth(method = "lm", se = FALSE, color = "coral2")


### Save plot as PDF file
ggsave(
  filename = "IBD_puella_caribbean.pdf",
  plot = a,
  width = 14,
  height = 12,
  units = "cm",
  device = cairo_pdf
)



### ============================================================================
### Exercise 2: Investigate smaller spatial scales

### Subset analysis to Guna Yala sampling sites only
guna_yala <- c("el porvenir", "cayos holandeses", "islas puyadas", "cayos ratones", 
               "cayos ingleses", "achutupu", "isla dupak", "goedup")

guna_dist <- all_dist %>%
  filter(popA %in% guna_yala & popB %in% guna_yala)

(g <- ggplot(data = guna_dist, aes(x = km, y = fst)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, color = "coral2") +
  ylim(-0.001, 0.015) +
  labs(title = "H. puella, Guna Yala", x = "Distance (km)", y = "Pairwise Fst") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(size = 14, vjust = 2),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 2.5)
  )
)

reg_g <- lm(formula = fst ~ km, data = guna_dist)
summary(reg_g)


### Subset analysis to sampling sites from Honduras and Belize
honbel <- c()

# Hint: sites are numbered 1 to 5 in Puebla et al. 2009, Fig. 1
  
  

### ============================================================================
### Exercise 3: Test for isolation by distance

install.packages("ape")
library(ape)


### All sites
ageo <- as.matrix(geo[, -1])   # convert to matrix (discard first column)
afst <- as.matrix(fst[, -1])

ageo[upper.tri(ageo)] <- t(ageo)[upper.tri(ageo)]   # copy bottom matrix triangle to top
afst[upper.tri(afst)] <- t(afst)[upper.tri(afst)]

mantel.test(ageo, afst, permutations = 100000)   # perform Mantel test for similarity of two matrices


### Honduras and Belize
hgeo_tmp <- geo %>%
  filter(popB %in% honbel) %>%                           # drop rows not in list of sites
  select(popB, all_of(honbel), -`cayos de media luna`)   # retain only popB column, and columns in list of sites; drop last column

hfst_tmp <- fst %>%
  filter(popB %in% honbel) %>%
  select(popB, all_of(honbel), -`cayos de media luna`)

hgeo <- as.matrix(hgeo_tmp[, -1])
hfst <- as.matrix(hfst_tmp[, -1])

hgeo[upper.tri(hgeo)] <- t(hgeo)[upper.tri(hgeo)]
hfst[upper.tri(hfst)] <- t(hfst)[upper.tri(hfst)]


### Perform Mantel test
mantel.test()



### ============================================================================
### Exercise 1 solution

### Read in microsatellite data
caribbean <- read.genepop("../meg23_repo/data/msats/puella_caribbean.gen", ncode = 3)


### Calculate pairwise Fst (method following Nei 1978)
d <- genet.dist(caribbean, method = "Nei87")


### See steps above


### Plot geographical over genetic distance (minimum plot)
a <- ggplot(data = all_dist, aes(x = km, y = fst)) +
  geom_point()
a


### Plot geographical over genetic distance (improved plot)
a <- ggplot(data = all_dist, aes(x = km, y = fst)) +
  geom_point(alpha = 0.75) +
  ylim(-0.001, 0.015) +   # set range of y-axis
  labs(title = "H. puella, Caribbean", x = "Distance (km)", y = "Pairwise Fst") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(size = 14, vjust = 2),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 2.5)
  )
a



### ============================================================================
### Exercise 2 solution

### Subset analysis to sampling sites from Honduras, Belize, and Panama (without Guna Yala)
honbel <- c("carrie bow cay", "utila", "guanaja", "cayos becerros", 
         "cayos de media luna")

honbel_dist <- all_dist %>%
  filter(popA %in% honbel & popB %in% honbel)

(h <- ggplot(data = honbel_dist, aes(x = km, y = fst)) +
    geom_point(alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, color = "coral2") +
    ylim(-0.001, 0.015) +
    labs(title = "H. puella, Honduras / Belize", x = "Distance (km)", y = "Pairwise Fst") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 14, vjust = 2),
          axis.title.x = element_text(vjust = -1.5),
          axis.title.y = element_text(vjust = 2.5)
    )
)

reg_h <- lm(formula = fst ~ km, data = honbel_dist)
summary(reg_h)


### ============================================================================
### Exercise 3 solution

### Perform Mantel test
mantel.test(hgeo, hfst, permutations = 100000)


### Guna Yala only
# ggeo_tmp <- geo %>%
#   filter(popB %in% guna_yala) %>%
#   select(popB, all_of(guna_yala), -`goedup`)
# 
# gfst_tmp <- fst %>%
#   filter(popB %in% guna_yala) %>%
#   select(popB, all_of(guna_yala), -`goedup`)
# 
# ggeo <- as.matrix(ggeo_tmp[, -1])
# gfst <- as.matrix(gfst_tmp[, -1])
# 
# ggeo[upper.tri(ggeo)] <- t(ggeo)[upper.tri(ggeo)]
# gfst[upper.tri(gfst)] <- t(gfst)[upper.tri(gfst)]
# 
# mantel.test(ggeo, gfst, permutations = 100000)



### ============================================================================
### Bonus material: plot as panels

install.packages("gridExtra")
library(gridExtra)

p <- grid.arrange(g, 
                  h, 
                  a + geom_smooth(method = "lm", se = FALSE, color = "coral2"), 
                  nrow = 1)


### Save plot as PDF file
ggsave(
  filename = "IBD_panels.pdf",
  plot = p,
  width = 30,
  height = 10,
  units = "cm",
  device = cairo_pdf
)
