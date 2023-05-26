### =========================================================================###
### Exercises in Marine Ecological Genetics 2023                             ###
### 01. Introduction: R basics                                               ###
### =========================================================================###


### Simple math
3 + 5 / 2
(3 + 5) / 2


### Assigning values to variables
a <- 2   # assign value "2" to variable "a"
a        # recall value
a + 5

b <- a + 5

c <- "hello"
a + c


### Example data frame (i.e. table)
iris


### Basic R syntax: function(object, parameters)
head(iris)
help(head)   # Getting help with a function
head(iris, n = 5)


### Access specific data in a data frame (here, a column)
iris$Sepal.Width


# Optional exercise: calculate maximum, mean and standard deviation of sepal width
# Hint: search for solutions online (e.g. stackoverflow.com/questions) or
# in other/R_cheatsheet_Short.pdf


### Simple plotting example 1: histogram
hist(iris$Sepal.Width)


### Simple plotting example 2: scatter plot
plot(x = iris$Sepal.Length,
     y = iris$Petal.Length)

reg <- lm(formula = Petal.Length ~ Sepal.Length,
          data = iris)   # fit linear model to data

abline(reg, col = "red") # add linear model (regression line) to plot



### ============================================================================
### Links and resources

# R cheatsheet: https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# R for Beginners: YaRrr! The Pirate's Guide to R, https://bookdown.org/ndphillips/YaRrr/
# Advanced R: R for Data Science, https://r4ds.had.co.nz/index.html
