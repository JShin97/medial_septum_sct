.libPaths(c("u/home/j/jshin/project-xyang123/apps/R/4.3.0", .libPaths()))

library(tidyverse)

x <- 1:10
y <- rep(1:2, 5)
df <- data.frame("x" = x, "y" = y)
write.csv(df, file = "test.csv")
