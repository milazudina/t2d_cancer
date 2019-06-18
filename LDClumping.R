# colMin <- function(X, na.rm = TRUE) apply(X, MARGIN = 2, FUN = min, na.rm = na.rm)
# colMax <- function(X, na.rm = TRUE) apply(X, MARGIN = 2, FUN = max, na.rm = na.rm)
# 
# z_min_each_phen <- colMin(z_table[, 7:ncol(z_table)])
# min_of_min <- which.min(z_min_each_phen)  
# 
# z_max_each_phen <- colMax(z_table[, 7:ncol(z_table)])
# max_of_max <- which.max(z_max_each_phen)  
# 
# d <- density(z_table$z_t2d)
# plot(d)
library(dplyr)
z_table_bcac <- z_table[order(LD, z_breast_cancer), ]

res1 <- z_table %>% 
  group_by(LD) %>%
  slice(which.max(z_breast_cancer))
