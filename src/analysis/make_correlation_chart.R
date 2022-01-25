# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table, janitor, magrittr, fixest, 
  broom, tidyverse, tidylog, corrplot
)
options("tidylog.display" = NULL)
`%notin%` = Negate(`%in%`)

source("src/data_cleaining_scripts/cleaning_master.R")

variables
M = dat %>% 
  select(`Carbon Cycle`:`Alternative ethical approaches (not Discounted Utilitarianism)`) %>% 
  as_tibble() %>% 
  mutate(across(everything(), ~ as.numeric(.x))) %>% 
  replace(is.na(.), 0) %>% 
  cor()

colnames(M) <- c("Carbon Cycle", "Fast Food", "Milk", "Dairy (non-milk)", "Vitamins", "Produce", "Canned/dried veg.", "% in Poverty", "Median Income", "% Black", "% Hispanic", "% Old Homes")
rownames(M) <- c("Calcium products", "Fast Food", "Milk", "Dairy (non-milk)", "Vitamins", "Produce", "Canned/dried veg.", "% in Poverty", "Median Income", "% Black", "% Hispanic", "% Old Homes")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


png("outputs/corr_chart.png", width = 10, height = 10, units = "in", res = 500)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(100),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black",
         # Combine with significance
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

