---
title: "CGS 4144 data exploration"
output: pdf_document
date: "2025-09-22"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
#installing necessary packages / libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require(readr, quietly = TRUE))  install.packages("readr")
library(readr)
if (!require(dplyr, quietly = TRUE))  install.packages("dplyr")
library(dplyr)
if (!require(ggplot2, quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
if (!require(Rtsne, quietly = TRUE))  install.packages("Rtsne")
library(Rtsne)
if (!require(uwot, quietly = TRUE))   install.packages("uwot")
library(uwot)

#reading raw TSV file
tsv_path <- "C:/Users/anae2/Preeclampsia_Vs_Normal/refine_bio_data/GSE60438/GSE60438.tsv"
raw_data <- read_tsv(tsv_path)
gene_ids <- raw_data[[1]]
sample_names <- colnames(raw_data)[-1]

expr_mat <- raw_data %>%
  select(-1) %>%                          # drop the Gene column
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix()

rownames(expr_mat) <- gene_ids
colnames(expr_mat) <- sample_names

#plotting per gene variation in expression
#compute per-gene stats across samples
gene_stats <- data.frame(
  gene = rownames(expr_mat),
  mean = rowMeans(expr_mat, na.rm = TRUE),
  var  = apply(expr_mat, 1, var, na.rm = TRUE)
)
gene_stats$sd <- sqrt(gene_stats$var)
gene_stats$cv <- gene_stats$sd / (abs(gene_stats$mean) + 1e-8)

#plotting mean vs variance
ggplot(gene_stats, aes(x = mean, y = var)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Mean vs Variance of Gene Expression", x = "Mean Expression (log10)", y = "Variance (log10)")

  #plotting distribution of CV
ggplot(gene_stats, aes(x = cv)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Coefficient of Variation (CV)", x = "CV", y = "Frequency")


# plotting Distribution of variance 
p_hist <- ggplot(gene_stats, aes(x = var)) +
  geom_histogram(bins = 50, color = "red") +
  labs(title = "Distribution of per-gene variance", x = "Variance", y = "Genes") +
  theme_minimal()

#PCA plot
pca <- prcomp(expr_mat, center = FALSE, scale. = FALSE)
pca_df <- data.frame(sample = rownames(expr_mat), PC1 = pca$x[,1], PC2 = pca$x[,2])
library(ggplot2)
ggplot(pca_df, aes(PC1, PC2)) + geom_point() + theme_minimal()

#t-SNE plot
set.seed(42)
tsne_out <- Rtsne(t(expr_mat), perplexity = 30, verbose = TRUE, max_iter = 500)
tsne_df <- data.frame(sample = colnames(expr_mat), tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2])
ggplot(tsne_df, aes(tSNE1, tSNE2)) + geom_point() + theme_minimal()

#UMAP plot
set.seed(42)
library(uwot)
um <- umap(expr_mat, n_neighbors = min(15, max(2, nrow(expr_mat)-1)), metric = "cosine")
umap_df <- data.frame(sample = rownames(expr_mat), UMAP1 = um[,1], UMAP2 = um[,2])
ggplot(umap_df, aes(UMAP1, UMAP2)) + geom_point() + theme_minimal()
```

