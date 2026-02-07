library(dplyr)
library(ggplot2)
library(cowplot)

### R script for enrichment analysis of indels (2-10 bp) in marker genes
# Author: Alex Samano, 2025

# indel enrichment 
#KDE is used to inform the rejection sampling process by comparing candidate gene lengths against the empirical distribution (marker genes)

set.seed(111)  # for reproducibility

# load data 
marker_INDEL_genes <- read.table("biol450_markergenes_INDELs.txt", header = TRUE, stringsAsFactors = FALSE) 

bio450_INDEL_genes <- read.table("biol450_10genomes_euchrom_genes_INDELS_counts.txt", header = F,col.names = c("gene","length","INDEL_count"), stringsAsFactors = FALSE)
bio450_INDEL_genes$INDEL_pres <-as.numeric(bio450_INDEL_genes$INDEL_count>0)

dspr_INDEL_genes <- read.table("dspr_euchrom_genes_INDELS_counts.txt", header = F,col.names = c("gene","length","INDEL_count"), stringsAsFactors = FALSE)
dspr_INDEL_genes$INDEL_pres <-as.numeric(dspr_INDEL_genes$INDEL_count>0)

sum(bio450_INDEL_genes$INDEL_pres)

# estimate the density of marker gene lengths
marker_lengths <- marker_INDEL_genes$length
marker_kde <- density(marker_lengths, kernel = "gaussian", n = 2048) #gaussian or normal

# function to interpolate density at any length
interpolate_density <- approxfun(marker_kde$x, marker_kde$y, rule = 2)

# get observed count of INDEL affected marker genes
observed_indel_count <- sum(marker_INDEL_genes$indel)

# run monte carlo simulation 100K times
num_simulations <- 100000

### dspr genome INDELS as null
dspr_indel_counts <- numeric(num_simulations)
marker_n <- nrow(marker_INDEL_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- dspr_INDEL_genes[sample(nrow(dspr_INDEL_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  dspr_indel_counts[i] <- sum(sampled_genes$INDEL_pres)
}
#write to file

write.table(dspr_indel_counts, file = "dspr_INDELS_100Ksims.txt",
            row.names = FALSE,
            col.names = FALSE)

### bio450 genomes as null
bio450_INDEL_counts <- numeric(num_simulations)
marker_n <- nrow(marker_INDEL_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- bio450_INDEL_genes[sample(nrow(bio450_INDEL_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  bio450_INDEL_counts[i] <- sum(sampled_genes$INDEL_pres)
}


write.table(bio450_INDEL_counts, file = "biol450_INDELS_100Ksims.txt",
            row.names = FALSE,
            col.names = FALSE)





# stats


bio450_mcsims<- scan("biol450_INDELS_100Ksims.txt")
dspr_mcsims <- scan("dspr_INDELS_100Ksims.txt")

#calculate mean and median
bio450_median<-median(bio450_mcsims) # 26
dspr_median<-median(dspr_mcsims) # 30

# calculate enrichment percent
bio450_enrichment_percent <- ((observed_indel_count - bio450_median) / bio450_median) * 100
dspr_enrichment_percent <- ((observed_indel_count - dspr_median) / dspr_median) * 100

# p value calculation
bio450_p_value <- (sum(bio450_mcsims <= observed_indel_count) + 1) / (num_simulations + 1)
dspr_p_value <- (sum(dspr_mcsims <= observed_indel_count) + 1) / (num_simulations + 1)






# stats for comparison to candidate count

#bio450 genomes

# 30 marker genes have an indel
bio450_total_indel_markers <- 30


#dspr genomes

# 31 marker genes have an indel
dspr_total_indel_markers<- 31


# calculate enrichment percent
bio450_enrichment_percent_total <- ((bio450_total_indel_markers - bio450_median) / bio450_median) * 100
dspr_enrichment_percent_total <- ((dspr_total_indel_markers - dspr_median) / dspr_median) * 100

# p value calculation
bio450_p_value <- (sum(bio450_mcsims >= bio450_total_indel_markers) + 1) / (num_simulations + 1)
dspr_p_value <- (sum(dspr_mcsims >= dspr_total_indel_markers) + 1) / (num_simulations + 1)




# Plot null distribution

bio450_hist_data <- hist(bio450_mcsims, breaks = 30, plot = FALSE)
dspr_hist_data <-hist(dspr_mcsims, breaks = 30, plot = FALSE)

bio450_df <- data.frame(
  mids = bio450_hist_data$mids,
  density = bio450_hist_data$density,
  group = "BIO450"
)

dspr_df <- data.frame(
  mids = dspr_hist_data$mids,
  density = dspr_hist_data$density,
  group = "DSPR"
)

#combined_df <- rbind(bio450_df, dspr_df)
dspr_binwidth <- diff(dspr_hist_data$breaks)[1]  # assuming equal-width bins
bio450_binwidth <- diff(bio450_hist_data$breaks)[1]  # assuming equal-width bins






dspr_indel_plot <- ggplot(dspr_df, aes(x = mids, y = density, fill = group)) +
  geom_bar(stat = "identity", position = "identity", width = dspr_binwidth, alpha = 0.75) +
  theme_bw() +
  labs(x = "Simulated Value", y = "Density", title = "DSPR Indel MC Distribution") +
  scale_fill_manual(values = c("DSPR" = "skyblue"))+
  geom_vline(xintercept = observed_indel_count, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = dspr_total_indel_markers, color = "black", linetype = "dashed", size = 1) +
  scale_x_continuous(limits=c(0,40),breaks = seq(0, 40, by = 5))+
  scale_y_continuous(limits=c(0,0.15),breaks = seq(0, 0.15, by = 0.05))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")



bio450_indel_plot <- ggplot(bio450_df, aes(x = mids, y = density, fill = group)) +
  geom_bar(stat = "identity", position = "identity", width = bio450_binwidth, alpha = 0.75) +
  theme_bw() +
  labs(x = "Simulated Value", y = "Density", title = "BIO450 Indel MC Distribution") +
  scale_fill_manual(values = c("BIO450" = "skyblue"))+
  geom_vline(xintercept = observed_indel_count, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = bio450_total_indel_markers, color = "black", linetype = "dashed", size = 1) +
  scale_x_continuous(limits=c(0,40),breaks = seq(0, 40, by = 5))+
  scale_y_continuous(limits=c(0,0.15),breaks = seq(0, 0.15, by = 0.05))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")



plot_grid(dspr_indel_plot,bio450_indel_plot,ncol=1)


# 500 x 600




