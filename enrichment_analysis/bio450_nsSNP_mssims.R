library(dplyr)
library(ggplot2)

### R script for enrichment analysis of nonsynonymous SNPs in marker genes
# Author: Alex Samano, 2025

# snSNP enrichment 
#KDE is used to inform the rejection sampling process by comparing candidate gene lengths against the empirical distribution (marker genes)

set.seed(111)  # for reproducibility

# load data 
marker_nsSNP_genes <- read.table("biol450_markergenes_nsSNP.txt", header = TRUE, stringsAsFactors = FALSE) 

bio450_nsSNP_genes <- read.table("biol450_10genomes_euchrom_genes_nsSNPs_counts.txt", header = F,col.names = c("gene","length","nsSNP_count"), stringsAsFactors = FALSE)
bio450_nsSNP_genes$nsSNP_pres <-as.numeric(bio450_nsSNP_genes$nsSNP_count>0)

dspr_nsSNP_genes <- read.table("dspr_euchrom_genes_nsSNPs_counts.txt", header = F,col.names = c("gene","length","nsSNP_count"), stringsAsFactors = FALSE)
dspr_nsSNP_genes$nsSNP_pres <-as.numeric(dspr_nsSNP_genes$nsSNP_count>0)

sum(bio450_nsSNP_genes$nsSNP_pres)

# estimate the density of marker gene lengths
marker_lengths <- marker_nsSNP_genes$length
marker_kde <- density(marker_lengths, kernel = "gaussian", n = 2048) #gaussian or normal

# function to interpolate density at any length
interpolate_density <- approxfun(marker_kde$x, marker_kde$y, rule = 2)

# get observed count of nsSNP affected marker genes
observed_snp_count <- sum(marker_nsSNP_genes$nsSNP)

# run monte carlo simulation 100K times
num_simulations <- 100000

### dspr genome nsSNPs as null
dspr_snp_counts <- numeric(num_simulations)
marker_n <- nrow(marker_nsSNP_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- dspr_nsSNP_genes[sample(nrow(dspr_nsSNP_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  dspr_snp_counts[i] <- sum(sampled_genes$nsSNP_pres)
}
#write to file

write.table(dspr_snp_counts, file = "dspr_nsSNPs_100Ksims.txt",
            row.names = FALSE,
            col.names = FALSE)

### bio450 genomes as null
bio450_nsSNP_counts <- numeric(num_simulations)
marker_n <- nrow(marker_nsSNP_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- bio450_nsSNP_genes[sample(nrow(bio450_nsSNP_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  bio450_nsSNP_counts[i] <- sum(sampled_genes$nsSNP_pres)
}


write.table(bio450_nsSNP_counts, file = "biol450_nsSNPs_100Ksims.txt",
            row.names = FALSE,
            col.names = FALSE)

# stats for comparison to candidate count



bio450_mcsims<- scan("biol450_nsSNPs_100Ksims.txt")
dspr_mcsims <- scan("dspr_nsSNPs_100Ksims.txt")


#calculate mean and median
dspr_median<-median(dspr_mcsims) 
bio450_median<-median(bio450_mcsims) 

# calculate enrichment percent
bio450_enrichment_percent <- ((observed_snp_count - bio450_median) / bio450_median) * 100
dspr_enrichment_percent <- ((observed_snp_count - dspr_median) / dspr_median) * 100

# p value calculation
bio450_p_value <- (sum(bio450_mcsims <= observed_snp_count) + 1) / (num_simulations + 1)
dspr_p_value <- (sum(dspr_mcsims <= observed_snp_count) + 1) / (num_simulations + 1)


# stats for comparison to candidate count

#observed numbers over all genomes

#bio450 genomes
# 31 marker genes have a nsSNP
bio450_total_nsSNP_markers <- 31


#dspr genomes
# 33 marker genes have a nsSNP
dspr_total_nsSNP_markers<- 33


# calculate enrichment percent
bio450_enrichment_percent_total <- ((bio450_total_nsSNP_markers - bio450_median) / bio450_median) * 100
dspr_enrichment_percent_total <- ((dspr_total_nsSNP_markers - dspr_median) / dspr_median) * 100

# p value calculation
bio450_p_value <- (sum(bio450_mcsims <= bio450_total_nsSNP_markers) + 1) / (num_simulations + 1)
dspr_p_value <- (sum(dspr_mcsims <= dspr_total_nsSNP_markers) + 1) / (num_simulations + 1)




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

combined_df <- rbind(bio450_df, dspr_df)
binwidth <- diff(bio450_hist_data$breaks)[1]  # assuming equal-width bins






dspr_nsSNP_plot<- ggplot(dspr_df, aes(x = mids, y = density, fill = group)) +
  geom_bar(stat = "identity", position = "identity", width = binwidth, alpha = 0.75) +
  theme_bw() +
  labs(x = "Simulated Value", y = "Density", title = "DSPR nsSNP MC Distribution") +
  scale_fill_manual(values = c("DSPR" = "skyblue"))+
  geom_vline(xintercept = observed_snp_count, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = dspr_total_nsSNP_markers, color = "black", linetype = "dashed", size = 1) +
  scale_x_continuous(limits=c(0,45),breaks = seq(0, 45, by = 5))+
  scale_y_continuous(limits=c(0,0.35),breaks = seq(0, 0.35, by = 0.1))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")



bio450_nsSNP_plot<- ggplot(bio450_df, aes(x = mids, y = density, fill = group)) +
  geom_bar(stat = "identity", position = "identity", width = binwidth, alpha = 0.75) +
  theme_bw() +
  labs(x = "Simulated Value", y = "Density", title = "BIO450 nsSNP MC Distribution") +
  scale_fill_manual(values = c("BIO450" = "skyblue"))+
  geom_vline(xintercept = observed_snp_count, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = bio450_total_nsSNP_markers, color = "black", linetype = "dashed", size = 1) +
  scale_x_continuous(limits=c(0,45),breaks = seq(0, 45, by = 5))+
  scale_y_continuous(limits=c(0,0.2),breaks = seq(0, 0.2, by = 0.1))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")





plot_grid(dspr_nsSNP_plot,bio450_nsSNP_plot, ncol=1)

