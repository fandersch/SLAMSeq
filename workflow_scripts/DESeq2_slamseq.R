library(PerformanceAnalytics)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(gridExtra)


MAPlot <- function(export.deseq2, case, control, cutoff) {
  
  ###  Export MA-plots of differential gene expression calling
  
  #  Extract EntrezIDs for (top 20) deregulated genes for plotting
  
  dereg <- export.deseq2[which(export.deseq2$padj <= cutoff), ] # all dereg. genes
  dereg<-dereg[!is.na(dereg$gene_name),]
  down <- nrow(dereg[dereg$log2FC_deseq2 <= 0, ]) # downregulated genes
  up <- nrow(dereg[dereg$log2FC_deseq2 >= 0, ]) # upregulated genes
  top.dereg <- dereg[order(abs(dereg$log2FC_deseq2), decreasing = T)[1:20], 1]
  top.dereg <- top.dereg[!is.na(top.dereg)]
  
  #  Create generic data frame for plotting by ggplot2.
  #  d - gray-scale density coloring for scatter plots.
  
  x <- log10(export.deseq2$avg.RPM.ctrl) # x-axis: baseline mRNA expression
  y <- export.deseq2$log2FC_deseq2 # y-axis: log2 fold-change treated/ctrl
  gene_name <- export.deseq2$gene_name
  d <- densCols(x, y, nbin = 100,
                colramp = colorRampPalette((brewer.pal(9,"Greys")[-c(1:5)])))
  df <- data.frame(x, y, d, gene_name)
  
  ### Generate unlabeled plot with accurate marginal density for use in figures.
  
  # Generate basic MA-like plot with density coloring.
  
  p <- ggplot(df, aes(x = x, y = y, label = gene_name)) +
    theme_classic() +
    scale_color_identity() +
    labs(x = "average expression total mRNA (RPM)", y = "log2FC")
  
  p.basic <- p +
    geom_point(aes(x, y, col = d), size = 1.3, shape =16) +
    geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
    geom_hline(yintercept = 0, size = 0.8) +
    geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)
  
  #  Generate marginal density plot of fold-changes.
  
  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha = .5, fill = "gray40", size = 0.8) +
    labs(x = "", y = "") +
    coord_flip() +
    geom_vline(xintercept = 0) +
    theme_classic()
  
  grid.arrange(p.basic, p.right,
               ncol = 2, widths = c(4, 1)) # assemble plots
  
  
  ###  Generate summary of MA-plots with additional information & highlights.
  #  NB: Formatting of axes may shift marginal density plot from scatter-plot.
  #  Only use exported plots for visual inspection.
  
  #  Generate basic MA-like plot with density coloring.
  
  p <- ggplot(df, aes(x = x, y = y, label = gene_name)) +
    theme_classic() +
    scale_color_identity() +
    labs(x = "average expression total mRNA (RPM)", y = "log2FC") +
    ggtitle(paste0(case," vs ",control),
            paste("n =", nrow(df),
                  "// n(up) =", up,
                  "// n(down) =", down)) +
    theme(axis.line = element_line(size = 0.5),
          axis.text	 = element_text(size = 12),
          axis.ticks	= element_line(size = 0.5))
  
  p.basic <- p +
    geom_point(aes(x, y, col = d), size = 1.3, shape =16) +
    geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
    geom_hline(yintercept = 0, size = 0.8) +
    geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)
  
  #  Generate marginal density plot of fold-changes.
  
  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha=.5, fill="gray40", size = 0.8) +
    coord_flip() +
    geom_vline(xintercept = 0) +
    theme_classic() +
    ggtitle("","") +
    theme(legend.position = "none",
          text = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          axis.line.x = element_line(color = "white"))
  
  #  Generate MA-plot with highlights and labeling of significantly dereg. genes.
  
  if(nrow(dereg) > 0){
    
    p.highlight <- p +
      geom_point(data = df[!(df$gene_name %in% dereg$gene_name),],
                 aes(x, y, col = "gray60"), size = 1.3, shape =16) +
      geom_point(data = df[df$gene_name %in% dereg$gene_name,],
                 aes(x, y, col = "red1"), size = 1.3, shape =16) +
      geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
      geom_hline(yintercept = 0, size = 0.8) +
      geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)
    
    p.highlight.2 <- p.highlight +
      geom_label_repel(data = df[df$gene_name %in% top.dereg,], max.overlaps = 20, size=2)
    
    #  Export plot.
    
    grid.arrange(p.basic, p.right,
                 ncol = 2, widths = c(4,1)) # assemble plots
    grid.arrange(p.highlight, p.right,
                 ncol = 2, widths = c(4,1)) # add highlights
    grid.arrange(p.highlight.2,p.right,
                 ncol = 2, widths = c(4,1)) # add labels
  }
}

args = commandArgs(trailingOnly=TRUE)

sample_info <- read_tsv(args[1])
output_dir <- paste0(args[2], "/")
spikein<- args[3]

sample_ids_filtered <- c()
if(!is.null(args[4])){
  sample_ids_filtered <- str_split(args[4], pattern = ",") %>% unlist %>% as.numeric
}

dir.create(output_dir)

#get sampleIDs
species <- sample_info$species %>% unique
sample_groups <- sample_info$group %>% unique

#filter out bad samples
sample_info <- sample_info %>%
  filter(!sample_id %in% sample_ids_filtered)

#input
if(spikein == "spikein"){
  cts_total <- read.csv("slamdunk/counts_collapsed/counts_scaled.txt",sep="\t")
  cts_tc <- read.csv("slamdunk/counts_collapsed/tcCounts_scaled.txt",sep="\t")
  
  #RPM
  cts_total_stripped <- cts_total[,4:ncol(cts_total)]
  ncols <- cts_total_stripped %>% ncol
  scaling_factors <- (cts_total_stripped %>% colSums())/1000000 %>% as.numeric()
  cts_total_stripped_scaled<-cts_total_stripped
  for(i in 1:ncols){
    cts_total_stripped_scaled[,i]<-(cts_total_stripped[,i]/scaling_factors[i]) %>% round(digits = 3)
  }
  RPM<- cts_total[,1:3] %>% cbind(cts_total_stripped_scaled)
  
  #RPMu
  cts_total_stripped <- cts_total[,4:ncol(cts_total)]
  cts_tc_stripped <- cts_tc[,4:ncol(cts_tc)]
  
  ncols <- cts_total_stripped %>% ncol
  cts_tc_stripped_scaled<-cts_tc_stripped
  for(i in 1:ncols){
    nonTcReadCount <- cts_total_stripped[,i] - cts_tc_stripped[,i]
    cts_tc_stripped_scaled[,i] <- ((cts_tc_stripped[,i] / sum(nonTcReadCount)) * 10^6) %>% round(digits = 3)
  }
  RPMu <- cts_tc[,1:3] %>% cbind(cts_tc_stripped_scaled)
  
}else{
  cts_total <- read.csv("slamdunk/counts_collapsed/counts.txt",sep="\t")
  cts_tc <- read.csv("slamdunk/counts_collapsed/tcCounts.txt",sep="\t")
  RPM <- read.csv("slamdunk/counts_collapsed/RPMs.txt",sep="\t")
  RPMu <- read.csv("slamdunk/counts_collapsed/RPMus.txt",sep="\t")
}

gene_symbols <- cts_total %>% dplyr::select(entrez, symbol)

cts_total<-as.matrix(cts_total[,4:ncol(cts_total)])
rownames(cts_total) <- gene_symbols$entrez
class(cts_total) <- "numeric"
colnames(cts_total) <- str_split(colnames(cts_total), "_") %>% lapply(`[[`, 1)
colnames(cts_total) <- colnames(cts_total) %>% str_remove(pattern="^X")


cts_tc<-as.matrix(cts_tc[,4:ncol(cts_tc)])
rownames(cts_tc) <- gene_symbols$entrez
class(cts_tc) <- "numeric"
colnames(cts_tc) <- str_split(colnames(cts_tc), "_") %>% lapply(`[[`, 1)
colnames(cts_tc) <- colnames(cts_total) %>% str_remove(pattern="^X")

colnames(RPM) <- str_split(colnames(RPM), "_") %>% lapply(`[[`, 1)
colnames(RPM) <- colnames(RPM) %>% str_remove(pattern="^X")
colnames(RPMu) <- str_split(colnames(RPMu), "_") %>% lapply(`[[`, 1)
colnames(RPMu) <- colnames(RPMu) %>% str_remove(pattern="^X")

###################
# PCA all samples #
###################
output_dir_group<- paste0(output_dir, "/PCA_all/")
#specify output directory
dir.create(output_dir_group, recursive = T)

#make ordered metadata matrix
coldata <- sample_info %>%
  dplyr::select(sample_id, condition) %>%
  distinct %>%
  arrange(sample_id) %>%
  column_to_rownames(var = "sample_id")

coldata$condition <- factor(coldata$condition)
cols <- rownames(coldata)

#make dds from cts_total, coldata
dds <- DESeqDataSetFromMatrix(countData = round(cts_tc[,as.character(c(cols))]),
                              colData = coldata,
                              design = ~ condition)

dds.total <- DESeqDataSetFromMatrix(countData = round(cts_total[,as.character(c(cols))]),
                                    colData = coldata,
                                    design = ~ condition)

dds.total <- dds.total[ rowSums(counts(dds)) > 10, ] #  remove uninformative rows

dds <- dds[ rowSums(counts(dds)) > 10, ] #  remove uninformative rows
#  Run deseq main command.
if(spikein == "spikein"){
  factors <- rep(1, times=length(cols))
  names(factors) <- cols
  
  #total counts
  sizeFactors(dds.total) <- factors #Don't let Deseq2 rescale the libraries, since they are already pre-scaled from the spike-in normalization
  dds.total <- DESeq(dds.total)
  #tc counts
  sizeFactors(dds) <- factors #Don't let Deseq2 rescale the libraries, since they are already pre-scaled from the spike-in normalization
}else{
  dds.total <- DESeq(dds.total)
  sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts
}

# Input is a matrix of log transformed values
vst <- varianceStabilizingTransformation(dds, blind=T)
vst_mat <- assay(vst)
vst_mat <- vst_mat %>% cbind(rowVars(vst_mat))
vst_mat_ordered <- vst_mat[order(vst_mat[,ncol(vst_mat)], decreasing = T),]

pca <- prcomp(t(vst_mat_ordered[1:500,1:ncol(vst_mat)-1]))
sdev<-pca$sdev
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(coldata, pca$x) %>% 
  cbind(sdev) %>% 
  rownames_to_column(var="sample_id")

write_tsv(df, file=paste0(output_dir_group, "matrix_PCA.txt"), col_names = T)

PoV = df$sdev ^ 2 / sum(pca$sdev ^ 2) * 100

#plot pca with labels
ggplot(df, aes(x= PC1, y= PC2, colour=condition, label=paste0(sample_id,"_",condition))) +
  xlab(paste("PC1 (", round(PoV[1],digits=2), " % variance)",sep="")) +
  ylab(paste("PC2 (", round(PoV[2],digits=2), " % variance)",sep="")) +
  geom_point() + 
  geom_text_repel(aes(label=paste0(sample_id,"_",condition)), size=2, max.overlaps=100, force = 2) + 
  theme(aspect.ratio=3/4)

ggsave(paste0(output_dir_group, "PCA_sample_id_PC1_PC2.png"), height=10,width=10)


#get sample groups
df_LFC_all_group <- NULL
df_LFC_all <- NULL
# dev.off()
for(sample_group in sample_groups){
  print(sample_group)
  gc()
  
  #get sample conditions per group
  sample_conditions <- sample_info %>%
    filter(group==sample_group) %>%
    dplyr::select(control, condition) %>%
    distinct
  
  #get tretment conditions
  sample_group_treatments <- sample_conditions %>%
    filter(control==0) %>%
    .$condition
  
  #get control
  sample_group_control <- sample_conditions %>%
    filter(control==1) %>%
    .$condition
  
  #make ordered metadata matrix
  coldata <- sample_info %>%
    filter(group==sample_group) %>% 
    arrange(sample_id) %>%
    column_to_rownames(var = "sample_id") %>%
    dplyr::select(condition)
  
  coldata$condition <- factor(coldata$condition)
  cols <- rownames(coldata)
  
  #make dds from cts_total, coldata
  dds <- DESeqDataSetFromMatrix(countData = round(cts_tc[,as.character(c(cols))]),
                                colData = coldata,
                                design = ~ condition)
  
  dds.total <- DESeqDataSetFromMatrix(countData = round(cts_total[,as.character(c(cols))]),
                                      colData = coldata,
                                      design = ~ condition)
  
  dds.total <- dds.total[ rowSums(counts(dds)) > 10, ] #  remove uninformative rows
  
  dds <- dds[ rowSums(counts(dds)) > 10, ] #  remove uninformative rows
  #  Run deseq main command.
  if(spikein == "spikein"){
    factors <- rep(1, times=length(cols))
    names(factors) <- cols
    
    #total counts
    sizeFactors(dds.total) <- factors #Don't let Deseq2 rescale the libraries, since they are already pre-scaled from the spike-in normalization
    dds.total <- DESeq(dds.total)
    #tc counts
    sizeFactors(dds) <- factors #Don't let Deseq2 rescale the libraries, since they are already pre-scaled from the spike-in normalization
  }else{
    dds.total <- DESeq(dds.total)
    sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts
  }
  
  
  if(!isEmpty(sample_group_control)){
    refs=sample_group_control
  }else{
    refs=levels(coldata$condition)
  }
  
  for(ref in refs){
    #specify output directory
    output_dir_group<- paste0(output_dir, sample_group, "/vs_", ref,"/")
    output_dir_group_qc<-paste0(output_dir, sample_group, "/vs_", ref, "/qc/")
    dir.create(output_dir_group, recursive = T)
    dir.create(output_dir_group_qc, recursive = T)
    
    dds$condition <- relevel(dds$condition, ref = ref)
    
    DEA <- function(dds) {
      tryCatch(
        {
          dds <- DESeq(dds) 
          return(dds)
        },
        error=function(cond) {
          message("Here's the original error message:")
          message(cond)
          dds <- DESeq(dds, fitType="mean")
          write("Estimate dispersion: mean - use the mean of gene-wise dispersion estimates.", file = paste0(output_dir_group, "/estimate_dispersion_parameter.txt"))
          return(dds)
        },
        warning=function(cond) {
          message("Here's the original warning message:")
          message(cond)
          dds <- DESeq(dds, fitType="mean") 
          write("Estimate dispersion: mean - use the mean of gene-wise dispersion estimates.", file = paste0(output_dir_group, "/estimate_dispersion_parameter.txt"))
          return(dds)
        },
        finally={
        }
      )
    }
    
    dds <- DEA(dds)
    
    #  Export PCA plot (default deseq PCA on 500 most variable genes).
    # Input is a matrix of log transformed values
    vst <- varianceStabilizingTransformation(dds, blind=T)
    vst_mat <- assay(vst)
    vst_mat <- vst_mat %>% cbind(rowVars(vst_mat))
    vst_mat_ordered <- vst_mat[order(vst_mat[,ncol(vst_mat)], decreasing = T),]
    
    pca <- prcomp(t(vst_mat_ordered[1:500,1:ncol(vst_mat)-1]))
    sdev<-pca$sdev
    # Create data frame with metadata and PC3 and PC4 values for input to ggplot
    df <- cbind(coldata, pca$x) %>% 
      cbind(sdev) %>% 
      rownames_to_column(var="sample_id")
    
    write_tsv(df, file=paste0(output_dir_group_qc, "matrix_PCA.txt"), col_names = T)
    
    PoV = df$sdev ^ 2 / sum(pca$sdev ^ 2) * 100
    
    #plot pca with labels
    ggplot(df, aes(x= PC1, y= PC2, colour=condition, label=sample_id)) +
      xlab(paste("PC1 (", round(PoV[1],digits=2), " % variance)",sep="")) +
      ylab(paste("PC2 (", round(PoV[2],digits=2), " % variance)",sep="")) +
      geom_point() + 
      geom_label_repel(aes(label=sample_id), size=5) + 
      theme(aspect.ratio=3/4)
  
    ggsave(paste0(output_dir_group_qc, "PCA_PC1_PC2.png"))
  
    resultsNames(dds)
    norm_tc_counts <- counts(dds, normalized=T) %>%
      round(2) %>%
      as.data.frame() %>%
      rownames_to_column("entrez") %>%
      mutate(entrez = as.numeric(entrez)) %>%
      left_join(gene_symbols) %>%
      dplyr::select(entrez, symbol, everything()) %>%
      arrange(symbol)
      
    write_tsv(norm_tc_counts, file=paste0(output_dir_group, "normalized_tc_read_counts.tsv"))
    
    norm_counts <- counts(dds.total, normalized=T) %>%
      round(2) %>%
      as.data.frame() %>%
      rownames_to_column("entrez") %>%
      mutate(entrez = as.numeric(entrez)) %>%
      left_join(gene_symbols) %>%
      dplyr::select(entrez, symbol, everything()) %>%
      arrange(symbol)
    
    write_tsv(norm_counts, file=paste0(output_dir_group, "normalized_read_counts.tsv"))
    
    #sample distance heatmap
    vst <- varianceStabilizingTransformation(dds, blind=T)
    vst_mat <- assay(vst)
    
    variances <- apply(vst_mat, 1, var)
    mean <- apply(vst_mat, 1, mean)
  
    # Determine the upper quartile variance/mean cutoff value
    upper_var <- quantile(variances, 0.25)
    upper_mean <- quantile(mean, 0.25)
  
    # Filter the data choosing only genes that are variable
    df_by_var_mean <- data.frame(vst_mat) %>%
        dplyr::filter(variances > upper_var, mean > upper_mean)
    
    sampleDists <- dist(t(vst_mat), method = "euclidean")
    sampleDistMatrix <- as.matrix(sampleDists)
    
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pdf(file = paste0(output_dir_group_qc,"sample_distance.pdf"), width = 20, height=20)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
    dev.off()
    
    colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
  
    pdf(file = paste0(output_dir_group_qc,"clustering_heatmap.pdf"), width = 20, height=20)
    pheatmap(
      df_by_var_mean,
      cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
      cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
      show_rownames = FALSE, # There are too many genes to clearly show the labels
      col=colors,
      scale = "row" # Scale values in the direction of genes (rows)
    )
    dev.off()
    
    #correlation counts replicates
    png(paste0(output_dir_group_qc,"corr_replicates.png"),width = 250, height = 250, units='mm', res = 200)
    
    norm_counts_analytics<- norm_tc_counts[,c(cols)] %>%
      as.data.frame() %>%
      mutate(rowMean = rowMeans(norm_tc_counts[,c(cols)])) %>%
      filter(rowMean>0) %>%
      dplyr::select(-rowMean)
    
    chart.Correlation(norm_counts_analytics, histogram=TRUE, pch=19)
    dev.off()
    
    if(length(ref)>0){
      for(i in 2:length(resultsNames(dds))){
        
        contrast <- str_remove(resultsNames(dds)[i], "condition_")
        
        print(contrast)
        sample_group_treatment <- contrast %>% str_replace_all(paste0("_vs_", ref), "")
        #results
        # res <- results(dds, name=contrast, lfcThreshold = 0.585, alpha=0.05) 
        res <- results(dds, name=resultsNames(dds)[i]) 
      
        #volcano plot
        # create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
        # this can be achieved with nested ifelse statements
        p_thresh <- res %>%
          as.data.frame() %>%
          filter(padj < 0.1) %>%
          arrange(desc(pvalue)) %>%
          .$pvalue
        
        p_thresh<-p_thresh[1]
        keyvals <- ifelse(
          res$log2FoldChange < -0.5 & res$pvalue < p_thresh, 'blue',
          ifelse(res$log2FoldChange > 0.5 & res$pvalue < p_thresh, 'red',
                 'black'))
        keyvals[is.na(keyvals)] <- 'black'
        names(keyvals)[keyvals == 'red'] <- 'up'
        names(keyvals)[keyvals == 'black'] <- 'neutral'
        names(keyvals)[keyvals == 'blue'] <- 'down'
        
        labels <- data.frame(entrez=rownames(res)) %>%
          mutate(entrez=as.numeric(entrez)) %>%
          left_join(gene_symbols)
        
        EnhancedVolcano(res,
                        lab = labels$symbol,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        colCustom=keyvals,
                        labSize = 2.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.2,
                        max.overlaps=25,
                        arrowheads=F,
                        pCutoff = p_thresh[1],
                        FCcutoff = 0.5
        )
        
        ggsave(paste0(output_dir_group,contrast, "_volcano.pdf"), width = 12, height=12)
  
        col_treatment <- sample_info %>%
          filter(group==sample_group) %>%
          filter(condition == sample_group_treatment) %>%
          .$sample_id %>% 
          as.character
        
        col_control <- sample_info %>%
          filter(group==sample_group) %>%
          filter(condition == ref) %>%
          .$sample_id %>% 
          as.character
        
        df <- norm_tc_counts[,c("entrez", "symbol", col_treatment, col_control)] %>% as.data.frame
        df$mean_norm_tc_count_treated <- rowMeans(norm_tc_counts[,c(col_treatment)]) %>% unlist
        df$mean_norm_tc_count_control <- rowMeans(norm_tc_counts[,c(col_control)]) %>% unlist
        
        RPMu_add <- RPMu[,c("entrez", "symbol", col_treatment, col_control)] %>% as.data.frame
        RPMu_add$mean_RPMu_treated <- rowMeans(RPMu[,c(col_treatment)]) %>% unlist
        RPMu_add$mean_RPMu_control <- rowMeans(RPMu[,c(col_control)]) %>% unlist
        pseudocount <- RPMu[,c(col_treatment, col_control)] %>% as.matrix() %>% rowMins()
        pseudocount <- min(pseudocount[pseudocount!=0])
        RPMu_add$LFC_mean_RPMu <- log2((unlist(rowMeans(RPMu[,c(col_treatment)]))+pseudocount)/(unlist(rowMeans(RPMu[,c(col_control)]))+pseudocount))
        RPMu_add <- RPMu_add %>% dplyr::select("entrez", "symbol", mean_RPMu_treated, mean_RPMu_control, LFC_mean_RPMu)
  
        RPM_add <- RPM[,c("entrez", "symbol", col_treatment, col_control)] %>% as.data.frame
        RPM_add$mean_RPM_treated <- rowMeans(RPM[,c(col_treatment)]) %>% unlist
        RPM_add$mean_RPM_control <- rowMeans(RPM[,c(col_control)]) %>% unlist
        RPM_add <- RPM_add %>% dplyr::select("entrez", "symbol", mean_RPM_treated, mean_RPM_control)
    
        df <- df %>%
          full_join(
            as.data.frame(res) %>% 
              rownames_to_column(var = "entrez") %>%
              mutate(entrez=as.numeric(entrez)),
            by= "entrez") %>% 
          full_join(RPMu_add) %>%
          full_join(RPM_add) %>%
          mutate_at(vars(c("baseMean", "log2FoldChange", "lfcSE", "mean_norm_tc_count_treated", "mean_norm_tc_count_control",  
                           "mean_RPMu_treated", "mean_RPMu_control", "LFC_mean_RPMu", "mean_RPM_treated", "mean_RPM_control",
                           all_of(col_treatment), all_of(col_control))), round, 2) %>%
          dplyr::select(symbol, entrez, everything()) %>%
          arrange(pvalue)
          
        #Exporting results to CSV files
        write_tsv(df, file=paste0(output_dir_group, contrast, "_DESeq2_all.tsv"))
        
        
        #MA-plot
        export.deseq2 <- df  %>%
          dplyr::select(gene_name=symbol, log2FC_deseq2=log2FoldChange, padj, avg.RPM.ctrl=mean_RPM_control)
        
        pdf(file.path(output_dir_group, paste0(contrast, "_MAPlot.pdf")), width = 12, height=12)
        MAPlot(export.deseq2, sample_group_treatment, sample_group_control, 0.1)
        dev.off()
        
        
        #pepare LFC matrix 
        df_LFC <- df %>%
          dplyr::select(symbol, entrez, log2FoldChange, pvalue, padj, 
                        mean_RPMu_treated, mean_RPMu_control, mean_RPM_treated, mean_RPM_control
                        )
        
        names(df_LFC)[names(df_LFC) == "log2FoldChange"] <- paste0(contrast, "_LFC")
        names(df_LFC)[names(df_LFC) == "pvalue"] <- paste0(contrast, "_pvalue")
        names(df_LFC)[names(df_LFC) == "padj"] <- paste0(contrast, "_padj")
        names(df_LFC)[names(df_LFC) == "mean_RPMu_treated"] <- paste0(contrast, "_mean_RPMu_treated")
        names(df_LFC)[names(df_LFC) == "mean_RPMu_control"] <- paste0(contrast, "_mean_RPMu_control")
        names(df_LFC)[names(df_LFC) == "mean_RPM_treated"] <- paste0(contrast, "_mean_RPM_treated")
        names(df_LFC)[names(df_LFC) == "mean_RPM_control"] <- paste0(contrast, "_mean_RPM_control")
        
        if(is.null(df_LFC_all_group)){
          df_LFC_all_group <- df_LFC
        }else{
          df_LFC_all_group <- df_LFC_all_group %>%
            full_join(df_LFC) 
        }
        
        if(is.null(df_LFC_all)){
          df_LFC_all <- df_LFC
        }else{
          df_LFC_all <- df_LFC_all %>%
            full_join(df_LFC) 
        }
      }
      # write_tsv(df_LFC_all_group, file=paste0(output_dir_group,"DESeq2_compact_", sample_group,".tsv"), na = "")
      df_LFC_all_group<-NULL
    }
  }
}

write_tsv(df_LFC_all %>% arrange(symbol), file=paste0(output_dir,"DESeq2_all.tsv"), na = "")

