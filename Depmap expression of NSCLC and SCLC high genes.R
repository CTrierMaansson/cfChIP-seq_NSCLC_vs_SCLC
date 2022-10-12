library(depmap)
library(ExperimentHub)
library(dplyr)
library(writexl)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/Depmap data")
colnames(CCLE_sample_info)
fun1 <- function(x){
  res <- paste(x,sep = "","_LUNG")
  return(res)
}
SCLC_cells <- read.table("SCLC cell lines depmap.txt", header = T)
SCLC_cells <- SCLC_cells$name
SCLC_cells <- as.list(SCLC_cells)
SCLC_cells <- lapply(SCLC_cells, FUN = fun1)

NSCLC_cells <- read.table("NSCLC cell lines depmap.txt", header = T)
NSCLC_cells <- NSCLC_cells$name
NSCLC_cells <- as.list(NSCLC_cells)
NSCLC_cells <- lapply(NSCLC_cells, FUN = fun1)



eh <- ExperimentHub()
query(eh, "depmap")
TPM <- eh[["EH7525"]]
lung_cells <- unique(TPM$cell_line)[grepl("LUNG",unique(TPM$cell_line))]
lung_TPM <- TPM %>% dplyr::filter(cell_line %in% lung_cells)
lung_TPM_AVENIO <- lung_TPM %>% filter(gene_name %in% grs$SYMBOL)

fun2 <- function(x){
  res <- lung_TPM_AVENIO %>% dplyr::filter(cell_line %in% x)
  if(length(res$gene_name)>0){
    return(res)
  }
  else{
    return(NA)
  }
}
fun3 <- function(x){
  rna_expression <- x$rna_expression
  return(rna_expression)
}
fun4 <- function(x){
  name <- x$cell_line[1]
  return(name)
}
fun5_SCLC <- function(x){
  val <- SCLC_expression_df[[x]]
  avg <- mean(val)
  return(avg)
}
fun5_NSCLC <- function(x){
  val <- NSCLC_expression_df[[x]]
  avg <- mean(val)
  return(avg)
}
SCLC_expression <- lapply(SCLC_cells,FUN = fun2)
SCLC_expression <- SCLC_expression[!is.na(SCLC_expression)]
length(SCLC_expression)
SCLC_expression_list <- lapply(SCLC_expression, FUN = fun3)
SCLC_expression_df <- as.data.frame(do.call(cbind,SCLC_expression_list))
rownames(SCLC_expression_df) <- SCLC_expression[[1]]$gene_name
colnames(SCLC_expression_df) <- unlist(lapply(SCLC_expression,FUN = fun4))
length(SCLC_expression_df)
SCLC_expression_df <- as.data.frame(t(SCLC_expression_df))
SCLC_means <- unlist(lapply(as.list(colnames(SCLC_expression_df)),FUN = fun5_SCLC))
SCLC_expression_df <- rbind(SCLC_expression_df,SCLC_means)
len_SCLC <- length(SCLC_expression_df$HECW1)
rownames(SCLC_expression_df)[len_SCLC] <- "Average"
SCLC_expression_df[len_SCLC,]


NSCLC_expression <- lapply(NSCLC_cells, FUN = fun2)
NSCLC_expression <- NSCLC_expression[!is.na(NSCLC_expression)]
length(NSCLC_expression)
NSCLC_expression_list <- lapply(NSCLC_expression, FUN = fun3)
NSCLC_expression_df <- as.data.frame(do.call(cbind,NSCLC_expression_list))
rownames(NSCLC_expression_df) <- NSCLC_expression[[1]]$gene_name
colnames(NSCLC_expression_df) <- unlist(lapply(NSCLC_expression,FUN = fun4))
length(NSCLC_expression_df)
NSCLC_expression_df <- as.data.frame(t(NSCLC_expression_df))
NSCLC_means <- unlist(lapply(as.list(colnames(NSCLC_expression_df)),FUN = fun5_NSCLC))
NSCLC_expression_df <- rbind(NSCLC_expression_df,NSCLC_means)
len_NSCLC <- length(NSCLC_expression_df$HECW1)
rownames(NSCLC_expression_df)[len_NSCLC] <- "Average"
NSCLC_expression_df[len_NSCLC,]


write_xlsx(df_SCLC, "SCLC depmap cells.xlsx")
write_xlsx(df_NSCLC, "NSCLC depmap cells.xlsx")
x #Enrichment data.frame return from e.score() for sample 1
y #Enrichment data.frame return from e.score() for sample 2
p #Name of sample 1
q #Name of sample 2
b #optional: list of genes returned by badgene(), otherwise type NULL
z #optional: cut-off of difference, otherwise NULL for top 15 differentially expressed genes
ChIPdif <- function(x,y,p,q,b = NULL,z = NULL){
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  `%ni%` <- Negate(`%in%`)
  if(length(x$genes) == length(y$genes)){
    y = y
    x = x
  }
  else{
    if (length(x$genes) > length(y$genes)){
      x <- x %>% filter(genes %in% y$genes)
    }
    if (length(y$genes) > length(x$genes)){
      y <- y %>% filter(genes %in% x$genes)
    }
  }
  y <- y[(match(x$genes, y$genes)),]
  df <- data.frame(genes = x$genes, x = x$enrichment, y = y$enrichment)
  df$log2FC <- log2(df$x/df$y)
  df <- df[order(df$log2FC),]
  if (is.null(b)){
    df <- df
  }
  else{
    df <- df %>% filter(genes %ni% b)
  } 
  df$reverse <- c(1:length(df$genes))
  return(df)
}
NSCLC_SCLC_ChIPdif <- ChIPdif(NSCLC_tot_enrichment, SCLC_tot_enrichment, "NSCLC", "SCLC", bad)


x #data.frame returned by ChIPdif
y #DepMap expression of cell lines associated to sample 1 (positive log2FC)
z #DepMap expression of cell lines associated to sample 2 (negative log2FC)
q #number of genes to be analyzed. Default is 15
cell_expression_enriched <- function(x,y,z,q = 15){
  library(dplyr)
  x <- x[x$genes %in% colnames(y),]
  s1_genes <- x %>% filter(log2FC > 0)
  len_s1_genes <- length(s1_genes$genes)
  s1_genes_filt <- s1_genes[(len_s1_genes-q+1):len_s1_genes,]
  s2_genes <- x %>% filter(log2FC < 0)
  len_s2_genes <- length(s2_genes$genes)
  s2_genes_filt <- s2_genes[1:q,]
  
  s1_expression_s1 <- y %>% dplyr::select(s1_genes_filt$genes)
  s1_average_s1 <- unlist(as.vector(s1_expression_s1[length(y[,1]),]))
  s2_expression_s1 <- z %>% dplyr::select(s1_genes_filt$genes)
  s2_average_s1 <- unlist(as.vector(s2_expression_s1[length(z[,1]),]))
  
  s1_expression_s2 <- y %>% dplyr::select(s2_genes_filt$genes)
  s1_average_s2 <- unlist(as.vector(s1_expression_s2[length(y[,1]),]))
  s2_expression_s2 <- z %>% dplyr::select(s2_genes_filt$genes)
  s2_average_s2 <- unlist(as.vector(s2_expression_s2[length(z[,1]),]))
  
  s1_high_res <- data.frame(s1 = s1_average_s1, 
                            s2 = s2_average_s1)
  s2_high_res <- data.frame(s1 = s1_average_s2,
                            s2 = s2_average_s2)
  return(list(s1_high = s1_high_res,
              s2_high = s2_high_res))

}
NSCLC_vs_SCLC_res <- cell_expression_enriched(NSCLC_SCLC_ChIPdif, NSCLC_expression_df, SCLC_expression_df, 15)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article")
write_xlsx(NSCLC_vs_SCLC_res$s1_high, "NSCLC_high.xlsx")
write_xlsx(NSCLC_vs_SCLC_res$s2_high, "SCLC_high.xlsx")




x #data.frame with gene expression of AVENIO genes in NSCLC cell lines
y #data.frame with gene expression of AVENIO genes in SCLC cell lines
dge_NSCLC_SCLC <- function(x,y){
  p <- c()
  NSCLC <- c()
  SCLC <- c()
  se_NSCLC <- c()
  se_SCLC <- c()
  for (i in 1:length(x)){
    NSCLC_exp <- x[,i]
    SCLC_exp <- y[,i]
    se_NSCLC[i] <- parameters::standard_error(NSCLC_exp)
    se_SCLC[i] <- parameters::standard_error(SCLC_exp)
    p[i] <- t.test(NSCLC_exp,SCLC_exp, paired = F, alternative = "two.sided")$p.value
    NSCLC[i] <- mean(NSCLC_exp)
    SCLC[i] <- mean(SCLC_exp)
  }
  q <- p.adjust(p, method = "fdr")
  df <- data.frame(genes = colnames(x),
                   mean_NSCLC = NSCLC,
                   se_NSCLC = se_NSCLC,
                   mean_SCLC = SCLC,
                   se_SCLC = se_SCLC,
                   p.value = p,
                   q.value = q,
                   Log2FC = NSCLC-SCLC)
  df <- df[order(df$Log2FC),]
  df <- df[df$genes %ni% bad,]
  df$number <- 1:length(df$genes)
  return(df)
}
dge_res <- dge_NSCLC_SCLC(NSCLC_expression_df,SCLC_expression_df)
length(dge_res$genes)

sort_FC <- dge_res$Log2FC[order(dge_res$genes)]
names(sort_FC) <- sort(dge_res$genes)
NSCLC_SCLC_ChIPdif_sort <-  NSCLC_SCLC_ChIPdif[order(NSCLC_SCLC_ChIPdif$genes),]
NSCLC_SCLC_ChIPdif_sort <- NSCLC_SCLC_ChIPdif_sort[NSCLC_SCLC_ChIPdif_sort$genes %in% dge_res$genes,]
sort_FC <- sort_FC[names(sort_FC) %in% NSCLC_SCLC_ChIPdif_sort$genes]
ratio_df <- data.frame(genes = NSCLC_SCLC_ChIPdif_sort$genes,
                       cell_FC =unlist(as.vector(sort_FC)),
                       cfChIP_FC = NSCLC_SCLC_ChIPdif_sort$log2FC)

write_xlsx(ratio_df, "ratio_df.xlsx")

dge_res
volcano <- function(x,y){
  library(ggplot2)
  library(ggrepel)
  p <- -log10(x$q.value)
  x$Log10p <- p
  which_m <- apply(x, MARGIN =1, FUN = function(z){
    NSCLC <- as.numeric(z["mean_NSCLC"])
    SCLC <- as.numeric(z["mean_SCLC"])
    return(max(NSCLC,SCLC))})
  x$which <- which_m
  SCLC_high <- length(x[x$Log2FC<0,]$Log2FC)
  NSCLC_high <- length(x[x$Log2FC>0,]$Log2FC)
  gg <- ggplot(data = x, aes(x = Log2FC, y = Log10p,
                       size = which, fill = Log2FC))+
    geom_point(shape = 21,
               stroke = 0,)+
    scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                         mid = "grey",
                         name = expression(log[2]("FC")),
                         limits= c(-7,7))+
    scale_size(name = expression(log[2]("TPM+1")),
               limits= c(0,9))+
    xlab(expression(log[2]("FC")))+
    ylab(expression(paste("-",log[10]("q-value"), sep = "")))+
    labs(title = y)+
    geom_vline(xintercept = 0,
               linetype = "solid",
               colour = "black",
               size = 1)+
    geom_label(x = -4, y = 1, label = paste("SCLC"),
               color = "#ffa10c", label.size = 0, size = 5,
               fill="white")+
    geom_label(x = 4, y = 1, label = paste("NSCLC"),
               color = "#6a00fc", label.size = 0, size = 5,
               fill="white")+
    theme_bw(base_size = 17)+
    scale_x_continuous(breaks = seq(-7, 7, by = 1), limits = c(-7,7))+
    geom_label_repel(
      aes(label=ifelse(Log10p > 10, as.character(genes),"")),
      segment.color="black",
      color="black",
      label.size = NA,
      size = 4,
      fill = alpha(c("white"),0),
      parse = F,
      max.overlaps = 10)+
    th
  return(gg)
}
volcano(dge_res, "NSCLC versus SCLC")
ggsave(filename = "NSCLC versus SCLC.png",
       width = 12275, height = 7800, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/Figures and tables/Supplementary/Supp 8",
       dpi = 1200,
       device = "png")

dge_res[dge_res$Log2FC>0.5,]
NSCLC_high <- NSCLC_SCLC_ChIPdif[NSCLC_SCLC_ChIPdif$genes %in% dge_res[dge_res$Log2FC>1,]$genes,]
SCLC_high <- NSCLC_SCLC_ChIPdif[NSCLC_SCLC_ChIPdif$genes %in% dge_res[dge_res$Log2FC<(-1),]$genes,]

write_xlsx(NSCLC_high, "NSCLC_FC.xlsx")
write_xlsx(SCLC_high, "SCLC_FC.xlsx")


ratio_df

cell_low <- ratio_df %>% filter(cell_FC < -1)
cell_high <- ratio_df %>% filter(cell_FC > 1)
cell_high_cfChIP_low <- cell_high %>% filter(cfChIP_FC < 0)
cell_high_cfChIP_high <- cell_high %>% filter(cfChIP_FC > 0)
cell_low_cfChIP_low <- cell_low %>% filter(cfChIP_FC < 0)
cell_low_cfChIP_high <- cell_low %>% filter(cfChIP_FC > 0)
contin_df <- data.frame(SCLC = c(length(cell_low_cfChIP_low$genes),
                                 length(cell_high_cfChIP_low$genes)),
                        NSCLC = c(length(cell_low_cfChIP_high$genes),
                                  length(cell_high_cfChIP_high$genes)))
rownames(contin_df) <- c("SCLC", "NSCLC")
contin_df
cell_low_cfChIP_high
length(ratio_df$genes)

ratio_df

x #ratio_df
ratio_plot <- function(x){
  fun1 <- function(z){
    cell_FC <- as.numeric(z["cell_FC"])
    cfChIP_FC <- as.numeric(z["cfChIP_FC"])
    if (cell_FC < (-1) | cell_FC > 1){
      if(cell_FC < (-1)){
        if(cfChIP_FC < 0){
          return("Agree")
        }
        else{
          return("Disagree")
        }
      }
      if(cell_FC > 1){
        if(cfChIP_FC > 0){
         return("Agree") 
        }
        else{
          return("Disagree")
        }
      }
    }
    else{
      return("filtered")
    }
  }
  res <- apply(x, FUN = fun1, MARGIN = 1)
  res <- unlist(res)
  x$group <- res
  cor.res <- cor.test(x$cfChIP_FC, x$cell_FC, method = "spearman")
  p_values <- cor.res$p.value
  if (p_values < 0.0001){
    p <- ", P < 0.0001"
  }
  else{
    p <- paste(", P =", p_values)
  }
  rhos <- as.numeric(cor.res$estimate)
  gg <- ggplot(data = x, aes(x = cell_FC, y = cfChIP_FC, colour = group))+
    geom_point(size = 2.5)+
    theme_bw()+
    scale_color_manual(name = "cfChIP results", 
                       values = c("#ffa10c", "#6a00fc", "#999999"))+
    labs(title = "cfChIP compared to CCLE expression",
           y = "cfChIP relative enrichment", x = "RNA log2FC",
         subtitle = paste("Spearman's rho =", round(rhos,3), 
                          as.character(p)))+
    scale_x_continuous(breaks = seq(-6,6,2), limits = c(-6,6))+
    scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.5,1.5))+
    geom_vline(xintercept=c(0,-1,1), 
               color=c("black", "black","black"),
               linetype=c("solid", "dashed", "dashed"), 
               size = c(1,0.5,0.5))+
    geom_hline(yintercept = 0,
               color = "black",
               linetype = "solid",
               size = 1)+
    geom_smooth(aes(x=cell_FC, y = cfChIP_FC), method = "lm", se = F, color = "black")+
    th
  return(gg)
}
ratio_plot(ratio_df)


x #ratio_df
ratio_contin <- function(x){
  library(ggpubr)
  library(tidyverse)
  rna_df_high <- x %>% filter(cell_FC > 1)
  ChIP_df_high <- x %>% filter(cfChIP_FC > 0) %>% filter(genes %in% rna_df_high$genes)
  rna_df_low <- x %>% filter(cell_FC < -1)
  ChIP_df_low <- x %>% filter(cfChIP_FC < 0) %>% filter(genes %in% rna_df_low$genes)
  ChIP_negative_high <- rna_df_high %>% filter(genes %ni% ChIP_df_high$genes)
  ChIP_negative_low <- rna_df_low %>% filter(genes %ni% ChIP_df_low$genes)
  df1 <- data.frame(discovery = factor(c("Agree", "Disagree",
                                         "Agree", "Disagree")),
                    activity = c(rep("NSCLC upregulated",2),
                                 rep("SCLC upregulated",2)),
                    numbergenes = c(length(ChIP_df_high$genes),
                                    length(ChIP_negative_high$genes),
                                    length(ChIP_df_low$genes),
                                    length(ChIP_negative_low$genes)))
  percent1 <- (df1$numbergenes[1])/(df1$numbergenes[1]+df1$numbergenes[2])*100
  percent2 <- (df1$numbergenes[3])/(df1$numbergenes[3]+df1$numbergenes[4])*100
  height1 <- length(ChIP_df_high$genes)+length(ChIP_negative_high$genes)
  height2 <- length(ChIP_df_low$genes)+length(ChIP_negative_low$genes)
  scales <- c(height1/5, height2/5)
  height1 <- height1 + min(scales)
  height2 <- height2 + min(scales)
  df1$discovery2 <- relevel(df1$discovery, 'Disagree')
  gg <- ggplot(data=df1, aes(x = fct_inorder(activity), y = numbergenes, fill = discovery2))+
    geom_bar(stat = "identity")+
    theme_bw()+
    labs(x = "", y="Number of genes")+
    scale_fill_manual("cfChIP results", values = c("#6a00fc","#ffa10c"))+
    geom_text(x = 1, y = height1, 
              label = paste(round(percent1,2),"%"),
              color = "black", size = 5)+
    geom_text(x = 2, y = height2, 
              label = paste(round(percent2,2),"%"),
              color = "black", size = 5)+
    ylim(0,(max(height1,height2)+2))+
    th+
    theme(axis.text.x = element_text(size = 14, face = "bold", colour = "black"))
  return(gg)
}
ratio_contin(ratio_df)
