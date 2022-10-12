####Function for avg. log2TPM calculations####
x #data.frame with gene identifiers as first column, and samples as subsequent columns
#Only genes present in AVENIO surveillance panel are alowed
avg_log2TPM <- function(x){
  library(dplyr)
  gene_column_name <- colnames(x)[1]
  x_num <- x %>% dplyr::select(!gene_column_name)
  x_num <- as.data.frame(sapply(x_num[colnames(x_num)], FUN = function(x){
    res <- log2(x+1)
    return(res)
  }))
  x_num$avg <- rowMeans(x_num)
  genes <- as.vector(x[gene_column_name])[[1]]
  x_num$genes <- genes
  return(x_num[order(x_num$genes),])
}
options(scipen = 100)
####NSCLC data####
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
df_NSCLC <- read.table("GSE179879_LungCohort_TPMcounts.txt", header = T)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179879
NSCLC_AVENIO <- df_NSCLC[df_NSCLC$Genes %in% grs$SYMBOL,]

avg_NSCLC_AVENIO <- avg_log2TPM(NSCLC_AVENIO)


####PBMC data####
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
TPM_df <- read.table("GSE107011_TPM.txt", header = T)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
PBMC_cols <- colnames(TPM_df)[grepl("PBMC",colnames(TPM_df))]
PBMC_cols
PBMC_df <- TPM_df %>% dplyr::select(PBMC_cols)
library(EnsDb.Hsapiens.v86)

PBMC_df$GENEID <- unlist(lapply(stringr::str_split(rownames(PBMC_df), "[.]"), "[[", 1))
symbols <- ensembldb::select(EnsDb.Hsapiens.v86, 
                             keys=PBMC_df$GENEID, 
                             keytype = "GENEID", 
                             columns = c("SYMBOL","GENEID"))
merged_PBMC <- merge(symbols,PBMC_df)
PBMC_AVENIO <- merged_PBMC %>% dplyr::select(SYMBOL,PBMC_cols) %>% dplyr::filter(SYMBOL %in% grs$SYMBOL)

avg_PBMC_AVENIO <- avg_log2TPM(PBMC_AVENIO)

####Putting data together####
avg_NSCLC_AVENIO <- avg_NSCLC_AVENIO %>% dplyr::filter(genes %in% avg_PBMC_AVENIO$genes)
avg_PBMC_AVENIO <- avg_PBMC_AVENIO %>% dplyr::filter(genes %in% avg_NSCLC_AVENIO$genes)
avg_collected <- data.frame(genes = avg_NSCLC_AVENIO$genes,
                       PBMC = avg_PBMC_AVENIO$avg,
                       NSCLC = avg_NSCLC_AVENIO$avg)
avg_collected$log2FC <- avg_collected$NSCLC-avg_collected$PBMC
avg_collected
avg_collected[order(avg_collected$log2FC),]
library(ggplot2)
####Log2FC plot####

gg <- ggplot(data=avg_collected, aes(x=c(1:length(avg_collected$genes)), y=log2FC))+
  geom_point(aes(colour = log2FC), size = 3)+
  scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="#808080", 
                         high="#6a00fc", name = "Log2FC", 
                         limits=c((min(avg_collected$log2FC)-0.1),(max(avg_collected$log2FC)+0.1)))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "Black", size = 1)+
  geom_hline(yintercept = 2.5, color = "Black", linetype = "dashed", size = 1)+
  geom_hline(yintercept = -2.5, color = "Black", linetype = "dashed", size = 1)+
  labs(title=paste("NSCLC compared to PBMC RNA-seq"), 
       x = "Genes", y="Log2FC")+
  xlim(0, length(avg_collected$log2FC))+
  ylim(round(min(avg_collected$log2FC))-1.5, round(max(avg_collected$log2FC))+1.5)+
  geom_label(x = 30, y = round(max(avg_collected$log2FC))+1.5, label = paste("Upregulated in NSCLC"),
             color = "#6a00fc", label.size = 0, size = 5)+
  geom_label(x = 30, y = round(min(df_final$log2FC))-1.5, label = paste("Upregulated in PBMC"),
             color = "#ffa10c", label.size = 0, size = 5)+
  geom_label_repel(
    aes(label=ifelse(log2FC < -2.5,
                     as.character(genes),"")),
    segment.color="#ffa10c",
    color="#ffa10c",
    label.size = NA,
    parse = F,
    size = 3.5,
    max.overlaps = 100,
    fill = alpha(c("white"),0))+
  geom_label_repel(
    aes(label=ifelse(log2FC > 2.5,
                     as.character(genes),"")),
    segment.color="#6a00fc",
    color="#6a00fc",
    label.size = NA,
    nudge_y = 0.3,
    parse = F,
    size = 3.5,
    max.overlaps = 100,
    fill = alpha(c("white"),0))+
  th+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

####Volcano plot####
df_PBMC_AVENIO <- avg_PBMC_AVENIO %>% dplyr::select(!c("avg","genes"))
df_NSCLC_AVENIO <- avg_NSCLC_AVENIO %>% dplyr::select(!c("avg","genes"))
df_PBMC_AVENIO_t <- as.data.frame(t(df_PBMC_AVENIO))
colnames(df_PBMC_AVENIO_t) <- avg_PBMC_AVENIO$genes
df_temp <- as.data.frame(sapply(df_PBMC_AVENIO_t[colnames(df_PBMC_AVENIO_t)],as.numeric))
rownames(df_temp) <- rownames(df_PBMC_AVENIO_t)
df_PBMC_AVENIO_t <- df_temp

df_NSCLC_AVENIO_t <- as.data.frame(t(df_NSCLC_AVENIO))
colnames(df_NSCLC_AVENIO_t) <- avg_NSCLC_AVENIO$genes
df_temp <- as.data.frame(sapply(df_NSCLC_AVENIO_t[colnames(df_NSCLC_AVENIO_t)],as.numeric))
rownames(df_temp) <- rownames(df_NSCLC_AVENIO_t)
df_NSCLC_AVENIO_t <- df_temp


collect_df <- rbind(df_PBMC_AVENIO_t,df_NSCLC_AVENIO_t)
samples <- sapply(rownames(collect_df), FUN = function(x){
  if(grepl("PBMC",x)){
    return("PBMC")
  }
  else{
    return("NSCLC")
  }
})
samples
collect_df$samples <- samples

x #df with samples as rows, genes as columns, but last column is called samples
#indicating whether the sample is a PBMC or a NSCLC sample

dif_gene_express <- function(x){
  library(dplyr)
  library(parameters)
  df <- x
  NSCLC_cases <- df %>% dplyr::filter(samples == "NSCLC")
  PBMC_cases <- df %>% dplyr::filter(samples == "PBMC")
  p <- c()
  NSCLC <- c()
  PBMC <- c()
  se_PBMC <- c()
  se_NSCLC <- c()
  for (i in 1:(length(df)-1)){
    x <- NSCLC_cases[,i]
    y <- PBMC_cases[,i]
    se_PBMC[i] <- parameters::standard_error(y)
    se_NSCLC[i] <- parameters::standard_error(x)
    p[i] <- t.test(y,x, paired = F, alternative = "two.sided")$p.value
    NSCLC[i] <- mean(x)
    PBMC[i] <- mean(y)
  }
  q <- p.adjust(p, method = "fdr")
  kk <- data.frame(gene = colnames(df)[1:length(df)-1],
                   mean_NSCLC = NSCLC,
                   se_NSCLC = se_NSCLC,
                   mean_PBMC = PBMC,
                   se_PBMC = se_PBMC,
                   p.value = p,
                   q.value = q,
                   Log2FC = NSCLC-PBMC)
  return(kk)
}
dif_gene_express(collect_df)
result_dif <- dif_gene_express(collect_df)
result_dif
sort(result_dif$Log2FC)
result_dif
volcano_dif <- function(x,y){
  library(ggplot2)
  library(ggrepel)
  p <- -log10(x$q.value)
  x$Log10p <- p
  x <- x %>% dplyr::filter(!is.na(q.value))
  x <- x %>% dplyr::filter(gene %ni% bad)
  which_m <- apply(x, MARGIN =1, FUN = function(z){
    NSCLC <- as.numeric(z["mean_NSCLC"])
    PBMC <- as.numeric(z["mean_PBMC"])
    return(max(NSCLC,PBMC))})
  x$which <- which_m
  PBMC_high <- length(x[x$Log2FC<0,]$Log2FC)
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
    geom_label(x = -4, y = 1, label = paste("PBMC"),
               color = "#ffa10c", label.size = 0, size = 5,
               fill="white")+
    geom_label(x = 4, y = 1, label = paste("NSCLC"),
               color = "#6a00fc", label.size = 0, size = 5,
               fill="white")+
    theme_bw(base_size = 17)+
    scale_x_continuous(breaks = seq(-7, 7, by = 1), limits = c(-7,7))+
    geom_label_repel(
      aes(label=ifelse(Log2FC > 4, as.character(gene),"")),
      segment.color="#6a00fc",
      color="#6a00fc",
      label.size = NA,
      size = 4,
      fill = alpha(c("white"),0),
      parse = F,
      max.overlaps = 100)+
    th
  return(gg)
}
sort_result_dif <- result_dif %>% dplyr::filter(gene %ni% bad)
volcano_dif(sort_result_dif, "NSCLC compared to PBMC RNA-seq")
grs$SYMBOL[grs$SYMBOL %ni% bad]

sort_result_dif <- sort_result_dif %>% dplyr::filter(!is.na(q.value))
sort_result_dif_NSCLC <- sort_result_dif %>% dplyr::filter(Log2FC>2)
sort_result_dif_healthy <- sort_result_dif %>% dplyr::filter(Log2FC < -2)
