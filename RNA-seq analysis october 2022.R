setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")

TPM <- read.table("Gene abundance 03112022.txt", header = T)
colnames(TPM)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/Data/Gene abundance files")

x #Data.frame with TPM data
y #vector of length 2 with root column names to be analyzed
dif_gene_express <- function(x,y){
    library(dplyr)
    library(parameters)
    df <- x
    group1_df <- df[grepl(y[1],colnames(df))]
    group2_df <- df[grepl(y[2],colnames(df))]
    rownames(group1_df) <- df$SYMBOL
    rownames(group2_df) <- df$SYMBOL
    group1_t <- as.data.frame(t(group1_df))
    group2_t <- as.data.frame(t(group2_df))
    p <- c()
    group1 <- c()
    group2 <- c()
    se_group1 <- c()
    se_group2 <- c()
    for (i in 1:(length(group1_t))){
        g1 <- group1_t[,i]
        g2 <- group2_t[,i]
        se_group1[i] <- parameters::standard_error(g1)
        se_group2[i] <- parameters::standard_error(g2)
        p[i] <- t.test(g2,g1, paired = F, alternative = "two.sided")$p.value
        group1[i] <- mean(g1)
        group2[i] <- mean(g2)
    }
    q <- p.adjust(p, method = "fdr")
    kk <- data.frame(gene = colnames(group1_t),
                     group1_name = group1,
                     group1_se = se_group1,
                     group2_name = group2,
                     group_2_se = se_group2,
                     p.value = p,
                     q.value = q,
                     Log2FC = log2(group1+1)-log2(group2+1))
    colnames(kk)[2] <- paste("mean_TPM_",strsplit(y[1],"_")[[1]], sep = "")
    colnames(kk)[3] <- paste("se_TPM_",strsplit(y[1],"_")[[1]], sep = "")
    colnames(kk)[4] <- paste("mean_TPM_",strsplit(y[2],"_")[[1]], sep = "")
    colnames(kk)[5] <- paste("se_TPM_",strsplit(y[2],"_")[[1]], sep = "")
    kk <- kk[order(-abs(kk$Log2FC)),]
    return(kk)
}
TPM_reduced <- TPM[TPM$SYMBOL %in% grs$SYMBOL,]
colnames(TPM_reduced)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
HCC827_vs_A549 <- dif_gene_express(TPM_reduced, c("HCC827_", "A549_"))
HCC827_vs_HCC827_MET <- dif_gene_express(TPM_reduced, c("HCC827_", "HCC827-MET_"))

HCC827_vs_A549 <- HCC827_vs_A549 %>% dplyr::filter(gene %ni% bad)
HCC827_vs_HCC827_MET <- HCC827_vs_HCC827_MET %>% dplyr::filter(gene %ni% bad)

volcano_dif <- function(x,y){
    library(ggplot2)
    library(ggrepel)
    p <- -log10(x$q.value)
    x$Log10p <- p
    x <- x %>% dplyr::filter(!is.na(q.value))
    x <- x %>% dplyr::filter(gene %ni% bad)
    which_m <- apply(x, MARGIN =1, FUN = function(z){
        g1 <- as.numeric(z[2])
        g2 <- as.numeric(z[4])
        return(log2(max(g1,g2)+1))
        })
    x$max_log2_TPM <- which_m
    group1_high <- length(x[x$Log2FC>0,]$Log2FC)
    group2_high <- length(x[x$Log2FC<0,]$Log2FC)
    gg <- ggplot(data = x, aes(x = Log2FC, y = Log10p,
                               size = max_log2_TPM, fill = Log2FC))+
        geom_point(shape = 21,
                   stroke = 0.5,)+
        scale_size_continuous(name = expression(bold(log[2]("TPM+1"))),
                   limits = c(0,15))+
        scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                             mid = "grey",
                             name = expression(bold(log[2]("FC"))))+
        guides(colour = guide_colourbar(order = 1),
               size = guide_legend(order = 2))+
        xlab(expression(bold(log[2]("FC"))))+
        ylab(expression(bold(paste("-",log[10]("q-value"), sep = ""))))+
        labs(title = y)+
        geom_vline(xintercept = c(0,-1,1),
                   linetype = c("solid","dashed","dashed"),
                   colour = c("black", "black", "black"),
                   size = c(1,1,1))+
        geom_label(x = 4, y = 0, label = paste(strsplit(colnames(x)[2],"_")[[1]][3]),
                   color = "#6a00fc", label.size = 0, size = 5,
                   fill="white")+
        geom_label(x = -4, y = 0, label = paste(strsplit(colnames(x)[4],"_")[[1]][3]),
                   color = "#ffa10c", label.size = 0, size = 5,
                   fill="white")+
        theme_bw(base_size = 17)+
        scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-10,10))+
        geom_label_repel(
            aes(label=ifelse(Log2FC > 1, as.character(gene),"")),
            segment.color="#6a00fc",
            color="#6a00fc",
            nudge_x = 1,
            label.size = NA,
            size = 2.5,
            fill = alpha(c("white"),0),
            parse = F,
            max.overlaps = 100)+
        geom_label_repel(
            aes(label=ifelse(Log2FC < -1, as.character(gene),"")),
            segment.color="#ffa10c",
            color="#ffa10c",
            nudge_x = -1,
            label.size = NA,
            size = 2.5,
            fill = alpha(c("white"),0),
            parse = F,
            max.overlaps = 100)+
        th
    return(gg)
}
gg1 <- volcano_dif(HCC827_vs_A549, "HCC827 compared to A549")
gg1
ggsave(filename = "A549 vs. HCC827.png",
       plot = gg1,
       width = 8513, height = 6338, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 2",
       dpi = 1200,
       device = "png")

gg2 <- volcano_dif(HCC827_vs_HCC827_MET, "HCC827 compared to HCC827-MET")

ggsave(filename = "HCC827-MET vs. HCC827.png",
       plot = gg2,
       width = 8513, height = 6338, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 2",
       dpi = 1200,
       device = "png")


nrow(TPM) 
colnames(TPM)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
umap_RNA_seq <- function(x,y){
    set.seed(8)
    library(dplyr)
    library(parameters)
    library(ggplot2)
    library(umap)
    df <- x
    group1_df <- df[grepl(y[1],colnames(df))]
    group2_df <- df[grepl(y[2],colnames(df))]
    group3_df <- df[grepl(y[3],colnames(df))]
    rownames(group1_df) <- df$SYMBOL
    rownames(group2_df) <- df$SYMBOL
    rownames(group3_df) <- df$SYMBOL
    group1_t <- as.data.frame(t(group1_df))
    group2_t <- as.data.frame(t(group2_df))
    group3_t <- as.data.frame(t(group3_df))
    combined_df <- rbind(group1_t,group2_t, group3_t)
    umap_res <- umap(combined_df, n_neighbors = 5)
    umap_df <- as.data.frame(umap_res$layout)
    umap_df[,3] <- c(rep(strsplit(y[1],"_")[[1]],3),
                         rep(strsplit(y[2],"_")[[1]],3),
                         rep(strsplit(y[3],"_")[[1]],3))
    gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V3))+
        geom_point(size = 5)+
        labs(title = "RNA-seq of cell lines",
             x = "UMAP-1", y = "UMAP-2")+
        scale_color_manual(name = "Cell line",
                             values = c("firebrick","#6a00fc","#ffa10c"))+
        theme_bw()+
        th
    return(gg)
}
gg_umap_RNA <- umap_RNA_seq(TPM,c("A549_", "HCC827_", "HCC827-MET_"))
gg_umap_RNA
ggsave(filename = "figure 2 UMAP RNA-seq.png",
       plot = gg_umap_RNA,
       width = 7563, height = 5200, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables",
       dpi = 1200,
       device = "png")

