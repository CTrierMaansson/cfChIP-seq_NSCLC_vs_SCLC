
library(ggplot2)
th <- theme(
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 10, vjust = 1, hjust = 0.37),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    legend.text = element_text(size = 12), 
    title = element_text(size = 14, face = "bold"))

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
bad <- badgene("Complete table of all targets.txt",25)
grs <- gr("AVENIO_genes.txt")
setwd("D:/Cell ChIP/Deduped")
A549_ChIP <- gene_count("A549_ChIP.bam", grs)
HCC827_ChIP <- gene_count("HCC827ChIP.bam", grs)
Clone3_ChIP <- gene_count("HCC827_ER_klon3_ChIP.bam", grs)

A549_ChIP_R1 <- gene_count("A549_R1_ChIP.bam", grs)
A549_ChIP_R2 <- gene_count("A549_R2_ChIP.bam", grs)
A549_ChIP_R3 <- gene_count("A549_R3_ChIP.bam", grs)
HCC827_ChIP_R1 <- gene_count("HCC827_R1_ChIP.bam", grs)
HCC827_ChIP_R2 <- gene_count("HCC827_R2_ChIP.bam", grs)
HCC827_ChIP_R3 <- gene_count("HCC827_R3_ChIP.bam", grs)
HCC827_MET_ChIP_R1 <- gene_count("HCC827-MET_R1_ChIP.bam", grs)
HCC827_MET_ChIP_R2 <- gene_count("HCC827-MET_R2_ChIP.bam", grs)
HCC827_MET_ChIP_R3 <- gene_count("HCC827-MET_R3_ChIP.bam", grs)
HCC827_MET_ChIP_R3

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
metrics("A549_R1_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("A549_R2_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("A549_R3_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827_R1_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827_R2_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827_R3_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827-MET_R1_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827-MET_R2_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")
metrics("HCC827-MET_R3_ChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "D:/Cell ChIP/Deduped")

setwd("D:/Cell input/Deduped")
A549_input <- gene_count("A549_input.bam", grs)
Clone3_input <- gene_count("HCC827_ER_klon3_input.bam", grs)
HCC827_input <- gene_count("HCC827input.bam", grs)
gc()
HCC827_bam <- bam("HCC827ChIP.bam", "HCC827ChIP.bam.bai")
A549_bam <- bam("A549_ChIP.bam", "A549_ChIP.bam.bai")
Clone3_bam <- bam("HCC827_ER_klon3_ChIP.bam", "HCC827_ER_klon3_ChIP.bam.bai")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/BAM files and manifests")
A549_ChIP_R1_bam <- bam("A549_R1_ChIP.bam", "A549_R1_ChIP.bam.bai")
A549_ChIP_R2_bam <- bam("A549_R2_ChIP.bam", "A549_R2_ChIP.bam.bai")
A549_ChIP_R3_bam <- bam("A549_R3_ChIP.bam", "A549_R3_ChIP.bam.bai")
HCC827_ChIP_R1_bam <- bam("HCC827_R1_ChIP.bam", "HCC827_R1_ChIP.bam.bai")
HCC827_ChIP_R2_bam <- bam("HCC827_R2_ChIP.bam", "HCC827_R2_ChIP.bam.bai")
HCC827_ChIP_R3_bam <- bam("HCC827_R3_ChIP.bam", "HCC827_R3_ChIP.bam.bai")
HCC827_MET_ChIP_R1_bam <- bam("HCC827-MET_R1_ChIP.bam", "HCC827-MET_R1_ChIP.bam.bai")
HCC827_MET_ChIP_R2_bam <- bam("HCC827-MET_R2_ChIP.bam", "HCC827-MET_R2_ChIP.bam.bai")
HCC827_MET_ChIP_R3_bam <- bam("HCC827-MET_R3_ChIP.bam", "HCC827-MET_R3_ChIP.bam.bai")

A549_bam_list <- GRangesList(R1 = A549_ChIP_R1_bam,
                             R2 = A549_ChIP_R2_bam,
                             R3 = A549_ChIP_R3_bam)
A549_bam_list
HCC827_bam_list <- GRangesList(R1 = HCC827_ChIP_R1_bam,
                               R2 = HCC827_ChIP_R2_bam,
                               R3 = HCC827_ChIP_R3_bam)
HCC827_MET_bam_list <- GRangesList(R1 = HCC827_MET_ChIP_R1_bam,
                                   R2 = HCC827_MET_ChIP_R2_bam,
                                   R3 = HCC827_MET_ChIP_R3_bam)  

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")

e.score_distribution("Complete table of all targets.txt",
                     c("NRAS", "KRAS","TP53", "ERBB2", "BRCA1", "MET", "BRAF"),
                     c("KIT","PDGFRA", "ROS1", "RET", "APC", "ALK", "BRCA2"),
                     HCC827_bam_list, "HCC827")
ggsave(filename = "H3K36me3 distribution HCC827.png",
       width = 9137, height = 5812, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 1",
       dpi = 1200,
       device = "png")
e.score_distribution("Complete table of all targets.txt",
                     c("NRAS", "KRAS","TP53", "ERBB2", "BRCA1"
                       , "APC", "EGFR", "MET", "BRAF"),
                     c("KIT","PDGFRA", "ROS1", "RET", "ALK", "BRCA2"),
                     A549_bam_list, "A549")
ggsave(filename = "H3K36me3 distribution A549.png",
       width = 9137, height = 5812, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 1",
       dpi = 1200,
       device = "png")
e.score_distribution("Complete table of all targets.txt",
                     c("NRAS", "KRAS","TP53", "ERBB2", "BRCA1"
                       , "EGFR", "MET", "BRAF"),
                     c("KIT","PDGFRA", "ROS1", "RET", "ALK", "BRCA2", "APC"),
                     HCC827_MET_bam_list, "HCC827-MET")
ggsave(filename = "H3K36me3 distribution HCC827-MET.png",
       width = 9137, height = 5812, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 1",
       dpi = 1200,
       device = "png")

coverages("HCC827.bedgraph", z = "AVENIO_genes.txt")

combined_bamlist <- GRangesList(R1 = A549_ChIP_R1_bam,
                                R2 = A549_ChIP_R2_bam,
                                R3 = A549_ChIP_R3_bam,
                                R4 = HCC827_ChIP_R1_bam,
                                R5 = HCC827_ChIP_R2_bam,
                                R6 = HCC827_ChIP_R3_bam,
                                R7 = HCC827_MET_ChIP_R1_bam,
                                R8 = HCC827_MET_ChIP_R2_bam,
                                R9 = HCC827_MET_ChIP_R3_bam)
e.score_distribution("Complete table of all targets.txt",
                     c("NRAS", "KRAS","TP53", "ERBB2", "BRCA1"
                       ,"MET", "BRAF", "BRCA2"),
                     c("KIT","PDGFRA", "ROS1", "RET"),
                     combined_bamlist, "A549, HCC827 and HCC827-MET")
ggsave(filename = "H3K36me3 distribution all cell lines.png",
       width = 9137, height = 5812, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 1",
       dpi = 1200,
       device = "png")

A549_enrichment <- e.score(A549_ChIP, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
HCC827_enrichment <- e.score(HCC827_ChIP, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
clone3_enrichment <- e.score(Clone3_ChIP, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))

A549_R1_enrichment <- e.score(A549_ChIP_R1, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
A549_R2_enrichment <- e.score(A549_ChIP_R2, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
A549_R3_enrichment <- e.score(A549_ChIP_R3, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_R1_enrichment <- e.score(HCC827_ChIP_R1, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_R2_enrichment <- e.score(HCC827_ChIP_R2, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_R3_enrichment <- e.score(HCC827_ChIP_R3, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_MET_R1_enrichment <- e.score(HCC827_MET_ChIP_R1, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_MET_R2_enrichment <- e.score(HCC827_MET_ChIP_R2, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
HCC827_MET_R3_enrichment <- e.score(HCC827_MET_ChIP_R3, "Coverage of AVENIO genes.txt", gr("AVENIO_genes.txt"))
bound_cell_ChIP <- bind_samples(list(A549_R1_enrichment,
                                     A549_R2_enrichment,
                                     A549_R3_enrichment,
                                     HCC827_R1_enrichment,
                                     HCC827_R2_enrichment,
                                     HCC827_R3_enrichment,
                                     HCC827_MET_R1_enrichment,
                                     HCC827_MET_R2_enrichment,
                                     HCC827_MET_R3_enrichment),
                                c("A549_R1", "A549_R2", "A549_R3",
                                  "HCC827_R1", "HCC827_R2", "HCC827_R3",
                                  "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3"),
                                c(rep("A549", 3), rep("HCC827", 3), rep("HCC827-MET",3)))
umap_ChIP_seq(bound_cell_ChIP)
ggsave(filename = "ChIP-seq UMAP.png",
       width = 7437, height = 7437, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables",
       dpi = 1200,
       device = "png")

mean_A549_enrichment <- average_enrichment(list(A549_R1_enrichment, 
                                                A549_R2_enrichment, 
                                                A549_R3_enrichment))
mean_HCC827_enrichment <- average_enrichment(list(HCC827_R1_enrichment, 
                                                  HCC827_R2_enrichment, 
                                                  HCC827_R3_enrichment))
mean_HCC827_MET_enrichment <- average_enrichment(list(HCC827_MET_R1_enrichment, 
                                                      HCC827_MET_R2_enrichment, 
                                                      HCC827_MET_R3_enrichment))
mean_A549_enrichment
mean_HCC827_enrichment
mean_HCC827_MET_enrichment
nrow(TPM)

gg_enrichment(A549_enrichment, "A549 H3K36me3 ChIP enrichment")
gg_enrichment(HCC827_enrichment, "HCC827 H3K36me3 ChIP enrichment")
gg_enrichment(clone3_enrichment, "Clone3 H3K36me3 ChIP enrichment")

gg_enrichment(A549_R1_enrichment, "A549 R1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(A549_R2_enrichment, "A549 R2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(A549_R3_enrichment, "A549 R3 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_R1_enrichment, "HCC827 R1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_R2_enrichment, "HCC827 R2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_R3_enrichment, "HCC827 R3 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_MET_R1_enrichment, "HCC827-MET R1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_MET_R2_enrichment, "HCC827-MET R2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_MET_R3_enrichment, "HCC827-MET R3 H3K36me3 ChIP enrichment", bad)




gg_enrichment(A549_enrichment, "A549 H3K36me3 ChIP enrichment", bad)
gg_enrichment(HCC827_enrichment, "HCC827 H3K36me3 ChIP enrichment", bad)
gg_enrichment(clone3_enrichment, "Clone3 H3K36me3 ChIP enrichment", bad)

-log10(0.05)

setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
TPM_AVENIO <- RNA_TPM("Gene abundance 03112022.txt",c("log2_A549_R1",
                                                      "log2_A549_R2",
                                                      "log2_A549_R3",
                                                      "log2_HCC827_R1",
                                                      "log2_HCC827_R2",
                                                      "log2_HCC827_R3",
                                                      "log2_HCC827-MET_R1",
                                                      "log2_HCC827-MET_R2",
                                                      "log2_HCC827-MET_R3"))

TPM_AVENIO
TPM <- read.table("Gene abundance 03112022.txt", header = T)
colnames(TPM)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")

TPM_reduced <- TPM[TPM$SYMBOL %in% grs$SYMBOL,]
colnames(TPM_reduced)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
HCC827_vs_A549 <- dif_gene_express(TPM_reduced, c("HCC827_", "A549_"))
HCC827_vs_HCC827_MET <- dif_gene_express(TPM_reduced, c("HCC827_", "HCC827-MET_"))
`%ni%` <- Negate(`%in%`)
HCC827_vs_A549 <- HCC827_vs_A549 %>% dplyr::filter(gene %ni% bad)
HCC827_vs_HCC827_MET <- HCC827_vs_HCC827_MET %>% dplyr::filter(gene %ni% bad)

volcano_dif(HCC827_vs_A549, "HCC827 compared to A549")
volcano_dif(HCC827_vs_HCC827_MET, "HCC827 compared to HCC827-MET")
umap_RNA_seq(TPM,c("A549_", "HCC827_", "HCC827-MET_"))
mean_A549_TPM <- RNA_mean(TPM_AVENIO, "log2_A549_")
mean_HCC827_TPM <- RNA_mean(TPM_AVENIO, "log2_HCC827_")
mean_HCC827_MET_TPM <- RNA_mean(TPM_AVENIO, "log2_HCC827-MET_")

mean_A549_TPM
mean_HCC827_MET_TPM
mean_HCC827_TPM

setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
RNA_plot_genes(RNA_TPM("Gene abundance 03112022.txt", c("A549_R1", "A549_R2", "A549_R3",
                                                        "HCC827_R1", "HCC827_R2", "HCC827_R3",
                                                        "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")),
               c("HCC827_R1", "HCC827_R2", "HCC827_R3"), "HCC827")



RNA_plot_genes(RNA_TPM("Gene abundance 03112022.txt", c("A549_R1", "A549_R2", "A549_R3",
                                                        "HCC827_R1", "HCC827_R2", "HCC827_R3",
                                                        "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")),
               c("A549_R1", "A549_R2", "A549_R3"), "A549")

RNA_plot_genes(RNA_TPM("Gene abundance 03112022.txt", c("A549_R1", "A549_R2", "A549_R3",
                                                        "HCC827_R1", "HCC827_R2", "HCC827_R3",
                                                        "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")),
               c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3"), "HCC827-MET")



RNA_plot_genes(RNA_TPM("Gene abundance 03112022.txt", c("A549_R1", "A549_R2", "A549_R3",
                                                        "HCC827_R1", "HCC827_R2", "HCC827_R3",
                                                        "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")),
               c("A549_R1", "A549_R2", "A549_R3",
                 "HCC827_R1", "HCC827_R2", "HCC827_R3",
                 "HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3"),
               "A549, HCC827, HCC827-MET")
ggsave(filename = "RNA expression 12 genes all cell lines.png",
       width = 12275, height = 7800, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary",
       dpi = 1200,
       device = "png")

versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2)
versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2, a = T,)
versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2, bad, T)
versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2, bad)
library(writexl)
write_xlsx(full_join(versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2),
                     versus(mean_A549_enrichment, mean_A549_TPM, "A549", 0.2, bad),
                     by = "genes"), 
           "A549 inclusion-exclusion.xlsx")
write_xlsx(full_join(versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2),
                     versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2, bad),
                     by = "genes"), 
           "HCC827 inclusion-exclusion.xlsx")
write_xlsx(full_join(versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2),
                     versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2, bad),
                     by = "genes"), 
           "HCC827-MET inclusion-exclusion.xlsx")

?full_join
versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2, "EGFR")
versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2, a = T,g = "EGFR")
versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2, bad, T)
versus(mean_HCC827_enrichment, mean_HCC827_TPM, "HCC827", 0.2, bad)

versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2, g = c("EGFR", "MET"))
versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2, a = T, g = c("EGFR", "MET"))
versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2, bad, a = T)
versus(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM, "HCC827-MET", 0.2, bad)

ggsave(filename = "A549 ChIP vs. RNA.png",
       width = 7437, height = 7437, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/figure 1",
       dpi = 1200,
       device = "png")

ChIPcorr_cell(mean_A549_enrichment, mean_HCC827_enrichment, "A549", "HCC827",bad)
ChIPcorr_cell(mean_A549_enrichment, mean_HCC827_enrichment, "A549", "HCC827")
ChIPcorr_cell(mean_A549_enrichment, mean_HCC827_enrichment, "A549", "HCC827",bad, z = 1)
ChIPcorr_cell(mean_A549_enrichment, mean_HCC827_enrichment, "A549", "HCC827", z = 1)

ChIPcorr_cell(mean_HCC827_MET_enrichment, mean_HCC827_enrichment, "HCC827-MET", "HCC827",bad)
ChIPcorr_cell(mean_HCC827_MET_enrichment, mean_HCC827_enrichment, "HCC827-MET", "HCC827")
ChIPcorr_cell(mean_HCC827_MET_enrichment, mean_HCC827_enrichment, "HCC827-MET", "HCC827",bad, z = 1)
ChIPcorr_cell(mean_HCC827_MET_enrichment, mean_HCC827_enrichment, "HCC827-MET", "HCC827", z = 1)




venn(A549_enrichment, HCC827_enrichment, TPM_AVENIO, "log2.HCC827", "log2.A549",
     bad, d = 1)
venn(clone3_enrichment, HCC827_enrichment, TPM_AVENIO, "log2.HCC827", "log2.Clone3",
     bad, d = 1)

bar_overlap(mean_A549_enrichment, mean_HCC827_enrichment, "HCC827_", "A549_",
            bad, d = 1)
bar_overlap(mean_HCC827_MET_enrichment, mean_HCC827_enrichment, "HCC827_", "HCC827-MET_",
            bad, d = 1)

get_ROC_data(mean_A549_enrichment, mean_A549_TPM,
             r = 0.2,
             b = bad)
get_ROC_data(mean_HCC827_enrichment, mean_HCC827_TPM,
             r = 0.2,
             b = bad)
get_ROC_data(mean_HCC827_MET_enrichment, mean_HCC827_MET_TPM,
             r = 0.2,
             b = bad)

gc()
setwd("D:/Lung cancer cfChIP/Deduped")
Adeno1 <- gene_count("A-1279-cfChIP.bam",grs)
Adeno2 <- gene_count("B-1288-cfChIP.bam", grs)
Adeno3 <- gene_count("C-1475-cfChIP.bam", grs)
Adeno4 <- gene_count("D-1578-cfChIP.bam", grs)
Plano1 <- gene_count("E-439-cfChIP.bam", grs)
Plano2 <- gene_count("F-1449-cfChIP.bam", grs)
Plano3 <- gene_count("I-645-cfChIP.bam",grs)
Plano4 <- gene_count("J-1663-cfChIP.bam",grs)
SCLC1 <- gene_count("G-514-cfChIP.bam", grs)
SCLC2 <- gene_count("H-1169-cfChIP.bam", grs)
SCLC3 <- gene_count("K-440-cfChIP.bam", grs)
SCLC4 <- gene_count("L-1100-cfChIP.bam", grs)
setwd("D:/Healthy cfChIP/Deduped")
healthy1 <- gene_count("Rask_kontrol1_cfChIP.bam",grs)
healthy2 <- gene_count("Rask_kontrol2_cfChIP.bam",grs)
healthy3 <- gene_count("Rask_kontrol3_cfChIP.bam",grs)
healthy4 <- gene_count("Rask_kontrol4_cfChIP.bam",grs)
healthies <- list(healthy1,healthy2,healthy3,healthy4)
healthy_reads <- healthy(healthies)
healthy_reads
Plano1
df_genecounts <- data.frame(genes = Adeno1$genes,
                            Adeno1 = Adeno1$readcounts,
                            Adeno2 = Adeno2$readcounts,
                            Adeno3 = Adeno3$readcounts,
                            Adeno4 = Adeno4$readcounts,
                            Plano1 = Plano1$readcounts,
                            Plano2 = Plano2$readcounts,
                            Plano3 = Plano3$readcounts,
                            Plano4 = Plano4$readcounts,
                            SCLC1 = SCLC1$readcounts,
                            SCLC2 = SCLC2$readcounts,
                            SCLC3 = SCLC3$readcounts,
                            SCLC4 = SCLC4$readcounts,
                            Healthy1 = healthy1$readcounts,
                            Healthy2 = healthy2$readcounts,
                            Healthy3 = healthy3$readcounts,
                            Healthy4 = healthy4$readcounts,
                            A549 = A549_ChIP$readcounts,
                            HCC827 = HCC827_ChIP$readcounts,
                            HCC827_MET_ChIP = Clone3_ChIP$readcounts)
library(writexl)
write_xlsx(df_genecounts, "Raw readcounts data.xlsx")

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
Adeno1_enrichment <- e.score(Adeno1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Adeno2_enrichment <- e.score(Adeno2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Adeno3_enrichment <- e.score(Adeno3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Adeno4_enrichment <- e.score(Adeno4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Plano1_enrichment <- e.score(Plano1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Plano2_enrichment <- e.score(Plano2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Plano3_enrichment <- e.score(Plano3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
Plano4_enrichment <- e.score(Plano4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
SCLC1_enrichment <- e.score(SCLC1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
SCLC2_enrichment <- e.score(SCLC2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
SCLC3_enrichment <- e.score(SCLC3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
SCLC4_enrichment <- e.score(SCLC4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"),healthy_reads)
healthy1_enrichment <- e.score(healthy1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
healthy2_enrichment <- e.score(healthy2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
healthy3_enrichment <- e.score(healthy3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
healthy4_enrichment <- e.score(healthy4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
list_healthy_enrichment <- list(healthy1_enrichment,healthy2_enrichment,healthy3_enrichment,healthy4_enrichment)
average_enrichment(list_healthy_enrichment)

Adeno1_enrichment_unnormalized <- e.score(Adeno1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno2_enrichment_unnormalized <- e.score(Adeno2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno3_enrichment_unnormalized <- e.score(Adeno3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno4_enrichment_unnormalized <- e.score(Adeno4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano1_enrichment_unnormalized <- e.score(Plano1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano2_enrichment_unnormalized <- e.score(Plano2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano3_enrichment_unnormalized <- e.score(Plano3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano4_enrichment_unnormalized <- e.score(Plano4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))


all_enrichment <- list(Adeno1_enrichment,Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment,
                       Plano1_enrichment,Plano2_enrichment,Plano3_enrichment, Plano4_enrichment,
                       SCLC1_enrichment,SCLC2_enrichment, SCLC3_enrichment, SCLC4_enrichment,
                       healthy1_enrichment,healthy2_enrichment,healthy3_enrichment,healthy4_enrichment)

all_samples <- bind_samples(list(Adeno1_enrichment, Adeno2_enrichment, 
                                 Adeno3_enrichment, Adeno4_enrichment,
                                 Plano1_enrichment, Plano2_enrichment,
                                 Plano3_enrichment, Plano4_enrichment,
                                 SCLC1_enrichment, SCLC2_enrichment,
                                 SCLC3_enrichment, SCLC4_enrichment,
                                 healthy1_enrichment, healthy2_enrichment, 
                                 healthy3_enrichment, healthy4_enrichment),
                            c("Adeno1", "Adeno2", "Adeno3", "Adeno4",
                              "Squamous1", "Squamous2", "Squamous3", "Squamous4",
                              "SCLC1", "SCLC2", "SCLC3", "SCLC4",
                              "Healthy1", "Healthy2", "Healthy3", "Healthy4"),
                            c(rep("Adeno",4),rep("Squamous",4),
                              rep("SCLC",4), rep("Healthy",4)))



gg_enrichment(Adeno1_enrichment, "Adeno1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Adeno2_enrichment, "Adeno2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Adeno3_enrichment, "Adeno3 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Adeno4_enrichment, "Adeno4 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Plano1_enrichment, "Plano1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Plano2_enrichment, "Plano2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Plano3_enrichment, "Plano3 H3K36me3 ChIP enrichment", bad)
gg_enrichment(Plano4_enrichment, "Plano4 H3K36me3 ChIP enrichment", bad)
gg_enrichment(SCLC1_enrichment, "SCLC1 H3K36me3 ChIP enrichment", bad)
gg_enrichment(SCLC2_enrichment, "SCLC2 H3K36me3 ChIP enrichment", bad)
gg_enrichment(SCLC3_enrichment, "SCLC3 H3K36me3 ChIP enrichment", bad)
gg_enrichment(SCLC4_enrichment, "SCLC4 H3K36me3 ChIP enrichment", bad)
ChIPcorr(Adeno1_enrichment, Adeno2_enrichment, "NAC.1", "NAC.2",bad)
ChIPcorr(Adeno1_enrichment, Adeno3_enrichment, "NAC.1", "NAC.3",bad)
ChIPcorr(Adeno1_enrichment, Adeno4_enrichment, "NAC.1", "NAC.4",bad)
ChIPcorr(Adeno2_enrichment, Adeno3_enrichment, "NAC.2", "NAC.3",bad)
ChIPcorr(Adeno2_enrichment, Adeno4_enrichment, "NAC.2", "NAC.4",bad)
ChIPcorr(Adeno4_enrichment, Adeno3_enrichment, "NAC.4", "NAC.3",bad)




Adeno_tot <- list(Adeno1_enrichment,Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment)
Plano_tot <- list(Plano1_enrichment,Plano2_enrichment,Plano3_enrichment,Plano4_enrichment)
SCLC_tot <- list(SCLC1_enrichment,SCLC2_enrichment,SCLC3_enrichment,SCLC4_enrichment)
NSCLC_tot <- list(Adeno1_enrichment,Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment,
                  Plano1_enrichment,Plano2_enrichment,Plano3_enrichment,Plano4_enrichment)
Adeno_tot_not1 <- list(Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment)
NSCLC_tot_unnormalized <- list(Adeno1_enrichment_unnormalized,Adeno2_enrichment_unnormalized,
                               Adeno3_enrichment_unnormalized,Adeno4_enrichment_unnormalized,
                               Plano1_enrichment_unnormalized,Plano2_enrichment_unnormalized,
                               Plano3_enrichment_unnormalized,Plano4_enrichment_unnormalized)
Adeno_tot_not4 <- list(Adeno1_enrichment, Adeno2_enrichment,Adeno3_enrichment)

Adeno_tot_enrichment <- average_enrichment(Adeno_tot)
Plano_tot_enrichment <- average_enrichment(Plano_tot)
SCLC_tot_enrichment <- average_enrichment(SCLC_tot)
NSCLC_tot_enrichment <- average_enrichment(NSCLC_tot)
NSCLC_tot_unnormalized_enrichment <- average_enrichment(NSCLC_tot_unnormalized)
Adeno_tot_not1_enrichment <- average_enrichment(Adeno_tot_not1)
Adeno_tot_not4_enrichment <- average_enrichment(Adeno_tot_not4)

ChIPcorr(Adeno_tot_enrichment, Plano_tot_enrichment, "Adeno avg", "Plano avg", bad)
ChIPcorr(Adeno_tot_enrichment, SCLC_tot_enrichment, "Adeno avg", "SCLC avg", bad)
ChIPcorr(Plano_tot_enrichment, SCLC_tot_enrichment, "Plano avg", "SCLC avg", bad)
ChIPcorr(NSCLC_tot_enrichment,SCLC_tot_enrichment, "NSCLC", "SCLC", bad)
ChIPcorr_fun(NSCLC_tot_unnormalized_enrichment,average_enrichment(list_healthy_enrichment), "NSCLC", "Healthy", bad)

EGFR_mut <- list(Adeno1_enrichment, Adeno4_enrichment)
EGFR_wt <- list(Adeno2_enrichment, Adeno3_enrichment,
                Plano1_enrichment, Plano2_enrichment,
                Plano3_enrichment, Plano4_enrichment)
EGFR_mut_enrichment <- average_enrichment(EGFR_mut)
EGFR_WT_enrichment <- average_enrichment(EGFR_wt)

ChIPcorr(EGFR_mut_enrichment, EGFR_WT_enrichment, "EGFR mutated", "EGFR WT", bad)

ggsave(filename = "cfChIP-seq comparison of EGFR mutated and EGFR WT NSCLC.png",
       width = 10312, height = 7437, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables",
       dpi = 1200,
       device = "png")


Adeno_tot <- list(Adeno1_enrichment,Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment)
Plano_tot <- list(Plano1_enrichment,Plano2_enrichment,Plano3_enrichment,Plano4_enrichment)
SCLC_tot <- list(SCLC1_enrichment,SCLC2_enrichment,SCLC3_enrichment,SCLC4_enrichment)
NSCLC_tot <- list(Adeno1_enrichment,Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment,
                  Plano1_enrichment,Plano2_enrichment,Plano3_enrichment,Plano4_enrichment)
Adeno_tot_not1 <- list(Adeno2_enrichment,Adeno3_enrichment,Adeno4_enrichment)
NSCLC_tot_unnormalized <- list(Adeno1_enrichment_unnormalized,Adeno2_enrichment_unnormalized,
                               Adeno3_enrichment_unnormalized,Adeno4_enrichment_unnormalized,
                               Plano1_enrichment_unnormalized, Plano2_enrichment_unnormalized,
                               Plano2_enrichment_unnormalized, Plano3_enrichment_unnormalized)




ChIPcorr(Adeno1_enrichment, Adeno_tot_not1_enrichment, "Adenocarcinoma 1", "Average Adenocarcinoma", bad)
ChIPcorr(Adeno4_enrichment, Adeno_tot_not4_enrichment, "Adenocarcinoma 4", "Average Adenocarcinoma", bad)



ChIPcorr_fun <- function(x,y,p,q,b = NULL,z = NULL){
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
            x <- x %>% dplyr::filter(genes %in% y$genes)
        }
        if (length(y$genes) > length(x$genes)){
            y <- y %>% dplyr::filter(genes %in% x$genes)
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
        df <- df %>% dplyr::filter(genes %ni% b)
    } 
    df$reverse <- c(1:length(df$genes))
    if (is.null(z)){
        gg <- ggplot(data = df, aes(x = reverse, y = log2FC))+
            geom_point(aes(colour=log2FC), size = 4)+
            labs(title=paste("cfChIP-seq comparison of", p, "and",q), 
                 x = "Genes", y= "Relative enrichment")+
            scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="grey", 
                                   high="#6a00fc", name = "Log2 difference", 
                                   limits=c(round(min(df$log2FC))-0.5,round(max(df$log2FC))+0.5))+
            geom_abline(intercept = 0, slope = 0, color = "Black", linetype = "solid", size = 1.5)+
            geom_label_repel(
                aes(label=ifelse(genes %in% c("TNFRSF21", "SLPI", "EGFR", "MET"), as.character(paste(genes)),"")),
                segment.color="#6a00fc",
                box.padding   = 0.35,
                point.padding = 0.1,
                label.padding = 0.2,
                color="#6a00fc",
                label.size = NA,
                parse = F,
                size = 3.5,
                max.overlaps = 100,
                fill = alpha(c("white"),0))+
            theme_bw()+
            geom_label(x = 40, y = (0.5), 
                       label = paste("Upregulated in",p),
                       color = "#6a00fc", label.size = 0, size = 5)+
            geom_label(x = 120, y = (-0.5), label = paste("Upregulated in",q),
                       color = "#ffa10c", label.size = 0, size = 5)+
            th
    }
    else{
        gg <- ggplot(data = df, aes(x = reverse, y = log2FC))+
            geom_point(aes(colour=log2FC), size = 4)+
            labs(title=paste("cfChIP-seq comparison of", p, "and",q), 
                 x = "Genes", y= "Relative enrichment",
                 caption = paste("Cut-off at Log2 difference = +/-",z))+
            scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="grey", 
                                   high="#6a00fc", name = "Log2 difference", 
                                   limits=c(round(min(df$log2FC))-0.5,round(max(df$log2FC))+0.5))+
            geom_abline(intercept = z, slope = 0, color = "Black", linetype = "dashed", size = 1)+
            geom_abline(intercept = -z, slope = 0, color = "Black", linetype = "dashed", size = 1)+
            geom_abline(intercept = 0, slope = 0, color = "Black", linetype = "solid", size = 1.5)+
            geom_label_repel(
                aes(label=ifelse(genes %in% c("TNFRSF21", "SLPI", "EGFR", "MET", "FBXL7", "PDZRN3",
                                              "ADAMTS12", "ROS1", "SOX9"), as.character(paste(genes)),"")),
                box.padding   = 0.35,
                point.padding = 0.1,
                label.padding = 0.2,
                nudge_x = 0,
                nudge_y= 0,
                segment.color="#6a00fc",
                color = "#6a00fc",
                label.size = NA,
                parse = F,
                size = 3.5,
                max.overlaps = 100,
                fill = alpha(c("white"),0))+
            theme_bw()+
            geom_label(x = 20, y = (z-0.25), 
                       label = paste("Upregulated in",p),
                       color = "#6a00fc", label.size = 0, size = 5)+
            geom_label(x = 120, y = (-z+0.25), label = paste("Upregulated in",q),
                       color = "#ffa10c", label.size = 0, size = 5)+
            th
    }
    return(gg)
}


setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")

cover_hist("Coverage of AVENIO genes.txt")
TSS_hist("Complete table of all targets.txt")
number_hist("Complete table of all targets.txt")
RNA_HCC827(TPM_AVENIO)



cell_samples <- bind_samples(list(HCC827_enrichment, clone3_enrichment, A549_enrichment),
                             c("HCC827_PAR", "HCC827_clone3", "A549"),
                             c("HCC827", "HCC827", "A549"))

histology_samples <- bind_samples(list(Adeno1_enrichment, Adeno2_enrichment, 
                                       Adeno3_enrichment, Adeno4_enrichment,
                                       Plano1_enrichment, Plano2_enrichment,
                                       Plano3_enrichment, Plano4_enrichment,
                                       SCLC1_enrichment, SCLC2_enrichment,
                                       SCLC3_enrichment, SCLC4_enrichment),
                                  c("Adeno1", "Adeno2", "Adeno3", "Adeno4",
                                    "Squamous1", "Squamous2", "Squamous3", "Squamous4",
                                    "SCLC1", "SCLC2", "SCLC3", "SCLC4"),
                                  c(rep("Adeno",4),rep("Squamous",4),rep("SCLC",4)))
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Healthy control cfChIP")

fragment_length(c("Rask_kontrol1_cfChIP.bam","Rask_kontrol2_cfChIP.bam",
                  "Rask_kontrol3_cfChIP.bam","Rask_kontrol4_cfChIP.bam"),
                c("Rask_kontrol1_input.bam","Rask_kontrol2_input.bam",
                  "Rask_kontrol3_input.bam","Rask_kontrol4_input.bam"),
                c("cfChIP", "Input"), "Fragment lengths healthy input and cfChIP")

fragment_length(c("A549_ChIP.bam","HCC827ChIP.bam",
                  "HCC827_ER_klon3_ChIP.bam"),
                c("A549_input.bam","HCC827input.bam",
                  "HCC827_ER_klon3_input.bam"),
                c("ChIP", "Input"), "Fragment lengths cell input and ChIP")

fragment_length(c("A-1279-cfChIP.bam","B-1288-cfChIP.bam",
                  "C-1475-cfChIP.bam","D-1578-cfChIP.bam",
                  "E-439-cfChIP.bam","F-1449-cfChIP.bam",
                  "I-645-cfChIP.bam","J-1663-cfChIP.bam",
                  "G-514-cfChIP.bam","H-1169-cfChIP.bam",
                  "K-440-cfChIP.bam","L-1100-cfChIP.bam"),
                c("Rask_kontrol1_cfChIP.bam","Rask_kontrol2_cfChIP.bam",
                  "Rask_kontrol3_cfChIP.bam","Rask_kontrol4_cfChIP.bam"),
                c("Cancer", "Healthy"), "Fragment lengths cancer and healthy cfChIP samples")


PCA(histology_samples, "PCA of different tumor subtypes")
tSNE(all_samples, "tSNE plot of all samples sorted")

metrics("A-1279-cfChIP.bam", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("B-1288-cfChIP.bam", "B-1288-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("C-1475-cfChIP.bam", "C-1475-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("D-1578-cfChIP.bam", "D-1578-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("E-439-cfChIP.bam", "E-439-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("F-1449-cfChIP.bam", "F-1449-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("G-514-cfChIP.bam", "G-514-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("H-1169-cfChIP.bam", "H-1169-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("I-645-cfChIP.bam","I-645-cfChIP.bam.bai","AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("J-1663-cfChIP.bam","J-1663-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("K-440-cfChIP.bam", "K-440-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)
metrics("L-1100-cfChIP.bam", "L-1100-cfChIP.bam.bai", "AVENIO_genes.txt", "Coverage of AVENIO genes.txt",
        "C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq", healthy_reads)

gen <- "CACNA1E"

NSCLC_tot_enrichment %>% filter(genes == gen)
SCLC_tot_enrichment %>% filter(genes == gen)