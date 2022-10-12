setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article")
SCLC <- read.table(file = 'SCLC.tsv', sep = '\t', header = TRUE)
SCLC
# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052
library(dplyr)
NSCLC <- read.table(file = "RPKM_NSCLC.txt", header = T)
NSCLC <- NSCLC %>% filter(Gene_type == "protein_coding")
NSCLC
# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190139

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
AVENIO <- gr("AVENIO_genes.txt")

SCLC_avenio <- SCLC %>% filter(gene %in% AVENIO$SYMBOL)
NSCLC_avenio <- NSCLC %>% filter(Gene_symbol %in% AVENIO$SYMBOL)
NSCLC_avenio

NSCLC_dd <- NSCLC %>% filter(Gene_symbol %in% 
                               c("EGFR", "CSMD1", "CYBB", "DMD", "KLHL31", "KPRP", "POLE"))
NSCLC_dd$mean <- rowMeans(NSCLC_dd[5:50])
NSCLC_dd$sd <- apply(NSCLC_dd[5:50],1,sd)
NSCLC_dd <- NSCLC_dd %>% select(Gene_symbol, mean, sd)
NSCLC_dd
library(writexl)
write_xlsx(NSCLC_dd, "NSCLC ddpcr genes.xlsx")
