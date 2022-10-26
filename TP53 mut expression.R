#Baseret på denne publikation https://pubmed.ncbi.nlm.nih.gov/32015526/
#Og data er hentet fra https://www.cbioportal.org/study/summary?id=luad_oncosg_2020
setwd("C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article")
TP53_all <- read.table("TP53__mRNA_expression_all.txt", header = T)
TP53_mut <- read.table("TP53__mRNA_expression_tp53mut.txt", header = T)
library(dplyr)
colnames(TP53_all)
TP53_WT <- TP53_all %>% 
    filter(!PatientID %in% TP53_mut$PatientID) %>% 
    filter(!is.na(mRNA_expression))
TP53_WT <- TP53_WT %>% 
    mutate(TP53_status = rep("WT", nrow(TP53_WT)))
length(TP53_WT$StudyID)

TP53_mut <- TP53_mut %>%
    filter(!is.na(mRNA_expression))
TP53_mut <- TP53_mut %>% 
    mutate(TP53_status = rep("mutant", nrow(TP53_mut)))
length(TP53_mut$StudyID)

TP53_expression <- bind_rows(TP53_WT, TP53_mut)

library(ggplot2)
library(ggpubr)

TP53_expression %>% 
    group_by(TP53_status) %>% 
    ggplot(aes(x=TP53_status, y=mRNA_expression, fill=TP53_status)) + 
    geom_boxplot(outlier.alpha = 0)+
    geom_jitter(alpha = 0.5, width = 0.1)+
    stat_compare_means(method = "t.test",
                       label.x = 1.5, 
                       label.y = c(5.5,4.5,5),
                       comparisons = list(c("mutant", "WT")))+
    scale_color_manual("TP53 status", values = c("#ffa10c","#6a00fc"))+
    scale_fill_manual("TP53 status", values = c("#ffa10c","#6a00fc"))+
    labs(x = "",
         y = "TP53 mRNA z-score")+
    theme_bw()+
    scale_x_discrete(labels = c("mutant" = "TP53-mutant \n n = 58",
                                "WT" = "TP53-WT \n n = 111"))+
    th+
    theme(axis.text.x = element_text(face = "bold"))

ggsave(filename = "TP53 expression in tumor samples.png",
       width = 7625, height = 5963, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables",
       dpi = 1200,
       device = "png")
