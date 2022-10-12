setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/Depmap data")
library(dplyr)
library(writexl)
`%ni%` <- Negate(`%in%`)
ABCC5 <- read.csv("ABCC5.csv")
DSC3 <- read.csv("DSC3.csv")
EGFR <- read.csv("EGFR.csv")
ERBB2 <- read.csv("ERBB2.csv")
HECW1 <- read.csv("HECW1.csv")
KIF19 <- read.csv("KIF19.csv")
KRAS <- read.csv("KRAS.csv")
PAX6 <- read.csv("PAX6.csv")
PIK3CA <- read.csv("PIK3CA.csv")
RIN3 <- read.csv("RIN3.csv")
SMAD4 <- read.csv("SMAD4.csv")
ZFPM2 <- read.csv("ZFPM2.csv")
MAP7D3 <- read.csv("MAP7D3.csv")

GRIK3 <- read.csv("GRIK3.csv")
CRMP1 <- read.csv("CRMP1.csv")
MYT1L <- read.csv("MYT1L.csv")
RALYL <- read.csv("RALYL.csv")
PDZRN3 <- read.csv("PDZRN3.csv")
MAP2 <- read.csv("MAP2.csv")
CACNA1E <- read.csv("CACNA1E.csv")


#Adeno
ABCC5_Adeno <- ABCC5 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(ABCC5.log2.TPM.1..Expression.22Q1.Public)
DSC3_Adeno <- DSC3 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(DSC3.log2.TPM.1..Expression.22Q1.Public)
EGFR_Adeno <- EGFR %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(EGFR.log2.TPM.1..Expression.22Q1.Public)
ERBB2_Adeno <- ERBB2 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(ERBB2.log2.TPM.1..Expression.22Q1.Public)
HECW1_Adeno <- HECW1 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(HECW1.log2.TPM.1..Expression.22Q1.Public)
KIF19_Adeno <- KIF19 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(KIF19.log2.TPM.1..Expression.22Q1.Public)
KRAS_Adeno <- KRAS %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(KRAS.log2.TPM.1..Expression.22Q1.Public)
PAX6_Adeno <- PAX6 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(PAX6.log2.TPM.1..Expression.22Q1.Public)
PIK3CA_Adeno <- PIK3CA %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(PIK3CA.log2.TPM.1..Expression.22Q1.Public)
RIN3_Adeno <- RIN3 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(RIN3.log2.TPM.1..Expression.22Q1.Public)
SMAD4_Adeno <- SMAD4 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma") %>% dplyr::select(SMAD4.log2.TPM.1..Expression.22Q1.Public)

#Squamous
ABCC5_Squamous <- ABCC5 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(ABCC5.log2.TPM.1..Expression.22Q1.Public)
DSC3_Squamous <- DSC3 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(DSC3.log2.TPM.1..Expression.22Q1.Public)
EGFR_Squamous <- EGFR %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(EGFR.log2.TPM.1..Expression.22Q1.Public)
ERBB2_Squamous <- ERBB2 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(ERBB2.log2.TPM.1..Expression.22Q1.Public)
HECW1_Squamous <- HECW1 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(HECW1.log2.TPM.1..Expression.22Q1.Public)
KIF19_Squamous <- KIF19 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(KIF19.log2.TPM.1..Expression.22Q1.Public)
KRAS_Squamous <- KRAS %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(KRAS.log2.TPM.1..Expression.22Q1.Public)
PAX6_Squamous <- PAX6 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(PAX6.log2.TPM.1..Expression.22Q1.Public)
PIK3CA_Squamous <- PIK3CA %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(PIK3CA.log2.TPM.1..Expression.22Q1.Public)
RIN3_Squamous <- RIN3 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(RIN3.log2.TPM.1..Expression.22Q1.Public)
SMAD4_Squamous <- SMAD4 %>% filter(all.Disease.Subtype == "Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma") %>% dplyr::select(SMAD4.log2.TPM.1..Expression.22Q1.Public)

#SCLC
ABCC5_SCLC <- ABCC5 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(ABCC5.log2.TPM.1..Expression.22Q1.Public)
DSC3_SCLC <- DSC3 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(DSC3.log2.TPM.1..Expression.22Q1.Public)
EGFR_SCLC <- EGFR %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(EGFR.log2.TPM.1..Expression.22Q1.Public)
ERBB2_SCLC <- ERBB2 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(ERBB2.log2.TPM.1..Expression.22Q1.Public)
HECW1_SCLC <- HECW1 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(HECW1.log2.TPM.1..Expression.22Q1.Public)
KIF19_SCLC <- KIF19 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(KIF19.log2.TPM.1..Expression.22Q1.Public)
KRAS_SCLC <- KRAS %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(KRAS.log2.TPM.1..Expression.22Q1.Public)
PAX6_SCLC <- PAX6 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(PAX6.log2.TPM.1..Expression.22Q1.Public)
PIK3CA_SCLC <- PIK3CA %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(PIK3CA.log2.TPM.1..Expression.22Q1.Public)
RIN3_SCLC <- RIN3 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(RIN3.log2.TPM.1..Expression.22Q1.Public)
SMAD4_SCLC <- SMAD4 %>% filter(all.Disease.Subtype == "Small Cell Lung Cancer (SCLC)") %>% dplyr::select(SMAD4.log2.TPM.1..Expression.22Q1.Public)
GRIK3_SCLC <- GRIK3 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(GRIK3.log2.TPM.1..Expression.22Q1.Public)
CRMP1_SCLC <- CRMP1 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(CRMP1.log2.TPM.1..Expression.22Q1.Public)
MYT1L_SCLC <- MYT1L %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(MYT1L.log2.TPM.1..Expression.22Q1.Public)
RALYL_SCLC <- RALYL %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(RALYL.log2.TPM.1..Expression.22Q1.Public)
PDZRN3_SCLC <- PDZRN3 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(PDZRN3.log2.TPM.1..Expression.22Q1.Public)
MAP2_SCLC <- MAP2 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(MAP2.log2.TPM.1..Expression.22Q1.Public)
CACNA1E_SCLC <- CACNA1E %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(CACNA1E.log2.TPM.1..Expression.22Q1.Public)
ZFPM2_SCLC <- ZFPM2 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(ZFPM2.log2.TPM.1..Expression.22Q2.Public)
MAP7D3_SCLC <- MAP7D3 %>% filter(X2.Lineage == "SCLC") %>% dplyr::select(MAP7D3.log2.TPM.1..Expression.22Q2.Public)

non_NSCLC <- c("Small Cell Lung Cancer (SCLC)","Mesothelioma","Carcinoid" )
#NSCLC


EGFR_NSCLC <- EGFR %>% filter(all.Disease.Subtype %ni% non_NSCLC) %>% dplyr::select(EGFR.log2.TPM.1..Expression.22Q1.Public)
KIF19_NSCLC <- KIF19 %>% filter(all.Disease.Subtype %ni% non_NSCLC) %>% dplyr::select(KIF19.log2.TPM.1..Expression.22Q1.Public)
PAX6_NSCLC <- PAX6 %>% filter(all.Disease.Subtype %ni% non_NSCLC) %>% dplyr::select(PAX6.log2.TPM.1..Expression.22Q1.Public)
RIN3_NSCLC <- RIN3 %>% filter(all.Disease.Subtype %ni% non_NSCLC) %>% dplyr::select(RIN3.log2.TPM.1..Expression.22Q1.Public)
SMAD4_NSCLC <- SMAD4 %>% filter(all.Disease.Subtype %ni% non_NSCLC) %>% dplyr::select(SMAD4.log2.TPM.1..Expression.22Q1.Public)
GRIK3_NSCLC <- GRIK3 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(GRIK3.log2.TPM.1..Expression.22Q1.Public)
CRMP1_NSCLC <- CRMP1 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(CRMP1.log2.TPM.1..Expression.22Q1.Public)
MYT1L_NSCLC <- MYT1L %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(MYT1L.log2.TPM.1..Expression.22Q1.Public)
RALYL_NSCLC <- RALYL %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(RALYL.log2.TPM.1..Expression.22Q1.Public)
PDZRN3_NSCLC <- PDZRN3 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(PDZRN3.log2.TPM.1..Expression.22Q1.Public)
MAP2_NSCLC <- MAP2 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(MAP2.log2.TPM.1..Expression.22Q1.Public)
CACNA1E_NSCLC <- CACNA1E %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(CACNA1E.log2.TPM.1..Expression.22Q1.Public)
ZFPM2_NSCLC <- ZFPM2 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(ZFPM2.log2.TPM.1..Expression.22Q2.Public)
MAP7D3_NSCLC <- MAP7D3 %>% filter(X2.Lineage == "NSCLC") %>% dplyr::select(MAP7D3.log2.TPM.1..Expression.22Q2.Public)

unique(SMAD4$all.Disease.Subtype)

Adeno_df <- data.frame(ABCC5 = ABCC5_Adeno$ABCC5.log2.TPM.1..Expression.22Q1.Public,
                       DSC = DSC3_Adeno$DSC3.log2.TPM.1..Expression.22Q1.Public,
                       EGFR = EGFR_Adeno$EGFR.log2.TPM.1..Expression.22Q1.Public,
                       ERBB2 = ERBB2_Adeno$ERBB2.log2.TPM.1..Expression.22Q1.Public,
                       HECW1 = HECW1_Adeno$HECW1.log2.TPM.1..Expression.22Q1.Public,
                       KIF19 = KIF19_Adeno$KIF19.log2.TPM.1..Expression.22Q1.Public,
                       KRAS = KRAS_Adeno$KRAS.log2.TPM.1..Expression.22Q1.Public,
                       PAX6 = PAX6_Adeno$PAX6.log2.TPM.1..Expression.22Q1.Public,
                       PIK3CA = PIK3CA_Adeno$PIK3CA.log2.TPM.1..Expression.22Q1.Public,
                       RIN3 = RIN3_Adeno$RIN3.log2.TPM.1..Expression.22Q1.Public,
                       SMAD4 = SMAD4_Adeno$SMAD4.log2.TPM.1..Expression.22Q1.Public)

Squamous_df <- data.frame(ABCC5 = ABCC5_Squamous$ABCC5.log2.TPM.1..Expression.22Q1.Public,
                       DSC = DSC3_Squamous$DSC3.log2.TPM.1..Expression.22Q1.Public,
                       EGFR = EGFR_Squamous$EGFR.log2.TPM.1..Expression.22Q1.Public,
                       ERBB2 = ERBB2_Squamous$ERBB2.log2.TPM.1..Expression.22Q1.Public,
                       HECW1 = HECW1_Squamous$HECW1.log2.TPM.1..Expression.22Q1.Public,
                       KIF19 = KIF19_Squamous$KIF19.log2.TPM.1..Expression.22Q1.Public,
                       KRAS = KRAS_Squamous$KRAS.log2.TPM.1..Expression.22Q1.Public,
                       PAX6 = PAX6_Squamous$PAX6.log2.TPM.1..Expression.22Q1.Public,
                       PIK3CA = PIK3CA_Squamous$PIK3CA.log2.TPM.1..Expression.22Q1.Public,
                       RIN3 = RIN3_Squamous$RIN3.log2.TPM.1..Expression.22Q1.Public,
                       SMAD4 = SMAD4_Squamous$SMAD4.log2.TPM.1..Expression.22Q1.Public)

SCLC_df <- data.frame(ABCC5 = ABCC5_SCLC$ABCC5.log2.TPM.1..Expression.22Q1.Public,
                       DSC = DSC3_SCLC$DSC3.log2.TPM.1..Expression.22Q1.Public,
                       EGFR = EGFR_SCLC$EGFR.log2.TPM.1..Expression.22Q1.Public,
                       ERBB2 = ERBB2_SCLC$ERBB2.log2.TPM.1..Expression.22Q1.Public,
                       HECW1 = HECW1_SCLC$HECW1.log2.TPM.1..Expression.22Q1.Public,
                       KIF19 = KIF19_SCLC$KIF19.log2.TPM.1..Expression.22Q1.Public,
                       KRAS = KRAS_SCLC$KRAS.log2.TPM.1..Expression.22Q1.Public,
                       PAX6 = PAX6_SCLC$PAX6.log2.TPM.1..Expression.22Q1.Public,
                       PIK3CA = PIK3CA_SCLC$PIK3CA.log2.TPM.1..Expression.22Q1.Public,
                       RIN3 = RIN3_SCLC$RIN3.log2.TPM.1..Expression.22Q1.Public,
                       SMAD4 = SMAD4_SCLC$SMAD4.log2.TPM.1..Expression.22Q1.Public,
                      GRIK3 = GRIK3_SCLC$GRIK3.log2.TPM.1..Expression.22Q1.Public,
                      CRMP1 = CRMP1_SCLC$CRMP1.log2.TPM.1..Expression.22Q1.Public,
                      MYT1L = MYT1L_SCLC$MYT1L.log2.TPM.1..Expression.22Q1.Public,
                      RALYL = RALYL_SCLC$RALYL.log2.TPM.1..Expression.22Q1.Public,
                      PDZRN3 = PDZRN3_SCLC$PDZRN3.log2.TPM.1..Expression.22Q1.Public,
                      MAP2 = MAP2_SCLC$MAP2.log2.TPM.1..Expression.22Q1.Public,
                      CACNA1E = CACNA1E_SCLC$CACNA1E.log2.TPM.1..Expression.22Q1.Public)
SCLC_df_new <- data.frame(ZFPM2 = ZFPM2_SCLC$ZFPM2.log2.TPM.1..Expression.22Q2.Public,
                          MAP7D3 = MAP7D3_SCLC$MAP7D3.log2.TPM.1..Expression.22Q2.Public)

NSCLC_df <- data.frame(EGFR = EGFR_NSCLC$EGFR.log2.TPM.1..Expression.22Q1.Public,
                      KIF19 = KIF19_NSCLC$KIF19.log2.TPM.1..Expression.22Q1.Public,
                      PAX6 = PAX6_NSCLC$PAX6.log2.TPM.1..Expression.22Q1.Public,
                      RIN3 = RIN3_NSCLC$RIN3.log2.TPM.1..Expression.22Q1.Public,
                      SMAD4 = SMAD4_NSCLC$SMAD4.log2.TPM.1..Expression.22Q1.Public,
                      GRIK3 = GRIK3_NSCLC$GRIK3.log2.TPM.1..Expression.22Q1.Public,
                      CRMP1 = CRMP1_NSCLC$CRMP1.log2.TPM.1..Expression.22Q1.Public,
                      MYT1L = MYT1L_NSCLC$MYT1L.log2.TPM.1..Expression.22Q1.Public,
                      RALYL = RALYL_NSCLC$RALYL.log2.TPM.1..Expression.22Q1.Public,
                      PDZRN3 = PDZRN3_NSCLC$PDZRN3.log2.TPM.1..Expression.22Q1.Public,
                      MAP2 = MAP2_NSCLC$MAP2.log2.TPM.1..Expression.22Q1.Public,
                      CACNA1E = CACNA1E_NSCLC$CACNA1E.log2.TPM.1..Expression.22Q1.Public)
NSCLC_df_new <- data.frame(ZFPM2 = ZFPM2_NSCLC$ZFPM2.log2.TPM.1..Expression.22Q2.Public,
                           MAP7D3 = MAP7D3_NSCLC$MAP7D3.log2.TPM.1..Expression.22Q2.Public)
ZFPM2_NSCLC
write_xlsx(Adeno_df, "Adeno Depmap expression.xlsx")
write_xlsx(Squamous_df, "Squamous Depmap expression.xlsx")
write_xlsx(SCLC_df, "SCLC Depmap expression.xlsx")
write_xlsx(NSCLC_df, "NSCLC Depmap expression.xlsx")
write_xlsx(SCLC_df_new, "SCLC new Depmap expression.xlsx")
write_xlsx(NSCLC_df_new, "NSCLC new Depmap expression.xlsx")


