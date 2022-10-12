setwd("C:/Users/Christoffer/OneDrive/1PhD/R files")
df <- read.table("Combat_filtered_exprs.txt", header = T)
info_df <- read.table("E-MTAB-6043.sdrf.txt", header = T)
#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6043/
df_AVENIO <- df[rownames(df) %in% grs$SYMBOL,]
df_AVENIO <- df_AVENIO %>% filter(rownames(df_AVENIO) %ni% bad)
avg <- rowMeans(df_AVENIO)
NSCLC_tot_enrichment
NSCLC_tot_enrichment_filt <- NSCLC_tot_enrichment %>% filter(genes %in% rownames(df_AVENIO))
NSCLC_tot_enrichment_filt
avg_df <- data.frame(genes = rownames(df_AVENIO), expression = avg)
avg_df <- avg_df[match(NSCLC_tot_enrichment_filt$genes,avg_df$genes),]
avg_df$enrichment <- NSCLC_tot_enrichment_filt$enrichment
res <- cor.test(avg_df$expression, avg_df$enrichment, method = "spearman")
p_values <- res$p.value
rhos <- res$estimate
ggplot(avg_df, aes(x = enrichment, y = expression))+
  geom_point(color="#ffa10c", size = 3)+
  theme_bw()+
  geom_smooth(aes(x=enrichment, y = expression), method = "lm", se = F, color = "black")+
  labs(title = paste("NSCLC RNA-seq correlation with cfChIP-seq"), 
       x = "cfChIP (Enrichment)", y = "Log2(TPM+1)",
        subtitle = paste("Spearman's rho =", round(rhos,3), 
        as.character(p_values), ", n =", length(avg_df$enrichment)))+
  th
avg_df$expression
avg_df

df_AVENIO
