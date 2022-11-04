library(karyoploteR)
library(pasillaBamSubset)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/BAM files and manifests")

genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)


####EGFR plot####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/EGFR track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
EGFR.region <- toGRanges("chr7:55,110,017-55,211,628")
kp <- plotKaryotype(zoom = EGFR.region, cex = 1)
kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
chip_bam <- c(A549_R1 = "A549_R1_ChIP.bam",
              A549_R2 = "A549_R2_ChIP.bam",
              A549_R3 = "A549_R3_ChIP.bam",
              HCC827_R1 = "HCC827_R1_ChIP.bam",
              HCC827_R2 = "HCC827_R2_ChIP.bam",
              HCC827_R3 = "HCC827_R3_ChIP.bam",
              HCC827_MET_R1 = "HCC827-MET_R1_ChIP.bam",
              HCC827_MET_R2 = "HCC827-MET_R2_ChIP.bam",
              HCC827_MET_R3 = "HCC827-MET_R3_ChIP.bam")
cell_color <- c(rep("firebrick",3),rep("#6a00fc",3),rep("#ffa10c",3))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(EGFR_maxes))
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    kpAxis(kp, ymin=0, ymax=max(EGFR_maxes), numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()



EGFR_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    EGFR_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
EGFR_maxes


####MET plot####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/MET track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
MET.region <- toGRanges("chr7:116,672,196-116,798,377")
kp <- plotKaryotype(zoom = MET.region, cex = 1)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
chip_bam <- c(A549_R1 = "A549_R1_ChIP.bam",
              A549_R2 = "A549_R2_ChIP.bam",
              A549_R3 = "A549_R3_ChIP.bam",
              HCC827_R1 = "HCC827_R1_ChIP.bam",
              HCC827_R2 = "HCC827_R2_ChIP.bam",
              HCC827_R3 = "HCC827_R3_ChIP.bam",
              HCC827_MET_R1 = "HCC827-MET_R1_ChIP.bam",
              HCC827_MET_R2 = "HCC827-MET_R2_ChIP.bam",
              HCC827_MET_R3 = "HCC827-MET_R3_ChIP.bam")
cell_color <- c(rep("firebrick",3),rep("#6a00fc",3),rep("#ffa10c",3))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(MET_maxes))
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    kpAxis(kp, ymin=0, ymax=max(MET_maxes), numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()

MET_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    MET_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
MET_maxes

####APC plot####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/APC track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
APC.region <- toGRanges("chr5:112,737,885-112,846,239")
kp <- plotKaryotype(zoom = APC.region, cex = 1)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
chip_bam <- c(A549_R1 = "A549_R1_ChIP.bam",
              A549_R2 = "A549_R2_ChIP.bam",
              A549_R3 = "A549_R3_ChIP.bam",
              HCC827_R1 = "HCC827_R1_ChIP.bam",
              HCC827_R2 = "HCC827_R2_ChIP.bam",
              HCC827_R3 = "HCC827_R3_ChIP.bam",
              HCC827_MET_R1 = "HCC827-MET_R1_ChIP.bam",
              HCC827_MET_R2 = "HCC827-MET_R2_ChIP.bam",
              HCC827_MET_R3 = "HCC827-MET_R3_ChIP.bam")
cell_color <- c(rep("firebrick",3),rep("#6a00fc",3),rep("#ffa10c",3))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(APC_maxes))
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    kpAxis(kp, ymin=0, ymax=max(APC_maxes), numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()

APC_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    APC_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
APC_maxes
TPM_AVENIO %>% filter(SYMBOL == "APC")


####EGFR-WT vs. EGFR-mut####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/EGFR cfChIP track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
EGFR.region <- toGRanges("chr7:55,110,017-55,211,628")
kp <- plotKaryotype(zoom = EGFR.region, cex = 1)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
setwd("D:/Lung cancer cfChIP/Deduped")
chip_bam <- c(NAC.1 = "A-1279-cfChIP.bam",
              NAC.4 = "D-1578-cfChIP.bam",
              NAC.2 = "B-1288-cfChIP.bam",
              NAC.3 = "C-1475-cfChIP.bam",
              NSC.1 = "E-439-cfChIP.bam",
              NSC.2 = "F-1449-cfChIP.bam",
              NSC.3 = "I-645-cfChIP.bam",
              NSC.4 = "J-1663-cfChIP.bam")
cell_color <- c(rep("#6a00fc",2),rep("#ffa10c",6))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(EGFR1_maxes))
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    kpAxis(kp, ymin=0, ymax=max(EGFR1_maxes), numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()



EGFR1_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    EGFR1_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
EGFR1_maxes


####KIF19####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/KIF19 cfChIP track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
KIF19.region <- toGRanges("chr17:74,340,226-74,355,820")
kp <- plotKaryotype(zoom = KIF19.region, cex = 1)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 5000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
setwd("D:/Lung cancer cfChIP/Deduped")
chip_bam <- c(NAC.1 = "A-1279-cfChIP.bam",
              NAC.2 = "B-1288-cfChIP.bam",
              NAC.3 = "C-1475-cfChIP.bam",
              NAC.4 = "D-1578-cfChIP.bam",
              NSC.1 = "E-439-cfChIP.bam",
              NSC.2 = "F-1449-cfChIP.bam",
              NSC.3 = "I-645-cfChIP.bam",
              NSC.4 = "J-1663-cfChIP.bam",
              SSC.1 = "G-514-cfChIP.bam",
              SSC.2 = "H-1169-cfChIP.bam",
              SSC.3 = "K-440-cfChIP.bam",
              SSC.4 = "L-1100-cfChIP.bam")
cell_color <- c(rep("#6a00fc",8),rep("#ffa10c",4))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMDensity(kp, data=bam.file, 
                           r0=at$r0, r1=at$r1, col = cell_color[i],
                           ymin = 0, ymax = max(KIF19_maxes),
                           window.size = width(KIF19.region)/100,
                           normalize = T)
    kpAxis(kp, data.panel = 1, tick.pos = c(0,1), 
           labels = c(0,format(round(max(KIF19_maxes), 3), scientific = T)),
           numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()



KIF19_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMDensity(kp, data=bam.file, 
                           r0=at$r0, r1=at$r1, col = cell_color[i],
                           window.size = width(KIF19.region)/100,
                           normalize = T)
    computed.ymax <- max(kp$latest.plot$computed.values$density)
    KIF19_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
KIF19_maxes


####RIN3####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/RIN3 cfChIP track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
RIN3.region <- toGRanges("chr14:92,683,781-92,688,994")
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kp <- plotKaryotype(zoom = RIN3.region, cex = 1)
kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
setwd("D:/Lung cancer cfChIP/Deduped")
chip_bam <- c(NAC.1 = "A-1279-cfChIP.bam",
              NAC.2 = "B-1288-cfChIP.bam",
              NAC.3 = "C-1475-cfChIP.bam",
              NAC.4 = "D-1578-cfChIP.bam",
              NSC.1 = "E-439-cfChIP.bam",
              NSC.2 = "F-1449-cfChIP.bam",
              NSC.3 = "I-645-cfChIP.bam",
              NSC.4 = "J-1663-cfChIP.bam",
              SSC.1 = "G-514-cfChIP.bam",
              SSC.2 = "H-1169-cfChIP.bam",
              SSC.3 = "K-440-cfChIP.bam",
              SSC.4 = "L-1100-cfChIP.bam")
cell_color <- c(rep("#6a00fc",8),rep("#ffa10c",4))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMDensity(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(RIN3_maxes),
                            window.size = width(RIN3.region)/100,
                            normalize = T)
    kpAxis(kp, data.panel = 1, tick.pos = c(0,1), 
           labels = c(0,format(round(max(RIN3_maxes), 3), scientific = T)),
           numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()



RIN3_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMDensity(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                           window.size = width(RIN3.region)/100,
                           normalize = T)
    computed.ymax <- max(kp$latest.plot$computed.values$density)
    RIN3_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
RIN3_maxes

####SLPI####
png(file = "C:/Users/Christoffer/OneDrive/1PhD/Manuskripter/Adeno, plano and SCLC article/For submission/Molecular oncology/Revised submission/figures and tables/supplementary/SLPI cfChIP track.png",
    width = 13200, height = 9725,
    units = "px", res = 1200)
SLPI.region <- toGRanges("chr20:45,252,239-45,254,564")
kp <- plotKaryotype(zoom = SLPI.region, cex = 1)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 5000, add.units = TRUE)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)
setwd("D:/Lung cancer cfChIP/Deduped")
chip_bam <- c(NAC.1 = "A-1279-cfChIP.bam",
              NAC.2 = "B-1288-cfChIP.bam",
              NAC.3 = "C-1475-cfChIP.bam",
              NAC.4 = "D-1578-cfChIP.bam",
              NSC.1 = "E-439-cfChIP.bam",
              NSC.2 = "F-1449-cfChIP.bam",
              NSC.3 = "I-645-cfChIP.bam",
              NSC.4 = "J-1663-cfChIP.bam",
              SSC.1 = "G-514-cfChIP.bam",
              SSC.2 = "H-1169-cfChIP.bam",
              SSC.3 = "K-440-cfChIP.bam",
              SSC.4 = "L-1100-cfChIP.bam")
cell_color <- c(rep("#6a00fc",8),rep("#ffa10c",4))

for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.2, r1=1, margin = 0.2)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i],
                            ymin = 0, ymax = max(SLPI_maxes))
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    kpAxis(kp, ymin=0, ymax=max(SLPI_maxes), numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}

dev.off()



SLPI_maxes <- c()
for(i in seq_len(length(chip_bam))) {
    bam.file <- paste0(chip_bam[i])
    at <- autotrack(i, length(chip_bam), r0=0.35, r1=1, margin = 0.3)
    kp <- kpPlotBAMCoverage(kp, data=bam.file, 
                            r0=at$r0, r1=at$r1, col = cell_color[i])
    computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
    SLPI_maxes[i] <- computed.ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1,
           cex = 0.4)
    kpAddLabels(kp, labels = names(chip_bam)[i], r0=at$r0, r1=at$r1, 
                cex=0.5, label.margin = 0.035)
}
SLPI_maxes

