setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/Data/Gene abundance files")
setwd("C:/Users/Christoffer/OneDrive/1PhD/R Kurser/R kursus")

####Load of BAM files####
q  #name of BAM file
z  #name of BAM index file
bam <- function(q,z){
  library(chromstaR)
  file <- readBamFileAsGRanges(
    q, bamindex = z,
    chromosomes = NULL,
    pairedEndReads = F,
    remove.duplicate.reads = T, 
    blacklist = NULL,
    what = "mapq")
  return(file)
}

####Raw counts#####
x #Name of BAM file
y #Granges object returned by gr() for the 197 AVENIO target genes
gene_count <- function(x,y){
  library(Rsamtools)
  bf <- BamFile(x)
  counts <- scanBam(bf, param = ScanBamParam(what = "pos", which = y))
  reads <- c()
  name <- c()
  for (i in 1:length(counts)){
    reads[i] <- length(counts[[i]]$pos)
    string <- names(counts)[i]
    splits <- strsplit(string, ":")
    chr <- splits[[1]][1]
    position <- as.numeric(strsplit(splits[[1]][2],"-")[[1]][1])
    g <- GRanges(chr, IRanges(start = position+10, end = position+110))
    gene <- y[y %over% g]
    name[i] <- gene$SYMBOL
  }
  df <- data.frame(genes = name, readcounts = reads)
  df <- df[order(df$genes),]
  return(df)
}


####Establishment of target region####
x #name on .txt file containing each gene, name of chromosome, start location, end location and strand
gr <- function(x){
  library(GenomicRanges)
  y <- read.table(x, header = T)
  y$start <- as.numeric(y$start) 
  y$end <- as.numeric(y$end) 
  if ("Regions" %in% colnames(y)){
    y$Regions <- as.numeric(y$Regions) 
  }
  granges <- GRanges(y$seqnames, 
                     IRanges(start = y$start, end = y$end),
                     strand = y$strand)
  if ("gene_start" %in% colnames(y)){
    y$gene_start <- as.numeric(y$gene_start)
    y$gene_end <- as.numeric(y$gene_end)
  }
  mcols(granges)$SYMBOL <- y$SYMBOL 
  return(granges)
}

####Load of RNA-seq data####
x #Name on combined TPM file
y #vector of names on samples
RNA_TPM <- function(x,y){
  library(dplyr)
  options(scipen=999)
  TPM <- read.table(x, header = T)
  df <- data.frame(SYMBOL = TPM$SYMBOL)
  for (i in 2:length(TPM)){
    res <- log2(TPM[,i]+1)
    df[ , ncol(df) + 1] <- res
    colnames(df)[ncol(df)] <- paste0(i)
  }
  setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
  AVENIO_gr <- gr("AVENIO_genes.txt")
  AVENIO_genelist <- AVENIO_gr$SYMBOL
  TPM_AVENIO <- df %>% filter(SYMBOL %in% AVENIO_genelist)
  colnames(TPM_AVENIO) <- c("SYMBOL",y)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
  return(TPM_AVENIO)
}
setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
TPM_AVENIO <- RNA_TPM("Gene abundance 25102022.txt",c("log2_A549_R1",
                                                      "log2_A549_R2",
                                                      "log2_A549_R3",
                                                      "log2_HCC827_R1",
                                                      "log2_HCC827_R2",
                                                      "log2_HCC827_R3",
                                                      "log2_HCC827-MET_R1",
                                                      "log2_HCC827-MET_R2",
                                                      "log2_HCC827-MET_R3"))
####mean TPM of cell lines####
x #data.frame returned by RNA_TPM()
y #colname of cell line to gen the mean log2(TPM+1)

RNA_mean <- function(x,y){
    df <- x[grepl(y,colnames(x))]
    return(data.frame(SYMBOL = x$SYMBOL,
                      log2_TPM = rowMeans(df)))
}


####Definition of the coverages in each gene of the AVENIO panel####
#KPRP and DLGAP2 have been put in manually. DLGAP2 is based on information given by Roche
#KPRP is made by manually looking on the HCC827 input file and determine where the coverage > 100
x #Name of bedgraph file used to study the coverages. Created using samtools. I use HCC827 input file
y #cutoff for minimum coverage value.
z #name on .txt file containing each gene, name of chromosome, start location, end location and strand

coverages <- function(x,y = 100,z) {
  library(GenomicRanges)
  library(dplyr)
  library(hwglabr2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  `%ni%` <- Negate(`%in%`)
  bg <- import_bedGraph(x)
  subsets <- bg[elementMetadata(bg)[,1]>y]
  reduced_subsets <- reduce(subsets)
  hg38_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  AVENIO <- read.table(z, header = T)
  AVENIO_gene_ids <- as.data.frame(org.Hs.egSYMBOL) %>% 
    filter(symbol %in% AVENIO$SYMBOL)
  AVENIO_gene_ranges <- hg38_genes[elementMetadata(hg38_genes)[,1] %in% 
                                     AVENIO_gene_ids$gene_id]
  symbols <- as.data.frame(org.Hs.egSYMBOL) %>% 
    filter(gene_id %in% AVENIO_gene_ranges$gene_id)
  symbols <- symbols[match(AVENIO_gene_ranges$gene_id, symbols$gene_id),]
  AVENIO_gene_ranges$symbols <- symbols$symbol
  tot <- c()
  for(i in 1:length(AVENIO_gene_ranges)){
    range <- AVENIO_gene_ranges[i]
    subjects <- subjectHits(findOverlaps(reduced_subsets,range))
    queris <- queryHits(findOverlaps(reduced_subsets,range))
    widt <- c()
    if (range$symbols != "DLGAP2"){
      for(j in 1:length(queris)){
        widt[j] <- width(reduced_subsets[queris[j]])
      }
      tot[i] <- sum(widt) 
    }
    else{
      tot[i] <- 1
    }
  }
  AVENIO_gene_ranges$coverage <- tot
  cover <- data.frame(SYMBOL = c(AVENIO_gene_ranges$symbols,"KPRP"), 
                      coverage = c(AVENIO_gene_ranges$coverage, 152761589-152759343))
  cover$coverage[187] = 875
  cover <- cover[order(cover$SYMBOL),]
  return(cover)
}

####Healthy normilization####
x #list of count tables returned by gene_count() from healhy controls
healthy <- function(x){
  library(GenomicRanges)
  library(GenomicAlignments)
  df <- data.frame(genes = x[[1]]$genes)
  for (i in 1:length(x)){
    df[ , ncol(df) + 1] <- x[[i]]$readcounts
    colnames(df)[ncol(df)] <- paste0(i)
  }
  av <- c()
  for (i in 1:length(df$genes)){
    av[i] <- mean(as.numeric(df[i,2:5]))
  }
  df$readcounts <- av
  return(df)
}


####Enrichment function####
x #Gene count table returned by gene_count()
y #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
z #Granges object returned by gr() for the 197 AVENIO target genes
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL
e.score <- function(x,y,z,h = NULL){
  library(GenomicAlignments)
  library(GenomicRanges)
  library(BiocParallel)
  library(dplyr)
  y <- read.table(y, header=T)
  y$coverage <- as.numeric(y$coverage)
  ChIP_reads <- x
  ChIP_reads <- ChIP_reads %>% filter(readcounts > 10)
  if(is.null(h)){
  }
  else {
    h <- h %>% filter(genes %in% ChIP_reads$genes)
    h <- h[match(ChIP_reads$genes, h$genes),]
    ChIP_reads$readcounts <- ChIP_reads$readcounts - h$readcounts 
  }
  y <- y[match(ChIP_reads$genes, y$SYMBOL),]
  len <- sum(ChIP_reads$readcounts)
  RPKM <- c()
  for (i in 1:length(ChIP_reads$genes)){
    RPKM[i] <- (ChIP_reads$readcounts[i]*1000*1000000)/(len*y$coverage[i])
  }
  ChIP_reads$e <- RPKM
  res <- data.frame(genes = ChIP_reads$genes, enrichment=ChIP_reads$e)
  res <- res[order(res$enrichment),]
  return(res)
}


####average enrichment between samples####
x #list of data.frames returned by e.score()
average_enrichment <- function(x){
  len <- length(x)
  std <- x[[1]]$genes
  df <- data.frame(genes = std)
  for (i in 1:len){
    x[[i]] <- x[[i]][match(std,x[[i]]$genes),]
    df[ , ncol(df) + 1] <- x[[i]]$enrichment
    colnames(df)[ncol(df)] <- paste0(i)
  }
  avg <- c()
  for (i in 1:length(df$genes)){
    avg[i] <- mean(as.numeric(df[i,2:length(df)]))
  }
  df$enrichment <- avg
  df <- df[order(df$enrichment),]
  return(df)
}



####bad genes####
x #name on .txt file containing each gene, the number of times each gene is sequenced and the distance each target is from TSS
y #cutoff of distance from TSS which is too long
badgene <- function(x,y = 25){
  library(dplyr)
  df <- read.table(x,header = T)
  df$Relative_distance <- as.numeric(df$Relative_distance)
  sr <- df %>% filter(Regions == 1) %>% filter(Relative_distance < y)
  gene_names <- sr$SYMBOL
  return(gene_names)
}

####Enrichment plot####
x #Ordered enrichment table returned by e.score()
y #Title of plot
b #optional: list of genes returned by badgene(), otherwise type NULL
gg_enrichment <- function(x, y, b = NULL){
  library(ggplot2)
  library(ggrepel)
  `%ni%` <- Negate(`%in%`)
  if (is.null(b)){
    x1 <- x
    t <- ""
  }
  else{
    x1 <- x %>% filter(genes %ni% b)
    t <- "sorted "
  }  
  x1$reverse <- c(1:length(x1$genes))
  gg <- ggplot(data = x1, aes(x=reverse, y=enrichment)) +
    geom_point(aes(colour=log2(enrichment)), size = 3)+
    scale_colour_gradient(low ="#ffa10c", high = "#6a00fc",
                          name = "Log2(Enrichment)", limits = c(min(log2(x1$enrichment)), max(log2(x1$enrichment))))+
    labs(title = paste(y,t), x = "Genes", y = "Enrichment")+
    scale_y_continuous(trans='log2')+
    geom_label_repel(
      aes(label = ifelse(reverse > (length(x1$reverse))-10, as.character(genes), "")),
      box.padding = 0.35,
      point.padding = 0.1,
      label.padding = 0.2,
      nudge_x = -30,
      nudge_y = max(log2(x1$enrichment))/20,
      segment.color = "black",
      label.size = 0.5,
      parse = F)+
    theme_bw()+
    th
  return(gg)
}


####Normalized correlation plot####
x #Enrichment data.frame return from e.score() for sample 1
y #Enrichment data.frame return from e.score() for sample 2
p #Name of sample 1
q #Name of sample 2
b #optional: list of genes returned by badgene(), otherwise type NULL
z #optional: cut-off of difference, otherwise NULL for top 15 differentially expressed genes

ChIPcorr <- function(x,y,p,q,b = NULL,z = NULL){
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
        aes(label=ifelse(reverse> (length(reverse))-15, as.character(paste(genes)),"")),
        segment.color="#6a00fc",
        color="#6a00fc",
        label.size = NA,
        nudge_y = 0,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      geom_label_repel(
        aes(label=ifelse(reverse < 16, as.character(paste(genes)),"")),
        segment.color="#ffa10c",
        color="#ffa10c",
        label.size = NA,
        nudge_y = 0,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      theme_bw()+
      geom_label(x = 40, y = (max(df$log2FC)-0.5), 
                 label = paste("Upregulated in",p),
                 color = "#6a00fc", label.size = 0, size = 5)+
      geom_label(x = 120, y = (min(df$log2FC)+0.5), label = paste("Upregulated in",q),
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
        aes(label=ifelse(log2FC>z, as.character(paste(genes)),"")),
        segment.color="#6a00fc",
        color="#6a00fc",
        label.size = NA,
        nudge_y = 0,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      geom_label_repel(
        aes(label=ifelse(log2FC < (-z), as.character(paste(genes)),"")),
        segment.color="#ffa10c",
        color="#ffa10c",
        label.size = NA,
        nudge_y = 0,
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

ChIPcorr(Adeno1_enrichment, Adeno_tot_not1_enrichment, "adenocarcinoma 1", "average adenocarcinoma", bad)
ChIPcorr(NSCLC_tot_enrichment, SCLC_tot_enrichment, "average NSCLC", "average SCLC", bad)

####Normalized correlation plot for cells####

x #Enrichment data.frame return from e.score() for sample 1
y #Enrichment data.frame return from e.score() for sample 2
p #x-axis title
q #y-axis title
b #optional: list of genes returned by badgene(), otherwise type NULL
z #cut-off of difference, otherwise NULL for top 15 differentially expressed genes
ChIPcorr_cell <- function(x,y,p,q,b = NULL,z = NULL){
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
  df <- data.frame(genes = x$genes, x = log2(x$enrichment), y = log2(y$enrichment))
  if (is.null(b)){
    df <- df
  }
  else{
    df <- df %>% filter(genes %ni% b)
  } 
  res <- summary(lm(y~x, data=df))
  a <- res$coefficients[2,1]
  c <- res$coefficients[1,1]
  df$norm <- df$y-((df$x*a)+c)
  df <- df[order(df$norm),]
  df$reverse <- c(1:length(df$genes))
  if(!is.null(z)){
    gg <- ggplot(data = df, aes(x = x, y = y))+
      geom_point(aes(colour=norm), size = 4)+
      scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="grey", 
                             high="#6a00fc", name = expression(log[2]("Difference")), 
                             limits=c(round(min(df$norm))-1,round(max(df$norm))+1))+
      labs(title=paste("Average ChIP-seq enrichment of", q, "and", p), 
           x = p, y=q)+
      theme_bw()+
      geom_smooth(method="lm", se=F, color = "black")+
      geom_abline(intercept = c+z, slope = a, color = "Black", linetype = "dashed", size = 1)+
      geom_abline(intercept = c-z, slope = a, color = "Black", linetype = "dashed", size = 1)+
      geom_label_repel(
        aes(label=ifelse(norm>z, as.character(paste(genes)),"")),
        segment.color="#6a00fc",
        color = "#6a00fc",
        label.size = 0.5,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      geom_label_repel(
        aes(label=ifelse(norm < (-z), as.character(paste(genes)),"")),
        segment.color="#ffa10c",
        color = "#ffa10c",
        label.size = 0.5,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      th
  }
  else{
    gg <- ggplot(data = df, aes(x = x, y = y))+
      geom_point(aes(colour=norm), size = 4)+
      scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="grey", 
                             high="#6a00fc", name = expression(log[2]("difference")), 
                             limits=c(round(min(df$norm))-1,round(max(df$norm))+1))+
      labs(title=paste("Average ChIP-seq enrichment of", q,"and", p), 
           x = p, y=q)+
      theme_bw()+
      geom_smooth(method="lm", se=F, color = "black")+
      geom_label_repel(
        aes(label=ifelse(reverse > (length(reverse))-15, as.character(paste(genes)),"")),
        segment.color="#6a00fc",
        color = "#6a00fc",
        label.size = NA,
        nudge_y = 0,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      geom_label_repel(
        aes(label=ifelse(reverse < 16, as.character(paste(genes)),"")),
        segment.color="#ffa10c",
        color = "#ffa10c",
        label.size = NA,
        nudge_y = 0,
        parse = F,
        size = 3.5,
        max.overlaps = 100,
        fill = alpha(c("white"),0))+
      th
  }
  return(gg)
}


####Plot of enrichment relative to TSS####
x #File with all Targets in AVENIO panel. Containing chromosome, start of target, end of target, gene SYMBOL
# Gene start, Gene end, strand, mean position of target and relative distance of target from TSS
y #List of gene symbols for active genes
z #List of gene symbols for inactive genes
b #name of Granges object returned by bam()
t #Name of cell line
e.score_distribution <- function(x,y,z,b, t){
  library(GenomicAlignments)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  targets <- gr(x)
  targets <- targets[targets$SYMBOL %in% c(y,z)]
  x1 <- read.table(x, header = T)
  xx <- x1[(x1$SYMBOL %in% c(y,z)),]
  targets$length <- xx$end - xx$start
  se_ChIP_R1 <- summarizeOverlaps(features = targets, reads = b$R1, 
                                  ignore.strand = T, mode = "Union")
  se_ChIP_R2 <- summarizeOverlaps(features = targets, reads = b$R2, 
                                  ignore.strand = T, mode = "Union")
  se_ChIP_R3 <- summarizeOverlaps(features = targets, reads = b$R3, 
                                  ignore.strand = T, mode = "Union")
  ChIP_reads <- data.frame(genes = targets$SYMBOL, 
                           readcounts = assays(se_ChIP_R1)$counts[1:length(assays(se_ChIP_R1)$counts)]+
                               assays(se_ChIP_R2)$counts[1:length(assays(se_ChIP_R2)$counts)]+
                               assays(se_ChIP_R3)$counts[1:length(assays(se_ChIP_R3)$counts)])
  ChIP_reads$coverage <- targets$length
  ChIP_reads <- ChIP_reads %>% filter(readcounts > 10)
  len <- as.numeric(sum(ChIP_reads$readcounts))
  RPKM <- c()
  for (i in 1:length(ChIP_reads$genes)){
    RPKM[i] <- (ChIP_reads$readcounts[i]*1000*1000000)/(len*ChIP_reads$coverage[i])
  }
  ChIP_reads$enrichment <- RPKM
  ChIP_reads$activity <- c(1:length(ChIP_reads$enrichment))
  for (i in 1:length(unique(ChIP_reads$genes))){
    norms <- c()
    c <- ChIP_reads %>% filter(genes %in% unique(ChIP_reads$genes)[i])
    avg <- median(c$enrichment)
    if (unique(ChIP_reads$genes)[i] %in% y) {
      activity <- "Active"
    }
    else {
      activity <- "Inactive"
    }
    ChIP_reads$activity[ChIP_reads$genes %in% unique(ChIP_reads$genes)[i]] <- activity
  }
  ChIP_reads$relative <- xx$Relative_distance
  active <- ChIP_reads %>% filter(activity == "Active")
  ac <- paste(c(sort(unique(active$genes))), collapse=', ')
  inactive <- ChIP_reads %>% filter(activity == "Inactive")
  ina <- paste(c(sort(unique(inactive$genes))), collapse=', ')
  cols <- c("green4", "firebrick")
  names(cols) <- c(paste("Active genes \nn = ", length(unique(active$genes)), "\n"),
                   paste("Inactive genes \nn = ", length(unique(inactive$genes)),"\n"))
  gg <- ggplot(ChIP_reads, aes(x = relative, y = enrichment),group = activity, fill = activity)+
    geom_smooth(data = ChIP_reads[ChIP_reads$activity=="Active",],
                method = "loess", span = 1.3, se = F,
                aes(colour = names(cols)[1], fill = names(cols)[1]))+
    geom_smooth(data = ChIP_reads[ChIP_reads$activity=="Inactive",],
                method = "loess", span = 1.3, se = F,
                aes(colour = names(cols)[2], fill = names(cols)[2]))+
    scale_color_manual(name = "Gene activity", 
                       values = c("green4", "firebrick"),
                       aesthetics = c("colour", "fill"))+
    theme_bw(base_size = 17)+
    
    geom_ribbon(data = ChIP_reads[ChIP_reads$activity=="Active",], 
                aes(ymin = 0,
                    ymax = predict(loess(enrichment ~ relative, 
                                         span=1.3))),
                alpha = 0.3,fill = "green4")+
    geom_ribbon(data = ChIP_reads[ChIP_reads$activity=="Inactive",], 
                aes(ymin = 0, 
                    ymax = predict(loess(enrichment ~ relative, 
                                         span=1.3))),
                alpha = 0.5,fill = "firebrick")+
    scale_x_continuous(limits=c(0,100))+
    labs(title =paste(t),
         x = "% Relative to TSS", y="H3K36me3 Enrichment")+
    th
  return(gg)
}
?summarizeOverlaps

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
####Comparison of ChIP-seq with RNA-seq on cell lines####
x #data.frame of enrichment in ChIP sample returned by e.score()
y #data.frame of mean log2(TPM+1) values from RNA-seq of the 197 AVENIO genes returned by ran_mean()
z #name of cell
r #Cutoff for active gene determined by RNA-seq
b #optional: list of genes returned by badgene(), otherwise type NULL
a #optional: TRUE if only active genes based on RNA should be analyzed, otherwise type FALSE
g #optional: Genes to be excluded due to being an outlier (otherwise NULL)
versus <- function(x,y,z,r,b = NULL,a = F,g = NULL){
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
  `%ni%` <- Negate(`%in%`)
  y1 <- y
  y1 <- y1 %>% filter(SYMBOL %in% x$genes)
  y1 <- y1[(match(x$genes, y1$SYMBOL)),]
  if (is.null(b)){
    y1 <- y1
    x1 <- x
    t <- ""
  }
  else{
    y1 <- y1 %>% filter(SYMBOL %ni% b)
    x1 <- x %>% filter(genes %ni% b)
    t <- ""
  }  
  if (is.null(g)){
    y1 <- y1
    x1 <- x1
  }  
  else{
    y1 <- y1 %>% filter(SYMBOL %ni% g)
    x1 <- x1 %>% filter(genes %ni% g)
  }
  if (a == TRUE){
    y2 <- y1 %>% filter(y1[,2] > r)
    x2 <- x1[match(y2$SYMBOL,x1$genes),]
    y2 <- y2[(match(x2$genes, y2$SYMBOL)),]
    df <- data.frame(genes = x2$genes, ChIP = x2$enrichment, RNA = y2[,2])
    tit <-"for active genes"
  }
  else{
    y1 <- y1 %>% filter(SYMBOL %in% x1$genes)
    y1 <- y1[(match(x1$genes, y1$SYMBOL)),]
    df <- data.frame(genes = x1$genes, ChIP = x1$enrichment, RNA = y1[,2])
    tit <- "for all genes"
  }
  res <- cor.test(df$RNA, df$ChIP, method = "spearman")
  p_values <- res$p.value
  if (p_values < 0.0001){
    p <- ", P < 0.0001"
  }
  else{
    p <- paste(", P =", p_values)
  }
  rhos <- as.numeric(res$estimate)
  gg <- ggplot(df, aes(x = ChIP, y = RNA))+
    geom_point(color="#ffa10c", size = 3)+
    theme_bw()+
    geom_smooth(aes(x=ChIP, y = RNA), method = "lm", se = F, color = "black")+
    labs(title = paste(z,"RNA-seq correlation with ChIP-seq",tit), 
         x = "Average ChIP (Enrichment)", y = expression(bold(paste("Average ",log[2]("TPM+1")))),
         subtitle = paste("Spearman's rho =", round(rhos,3), 
                          as.character(p), ", n =", length(df$ChIP)))+
    th
  return(gg)
}


#labs(title = paste(cell, " RNA-seq correlation with ChIP-seq for ",tit,t,"genes", sep = ""), 
#x = "ChIP (Enrichment)", y = "RNA Log2(TPM+1)",
#subtitle = paste("Spearman's rho =", round(rhos,3), as.character(p), ", n =", length(df$ChIP)))



####ROC data####

x #data.frame of enrichment in ChIP sample returned by e.score()
y #data.frame of log2(TPM+1) values from RNA-seq of the 197 AVENIO genes returned by mean_RNA()
r #Cutoff for active gene determined by RNA-seq
b #optional: list of genes returned by badgene(), otherwise type NULL

get_ROC_data <- function(x,y,r,b = NULL){
  library(dplyr)
  library(stringr)
  `%ni%` <- Negate(`%in%`)
  y1 <- y
  y1 <- y1 %>% filter(SYMBOL %in% x$genes)
  y1 <- y1[(match(x$genes, y1$SYMBOL)),]
  if (is.null(b)){
    y1 <- y1
    x1 <- x
  }
  else{
    y1 <- y1 %>% filter(SYMBOL %ni% b)
    x1 <- x %>% filter(genes %ni% b)
  }
  RNA_inactive <- y1 %>% filter(y1[,2] <= r)
  RNA_active <- y1 %>% filter(y1[,2] > r)
  ChIP_inactive <- x %>% filter(genes %in% RNA_inactive$SYMBOL)
  ChIP_active <- x %>% filter(genes %in% RNA_active$SYMBOL)
  return(list(active = ChIP_active,
              inactive = ChIP_inactive))
}

####RNA FC plot####
x #Name on column in TPM df
y #Name on column in TPM df
z #Name on TPM df
d #tick distance on Y-axis
b #optional: list of genes returned by badgene(), otherwise type NULL
FCplot <- function(x,y,z,b = NULL){
      library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
  `%ni%` <- Negate(`%in%`)
  z1 <- z[colnames(z) %in% c("SYMBOL",x,y)]
  index1 <- c(1:length(z1))[colnames(z1) %in% x]
  index2 <- c(1:length(z1))[colnames(z1) %in% y]
  z1$FC <- z1[,index1]-z1[,index2]
  xt <- gsub("log2.", '',x)
  yt <- gsub("log2.", '',y)
  if (is.null(b)){
    z1 <- z1
  }
  else{
    z1 <- z1 %>% filter(SYMBOL %ni% b)
  } 
  gg <- ggplot(data=z1, aes(x=c(1:length(FC)), y=FC))+
    geom_point(aes(colour = FC), size = 3)+
    scale_colour_gradient2(midpoint = 0, low="#ffa10c", mid ="#808080", 
                           high="#6a00fc", name = "Log2FC", 
                           limits=c((min(z1$FC)-0.1),(max(z1$FC)+0.1)))+
    theme_bw()+
    geom_hline(yintercept = 0, color = "Black", size = 1)+
    geom_hline(yintercept = 1, color = "Black", linetype = "dashed", size = 1)+
    geom_hline(yintercept = -1, color = "Black", linetype = "dashed", size = 1)+
    labs(title=paste(xt,"compared to",yt), 
         x = "Genes", y="Log2FC")+
    xlim(0, length(z1$FC))+
    ylim(round(min(z1$FC))-1.5, round(max(z1$FC))+1.5)+
    geom_label(x = 30, y = round(max(z1$FC))+1.5, label = paste("Upregulated in",xt),
               color = "#6a00fc", label.size = 0, size = 5)+
    geom_label(x = 30, y = round(min(z1$FC))-1.5, label = paste("Upregulated in",yt),
               color = "#ffa10c", label.size = 0, size = 5)+
    geom_label_repel(
      aes(label=ifelse(FC < -1,
                       as.character(SYMBOL),"")),
      segment.color="#ffa10c",
      color="#ffa10c",
      label.size = NA,
      parse = F,
      size = 3.5,
      max.overlaps = 100,
      fill = alpha(c("white"),0))+
    geom_label_repel(
      aes(label=ifelse(FC > 1,
                       as.character(SYMBOL),"")),
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
  return(gg)
}


####Histogram of coverage####
#.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
cover_hist <- function(x){
  library(dplyr)
  library(ggplot2)
  df <- read.table(x,header = T)
  gg <- ggplot(data = df)+
    geom_histogram(aes(coverage), binwidth = 500, bins = 20, fill="#69b3a2", color="#e9ecef")+
    theme_bw()+
    labs(title = "NGS coverage in each gene", x = "Coverage length (bp)", y="Frequency")+
    th
  return(gg)
}

####histogram of regions sequenced relative to TSS####
x #File with all Targets in AVENIO panel. Containing chromosome, start of target, end of target, gene SYMBOL
# Gene start, Gene end, strand, mean position of target and relative distance of target from TSS
TSS_hist <- function(x){
  library(dplyr)
  library(ggplot2)
  df <- read.table(x,header = T)
  sr <- df %>% filter(Regions == 1)
  gg <- ggplot(data = sr)+
    geom_histogram(aes(Relative_distance), binwidth = 5, bins = 20, fill="#69b3a2", color="#e9ecef")+
    theme_bw()+
    labs(title = "Distance from TSS genes are captured", x = "Relative distance from TSS (%)", y="Frequency")+
    ylim(c(0,20))+
    th
  return(gg)
}

####histogram of number of regions sequenced####
x #File with all Targets in AVENIO panel. Containing chromosome, start of target, end of target, gene SYMBOL
# Gene start, Gene end, strand, mean position of target and relative distance of target from TSS
number_hist <- function(x){
  library(dplyr)
  library(ggplot2)
  df <- read.table(x,header = T)
  df <- df[!duplicated(df$SYMBOL),]
  gg <- ggplot(data = df)+
    geom_histogram(aes(Regions), binwidth = 1, bins = 20, fill="#69b3a2", color="#e9ecef")+
    theme_bw()+
    labs(title = "Capture regions in each gene", x = "Number of capture regions", y="Frequency")+
    th
  return(gg)
}


####Transpose df####
x #data.frame returned by e.score()
y #name of sample
transpose <- function(x,y){
  x1 <- as.data.frame(t(x))
  colnames(x1) <- x1[1,1:length(x1)]
  x1 <- x1[-1,]
  rownames(x1) <- y
  return(x1)
}
####combine transposed data.frames####
x #list of data.frames returned by e.score()
y #vector of names for each data.frame.
z #vector of names for groups
bind_samples <- function(x,y,z){
  library(dplyr)
  len <- length(x)
  lens <- c()
  for (i in 1:len){
    tran <- transpose(x[[i]], y[i])
    x[[i]] <- tran
    lens[i] <- length(tran)
  }
  short <- which.min(lens)
  std <- colnames(x[[short]])
  for (i in 1:len){
    x[[i]] <- x[[i]] %>% dplyr::select(all_of(std))
  }
  res <- bind_rows(x)
  res <- as.data.frame(sapply(res, as.numeric))
  rownames(res) <- y
  res$groups <- z
  return(res)
}

####PCA####
x #data.frame containing rows of cfChIP-seq enrichment samples. columns contain gene names. Returned by bind_samples()
y #Title of plot
PCA <- function(x,y){
  library(ggbiplot)
  group <- x$groups
  x <- x %>% select(-bad)
  x.pca <- prcomp(x[1:(length(x)-1)],center = T, scale = T)
  gg <- ggbiplot(x.pca, ellipse = T, obs.scale = 1, var.scale = 1,var.axes=FALSE,
                 labels=rownames(x), groups=group, labels.size = 4,
                 varname.size = 5)+
    scale_colour_manual(name="Subtype", values= c("forest green", "red3", "dark blue", "#ffa10c"))+
    ggtitle(y)+
    theme_bw()+
    theme(legend.position = "bottom")+
    th
  return(gg)
}

#https://www.datacamp.com/community/tutorials/pca-analysis-r


####tSNE plot####
x #data.frame containing rows of cfChIP-seq enrichment samples. columns contain gene names. Returned by bind_samples()
z #name of plot
tSNE <- function(x,z){
  library(Rtsne)
  library(dplyr)
  library(ggrepel)
  set.seed(142)
  tSNE_fit <- x %>% select(-bad) %>% 
    select(where(is.numeric)) %>%
    Rtsne(perplexity = floor((nrow(x) - 1) / 3),
          dims = 2)
  tSNE_df <- tSNE_fit$Y %>%
    as.data.frame()
  colnames(tSNE_df) <- c("tSNE1", "tSNE2")
  rownames(tSNE_df) <- rownames(x)
  tSNE_df$group <- x$groups
  gg <- ggplot(data=tSNE_df, aes(x = tSNE1, y = tSNE2, color = group))+
    geom_point(size = 5)+
    labs(title = z)+
    theme_bw()+
    scale_color_manual(name = "Tumor histology", 
                       values = c("#ffa10c", "red3","Darkolivegreen", "#6a00fc"),
                       aesthetics = c("colour", "fill"))+
    geom_label_repel(
      aes(label=ifelse(group == "Adeno",
                       as.character(rownames(x)),"")),
      box.padding   = 0.35,
      point.padding = 0.1,
      label.padding = 0.2,
      nudge_x = 0,
      nudge_y=2,
      segment.color="#ffa10c",
      color="#ffa10c",
      label.size = 0.7,
      parse = F,
      size = 4.5,
      max.overlaps = 100)+
    geom_label_repel(
      aes(label=ifelse(group == "SCLC",
                       as.character(rowname,s(x)),"")),
      nudge_x = 0,
      nudge_y=2,
      box.padding   = 0.35,
      point.padding = 0.1,
      label.padding = 0.2,
      segment.color="Darkolivegreen",
      color="Darkolivegreen",
      label.size = 0.7,
      parse = F,
      size = 4.5,
      max.overlaps = 100)+
    geom_label_repel(
      aes(label=ifelse(group == "Squamous",
                       as.character(rownames(x)),"")),
      box.padding   = 0.35,
      point.padding = 0.1,
      label.padding = 0.2,
      nudge_x = 0,
      nudge_y=2,
      segment.color="#6a00fc",
      color="#6a00fc",
      label.size = 0.7,
      parse = F,
      size = 4.5,
      max.overlaps = 100)+
    geom_label_repel(
      aes(label=ifelse(group == "Healthy",
                       as.character(rownames(x)),"")),
      box.padding   = 0.35,
      point.padding = 0.1,
      label.padding = 0.2,
      nudge_x = 0,
      nudge_y=2,
      segment.color="red3",
      color="red3",
      label.size = 0.7,
      parse = F,
      size = 4.5,
      max.overlaps = 100)+
    th
  return(gg)
}

####Venn diagram of ChIP compared to RNA####
x #Enrichment data.frame return from e.score() for sample 1
y #Enrichment data.frame return from e.score() for sample 2
z #Name on TPM df
p #Name on column in TPM df of sample 2
q #Name on column in TPM df of sample 1
b #optional: list of genes returned by badgene(). Default is NULL
d #Cutoff of difference

venn <- function(x,y,z,p,q, b = NULL, d){
  library(dplyr)
  library(ggVennDiagram)
  library(ggpubr)
  
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
  df <- data.frame(genes = x$genes, x = log2(x$enrichment), y = log2(y$enrichment))
  if (is.null(b)){
    df <- df
    t <- ""
  }
  else{
    df <- df %>% filter(genes %ni% b)
    t <- "sorted "
  } 
  res <- summary(lm(y~x, data=df))
  a <- res$coefficients[2,1]
  c <- res$coefficients[1,1]
  df$norm <- df$y-((df$x*a)+c)
  df <- df[order(df$norm),]
  
  z1 <- z[colnames(z) %in% c("SYMBOL",p,q)]
  index1 <- c(1:length(z1))[colnames(z1) %in% p]
  index2 <- c(1:length(z1))[colnames(z1) %in% q]
  z1$FC <- z1[,index1]-z1[,index2]
  xt <- gsub("log2.", '',p)
  yt <- gsub("log2.", '',q)
  if (is.null(b)){
    z1 <- z1
    t <- ""
  }
  else{
    z1 <- z1 %>% filter(SYMBOL %ni% b)
    t <- "sorted "
  } 
  rna_df_high <- z1 %>% filter(FC > d)
  ChIP_df_high <- df %>% filter(norm > 0) %>% filter(genes %in% rna_df_high$SYMBOL)
  rna_df_low <- z1 %>% filter(FC < -d)
  ChIP_df_low <- df %>% filter(norm < 0) %>% filter(genes %in% rna_df_low$SYMBOL)
  high_genes <- list(rna_df_high$SYMBOL, ChIP_df_high$genes)
  low_genes <- list(rna_df_low$SYMBOL, ChIP_df_low$genes)
  gg_high <- ggVennDiagram(high_genes, set_color = "black", 
                           edge_lty = 1, edge_size = 1, 
                           label_percent_digit = 2,
                           category.names = c("RNA expression", "ChIP enrichment"),
                           color = "black", label_size = 5, label = "count")+
    scale_fill_gradient(low = "Yellow3", high = "Darkcyan", name = "Genes")+
    labs(title = paste("Genes upregulated in", xt, 
                       "compared to",yt,"\ndetermined by RNA-seq or ChIP-seq"), x="", y="")+
    theme_void(base_size = 12)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank())
  gg_low <- ggVennDiagram(low_genes, set_color = "black", 
                          edge_lty = 1, edge_size = 1, 
                          label_percent_digit = 2,
                          category.names = c("RNA expression", "ChIP enrichment"),
                          color = "black", label_size = 5, label = "count")+
    scale_fill_gradient(low = "Yellow3", high = "Darkcyan", name = "Genes")+
    labs(title = paste("Genes upregulated in", yt, 
                       "compared to",xt,"\ndetermined by RNA-seq or ChIP-seq"), x="", y="")+
    theme_void(base_size = 12)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank())
  gg <- ggarrange(gg_high, gg_low + rremove("x.text"),
                  ncol = 1, nrow = 2)
  return(gg)
}

venn(A549_enrichment, HCC827_enrichment, TPM_AVENIO, "log2.HCC827", "log2.A549",
     bad, d = 1)
venn(clone3_enrichment, HCC827_enrichment, TPM_AVENIO, "log2.HCC827", "log2.Clone3",
     bad, d = 1)


####Bar graph of ChIP compared to RNA####
x #Enrichment data.frame return from e.score() for sample 1
y #Enrichment data.frame return from e.score() for sample 2
p #Name on column in TPM df of sample 2
q #Name on column in TPM df of sample 1
b #optional: list of genes returned by badgene(). Default is NULL
d #Cutoff of difference

bar_overlap <- function(x,y,p,q, b = NULL, d){
  or_wd <- getwd()
  options(scipen=999)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  `%ni%` <- Negate(`%in%`)
  if(length(x$genes) == length(y$genes)){
    y = y
    x = x
    nix = 0
    niy = 0
  }
  else{
    if (length(x$genes) > length(y$genes)){
      niy <- x %>% filter(genes %ni% y$genes)
      niy <- length(niy$genes)
      x <- x %>% filter(genes %in% y$genes)
    }
    else{
      niy = 0
    }
    if (length(y$genes) > length(x$genes)){
      nix <- y %>% filter(genes %ni% x$genes)
      nix <- length(nix$genes)
      y <- y %>% filter(genes %in% x$genes)
    }
    else{
      nix = 0
    }
  }
  y <- y[(match(x$genes, y$genes)),]
  df <- data.frame(genes = x$genes, x = log2(x$enrichment), y = log2(y$enrichment))
  if (is.null(b)){
    df <- df
    t <- ""
  }
  else{
    df <- df %>% filter(genes %ni% b)
    t <- "sorted "
  } 
  res <- summary(lm(y~x, data=df))
  a <- res$coefficients[2,1]
  c <- res$coefficients[1,1]
  df$norm <- df$y-((df$x*a)+c)
  df <- df[order(df$norm),]
  setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
  TPM <- read.table("Gene abundance 25102022.txt", header = T)
  TPM <- TPM[TPM$SYMBOL %in% df$genes,]
  colnames(TPM)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
  z1 <- dif_gene_express(TPM,c(p,q))
  z1$SYMBOL <- z1$gene
  z1$FC <- z1$Log2FC
  setwd(or_wd)
  xt <- gsub("log2.", '',p)
  yt <- gsub("log2.", '',q)
  if (is.null(b)){
    z1 <- z1
    t <- ""
  }
  else{
    z1 <- z1 %>% filter(SYMBOL %ni% b)
    t <- "sorted "
  }
  rna_df_high <- z1 %>% filter(FC > d)
  ChIP_df_high <- df %>% filter(norm > 0) %>% filter(genes %in% rna_df_high$SYMBOL)
  rna_df_low <- z1 %>% filter(FC < -d)
  ChIP_df_low <- df %>% filter(norm < 0) %>% filter(genes %in% rna_df_low$SYMBOL)
  ChIP_negative_high <- rna_df_high %>% filter(SYMBOL %ni% ChIP_df_high$genes)
  ChIP_negative_low <- rna_df_low %>% filter(SYMBOL %ni% ChIP_df_low$genes)
  df1 <- data.frame(discovery = factor(c("Agree", "Disagree",
                                  "Agree", "Disagree")),
                    activity = c(rep(paste(xt,"upregulated"),2),
                                 rep(paste(yt,"upregulated"),2)),
                    numbergenes = c(length(ChIP_df_high$genes)+nix,
                                    length(ChIP_negative_high$SYMBOL)-nix,
                                    length(ChIP_df_low$genes)+niy,
                                    length(ChIP_negative_low$SYMBOL)-niy))
  percent1 <- (df1$numbergenes[1])/(df1$numbergenes[1]+df1$numbergenes[2])*100
  percent2 <- (df1$numbergenes[3])/(df1$numbergenes[3]+df1$numbergenes[4])*100
  height1 <- length(ChIP_df_high$genes)+length(ChIP_negative_high$SYMBOL)
  height2 <- length(ChIP_df_low$genes)+length(ChIP_negative_low$SYMBOL)
  scales <- c(height1/5, height2/5)
  height1 <- height1 + min(scales)
  height2 <- height2 + min(scales)
  df1$discovery2 <- relevel(df1$discovery, 'Disagree')
  gg <- ggplot(data=df1, aes(x = fct_inorder(activity), y = numbergenes, fill = discovery2))+
    geom_bar(stat = "identity")+
    theme_bw()+
    labs(x = "", y="Number of genes")+
    scale_fill_manual("ChIP results", values = c("#6a00fc","#ffa10c"))+
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

TPM_AVENIO

####Sample metrics####
x #Name of BAM file
y #name on .txt file containing each gene, name of chromosome, start location, end location and strand
c #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
b #BAM file folder
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL

metrics <- function(x,y,c,b,h = NULL){
  setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
  targets <- gr(y)
  setwd(b)
  readcounts <- gene_count(x,targets)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
  enrichment <- e.score(readcounts, c, targets,h)
  len <- sum(readcounts$readcounts)
  enrichment <- median(enrichment$enrichment)
  res <- list(total_reads = len, 
              median_enrichment = enrichment)
  return(res)
}


####Fragment lengths####
x #Vector of names for BAM files for group 1
y #Vector of names for BAM files for group 2
z #vector of names for group 1 and 2
p #Title of plot
fragment_length <- function(x,y,z,p){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  what <- c("isize")
  param <- ScanBamParam(which = grs, what = what)
  bam1 <- list()
  bam2 <- list()
  for (i in 1:length(x)){
    bam <- scanBam(x[i], param = param)
    bam_lengths <- unname(unlist(bam))
    bam_lengths <- bam_lengths[bam_lengths>0]
    bam1[[i]] <- bam_lengths
  }
  for (i in 1:length(y)){
    bam <- scanBam(y[i], param = param)
    bam_lengths <- unname(unlist(bam))
    bam_lengths <- bam_lengths[bam_lengths>0]
    bam2[[i]] <- bam_lengths
  }
  bam1 <- unlist(bam1)
  bam2 <- unlist(bam2)
  df <- data.frame(len = c(bam1,bam2),
                   Sample = c(rep(z[1],length(bam1)),rep(z[2],length(bam2))))
  mu <- ddply(df, "Sample", summarise, grp.median=median(len))
  gg <- ggplot(data = df, aes(x = len, color = Sample))+
    geom_density(size = 1)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
               linetype="dashed", size = 1)+
    theme_bw(base_size = 15)+
    scale_color_manual(values = c("#6a00fc","#ffa10c"))+
    labs(title = p, y = "Density", x = "Fragment length (bp)")+
    scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
    th
  return(gg)
}
install.packages("matrixStats")
####TPM of 15 genes####
x #data.frame of log2(TPM+1) values from RNA-seq of the 197 AVENIO genes returned by RNA_TPM()
x #name of colnames to be plotted
RNA_plot_genes <- function(x, y){
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(matrixStats)
  x1 <- x[c("SYMBOL", y)] %>% 
    filter(SYMBOL %in% c("NRAS", "RET", "KRAS",
                         "BRCA2","TP53", "ERBB2",
                         "BRCA1", "ALK","PDGFRA",
                         "KIT", "APC", "ROS1",
                         "EGFR", "MET", "BRAF"))
  
  activity <- c()
  sd <- c()
  x1$means <- rowMeans(x1[,-1])
  x1$sd <- rowSds(as.matrix(x1[,-1]))
  for (i in 1:length(x1$SYMBOL)){
    if(x1$means[i]>2){
      activity[i] <- "High"
    }
    else{
      activity[i] <- "Low"
    }
  }
  x1$activity <- activity
  x1 <- x1[order(x1$means),]
  gg <- ggplot(data = x1, aes(x = fct_inorder(SYMBOL), 
                              y = means, fill = activity))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin=means-sd, ymax=means+sd), width=.3)+
    labs(x = "", y = "Log2(TPM+1)")+
    scale_fill_manual("Gene activity", values = c("green4", "firebrick"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.7))+
    th
  return(gg)
}


####UMAP RNA-seq####
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
        theme(legend.position="bottom")+
        th
    return(gg)
}
gg_umap_RNA <- umap_RNA_seq(TPM,c("A549_", "HCC827_", "HCC827-MET_"))
####Volcano RNA-seq####
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
        scale_size_continuous(name = expression(log[2]("TPM+1")),
                              limits = c(0,15))+
        scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                             mid = "grey",
                             name = expression(log[2]("FC")))+
        guides(colour = guide_colourbar(order = 1),
               size = guide_legend(order = 2))+
        xlab(expression(log[2]("FC")))+
        ylab(expression(paste("-",log[10]("q-value"), sep = "")))+
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
####diff gene expression RNA_seq####
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
HCC827_vs_A549 <- dif_gene_express(TPM_AVENIO, c("HCC827_", "A549_"))
HCC827_vs_HCC827_MET <- dif_gene_express(TPM_AVENIO, c("HCC827_", "HCC827-MET_"))

####UMAP ChIP-seq####
x #transposed and bound samples returned by bind_samples()
umap_ChIP_seq <- function(x){
    library(umap)
    set.seed(8)
    df <- x %>% select(-groups)
    umap_res <- umap(df, n_neighbors = 5)
    umap_df <- as.data.frame(umap_res$layout)
    umap_df[,3] <- x$groups
    gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V3))+
        geom_point(size = 5)+
        labs(title = "ChIP-seq of cell lines",
             x = "UMAP-1", y = "UMAP-2")+
        scale_color_manual(name = "Cell line",
                           values = c("firebrick","#6a00fc","#ffa10c"))+
        theme_bw()+
        theme(legend.position="bottom")+
        th
    return(gg)
}


####diff gene enrichment ChIP####
x #mean enrichment of sample 1
y #mean enrichment of sample 2
p #list of data.frames returned by e.score for sample 1
q #list of data.frames returned by e.score for sample 2
b #optional, vector of bad genes
diff_gene_chip <- function(x,y,p,q, b = NULL){
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
    df <- data.frame(genes = x$genes, x = log2(x$enrichment), y = log2(y$enrichment))
    if (is.null(b)){
        df <- df
    }
    else{
        df <- df %>% filter(genes %ni% b)
    } 
    res <- summary(lm(y~x, data=df))
    a <- res$coefficients[2,1]
    c <- res$coefficients[1,1]
    df$norm <- df$y-((df$x*a)+c)
    df <- df[order(df$norm),]
    std <- p[[1]]$genes
    df1 <- data.frame(genes = p[[1]]$genes)
    for (i in 1:length(p)){
        p[[i]] <- p[[i]][match(std,p[[i]]$genes),]
        df1[ , ncol(df1) + 1] <- p[[i]]$enrichment
        colnames(df1)[ncol(df1)] <- paste0("Rx",i)
    }
    for (i in 1:length(q)){
        q[[i]] <- q[[i]][match(std,q[[i]]$genes),]
        df1[ , ncol(df1) + 1] <- q[[i]]$enrichment
        colnames(df1)[ncol(df1)] <- paste0("Ry",i)
    }
    if(!is.null(b)){
        df1 <- df1 %>% filter(genes %ni% b)
    }
    df1 <- na.omit(df1)
    options(scipen = 999)
    df <- df[match(df1$genes, df$genes),]
    df1 <- df1 %>% 
        mutate(norm_Rx1 = log2(Rx1)-((df$y*a)+c)) %>% 
        mutate(norm_Rx2 = log2(Rx2)-((df$y*a)+c)) %>% 
        mutate(norm_Rx3 = log2(Rx3)-((df$y*a)+c)) %>% 
        mutate(norm_Ry1 = log2(Ry1)-((df$y*a)+c)) %>% 
        mutate(norm_Ry2 = log2(Ry2)-((df$y*a)+c)) %>% 
        mutate(norm_Ry3 = log2(Ry3)-((df$y*a)+c)) 
    rownames(df1) <- df$genes
    df1 <- df1 %>% select(!genes)
    df1 <- df1[grepl("norm_",colnames(df1))]
    return(df1)
    group1_t <- as.data.frame(t(df1[1:length(x)]))
    group2_t <- as.data.frame(t(df1[(length(x)+1):(length(x)+length(y))]))
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
        p[i] <- t.test(g2, g1,
                       paired = FALSE,
                       alternative = "two.sided")$p.value
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
                     log2FC = log2(group1)-log2(group2))
    kk <- kk[order(-(kk$log2FC)),]
    return(kk)
    
}
diff_gene_chip(mean_HCC827_enrichment,
               mean_A549_enrichment,
               list(HCC827_R1_enrichment,
                    HCC827_R2_enrichment,
                    HCC827_R3_enrichment),
               list(A549_R1_enrichment,
                    A549_R2_enrichment,
                    A549_R3_enrichment),
               bad)
diff_gene_chip(list(HCC827_R1_enrichment,
                    HCC827_R2_enrichment,
                    HCC827_R3_enrichment),
               list(HCC827_MET_R1_enrichment,
                    HCC827_MET_R2_enrichment,
                    HCC827_MET_R3_enrichment),
               bad)


