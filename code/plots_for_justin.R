# Load R libs ####
library(data.table)
library(ggplot2)
library(limma)
library(pheatmap)

# Load MT annotation ####
anno <- read.delim("/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/mt_gencode.gtf")
anno$genes <- unlist(lapply(anno$V9, function(x)
  substr(x, regexpr("name", x)[1] + 5, gregexpr(";", x)[[1]][3] - 1)))

p_anno <- ggplot() +
  labs(y = "annotation", 
       x = "Position on chrM") +
  xlim(0, 16569) +
  scale_color_manual(values = c("black")) +
  geom_rect(
    data = anno[anno$V7 == "+",],
    aes(xmin = V4, xmax = V5, ymin = 0.1, ymax = 1,
        alpha = .1, color = 'black')) +
  geom_rect(
    data = anno[anno$V7 == "-",],
    aes(xmin = V4, xmax = V5, ymin = -1, ymax = -0.1,
        alpha = .1, color = 'black')) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.y = element_blank())


# Load DepMap ####
load("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/global.RData")


# Get coordinates ####
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genelist <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                  filters = "hgnc_symbol",
                  values = colnames(cnv),
                  mart = mart)
chrs <- c(as.character(seq(1:22)), "X", "Y")
genelist <- genelist[genelist$chromosome_name %in% chrs,]


# Load mutation data ####
mut <- fread("/Users/lukas/Downloads/OmicsSomaticMutations.csv", sep = ",")
mut_orig <- mut


# Load cleavage data ####
path <- "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/mt_cleavage_depmap/"
files <- list.files(path)
files <- files[grep("cleavage", files)]
files <- files[grep("SRR", files)]

dat <- lapply(files, function(x)
  read.delim(paste0(path, x)))
names(dat) <- gsub("_cleavage.tsv", "", fixed = T, basename(files))

sra <- read.csv(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/depmap_info.txt",
  row.names = 1)
sra <- sra[names(dat), ]


# Calculate start/end rates ####
coverage <- do.call(cbind, lapply(dat, function(x) x$total_num_reads))
starts <- do.call(cbind, lapply(dat, function(x) x$num_of_starts))
ends <- do.call(cbind, lapply(dat, function(x) x$num_of_ends))

start_rate <- do.call(cbind, lapply(dat, function(x)
  x$num_of_starts/x$total_num_reads))
end_rate <- do.call(cbind, lapply(dat, function(x)
  x$num_of_ends/x$total_num_reads))

names(cleave_eff) <- sra$Sample.Name[match(names(cleave_eff),
                                           rownames(sra))]
sample_info$cleave <- cleave_eff[
  match(sample_info$CCLE_Name, names(cleave_eff))]


# Define MT processing genes ####
genes <- c("ANGEL2",
           "ELAC1",
           "H1RNA",
           "FASTKD1",
           "FASTKD2",
           "FASTKD5",
           "FASTKD3",
           "FASTK",
           "TBRG4",
           "STAT3",
           "LRPPRC",
           "SLIRP",
           "PTPN1",
           "SUPV3L1",
           "PRORP",
           "TRMT10C",
           "SDR5C1",
           "NLRX1")


# Load FASTK knockout data ####
files <- list.files(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/mt_cleavage_GSE156260/", full.names = T)
dat <- lapply(files, read.delim)
names(dat) <- gsub("_1_preprocessed_cleavage.tsv", "", fixed = T, basename(files))

meta <- read.csv(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/GSE156260_info.txt")
meta <- meta[match(names(dat), meta$Run), ]
rownames(meta) <- meta$Run

start_rate_ko <- do.call(cbind, lapply(dat, function(x)
  x$num_of_starts/x$total_num_reads))
end_rate_ko <- do.call(cbind, lapply(dat, function(x)
  x$num_of_ends/x$total_num_reads))

coverage_ko <- do.call(cbind, lapply(dat, function(x) x$total_num_reads))
starts_ko <- do.call(cbind, lapply(dat, function(x) x$num_of_starts))
ends_ko <- do.call(cbind, lapply(dat, function(x) x$num_of_ends))


# Run PCA ####
tmp <- start_rate_ko[which(apply(start_rate_ko, 1, var) > 0), ]
nas <- apply(tmp, 1, function(x) mean(is.na(x)))
tmp <- tmp[nas == 0,]
tmp <- tmp[apply(tmp, 1, var) > 0,]
pca <- prcomp(t(tmp), scale. = T)
aframe <- data.frame(pca$x, meta)
p1 <- ggplot(aframe, aes(PC1, PC2, color = genotype)) +
  labs(title = "Start rate") +
  geom_point() + theme_classic()

tmp <- end_rate_ko[which(apply(end_rate_ko, 1, var) > 0), ]
nas <- apply(tmp, 1, function(x) mean(is.na(x)))
tmp <- tmp[nas == 0,]
tmp <- tmp[apply(tmp, 1, var) > 0,]
pca <- prcomp(t(tmp), scale. = T)
aframe <- data.frame(pca$x, meta)
p2 <- ggplot(aframe, aes(PC1, PC2, color = genotype)) +
  labs(title = "End rate") +
  geom_point() + theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 2)


# Calculate gene specific cleavage rates ####
gene_rates_ko <- do.call(cbind, lapply(1:nrow(anno), function(x){
  from <- anno$V4[x]
  to <- anno$V4[x] + 10
  start_rate <- apply(start_rate_ko[from:to,], 2, function(k)
    mean(k, na.rm = T))
  
  from <- anno$V5[x]
  to <- anno$V5[x] - 10
  end_rate <- apply(end_rate_ko[from:to,], 2, function(k)
    mean(k, na.rm = T))
  
  (start_rate + end_rate)/2
}))
nom <- unlist(lapply(anno$V9, function(x){
  substr(x,
         gregexpr("gene_name", x)[[1]] + 10,
         gregexpr(";", x)[[1]][3] - 1)
}))
colnames(gene_rates_ko) <- nom

aframe <- data.frame(
  gene_rates_ko,
  global = rowMeans(gene_rates_ko, na.rm = T),
  meta)

ggplot(aframe, 
       aes(genotype, global)) +
  geom_boxplot() +
  theme_classic()

ggplot(aframe[aframe$genotype %in% c("KD5 -/-", "KD4 -/-", "WT +/+"), ], 
       aes(genotype, MT.CO1)) +
  geom_boxplot() +
  theme_classic()

ggplot(aframe[aframe$genotype %in% c("KD5 -/-", "KD4 -/-", "WT +/+"), ], 
       aes(genotype, MT.CYB)) +
  geom_boxplot() +
  theme_classic()

gene_rates <- do.call(cbind, lapply(1:nrow(anno), function(x){
  from <- anno$V4[x]
  to <- anno$V4[x] + 10
  start_rate <- apply(start_rate[from:to,], 2, function(k)
    mean(k, na.rm = T))
  
  from <- anno$V5[x]
  to <- anno$V5[x] - 10
  end_rate <- apply(end_rate[from:to,], 2, function(k)
    mean(k, na.rm = T))
  
  (start_rate + end_rate)/2
}))
nom <- unlist(lapply(anno$V9, function(x){
  substr(x,
         gregexpr("gene_name", x)[[1]] + 10,
         gregexpr(";", x)[[1]][3] - 1)
}))
colnames(gene_rates) <- nom

correl <- cor(gene_rates, use = 'pairwise.complete')
pheatmap(correl,
         breaks = seq(-0.75, 0.75, length = 100),
         main = "Correlation of processing rates across cell lines")

gene_res <- apply(gene_rates, 2, function(cleave_eff){
  names(cleave_eff) <- sra$Sample.Name[match(names(cleave_eff),
                                             rownames(sra))]
  sample_info$cleave <- cleave_eff[
    match(sample_info$CCLE_Name, names(cleave_eff))]
  
  ok <- intersect(rownames(sample_info)[!is.na(sample_info$cleave)],
                  rownames(cnv))
  
  correl <- apply(cnv[ok,], 2, function(x)
    cor.test(x, sample_info$cleave[match(ok, sample_info$DepMap_ID)]))
  
  res <- do.call(rbind, lapply(correl, function(x)
    unlist(x[c("estimate", "p.value")])))
  res <- data.frame(res)
  colnames(res) <- c("coef", "pval")
  res
})


# Create heatmap ####
tmp <- do.call(rbind, lapply(gene_res, function(res_cnv){
  subm <- res_cnv[match(genes, rownames(res_cnv)), ]
  score <- -log10(subm[ , 'pval'])
  score[which(subm[, 'coef'] < 0)] <- (-1)*score[which(subm[, 'coef'] < 0)]
  score
}))
colnames(tmp) <- genes
tmp <- na.omit(t(tmp))

fastk_family <- c("FASTK", "FASTKD1", "FASTKD2", 
                  "FASTKD3", "FASTKD5", "TBRG4")
tmp <- tmp[fastk_family, -grep("MT-T", colnames(tmp))]

tmp <- tmp[, -which(colnames(tmp) %in% c("MT-RNR1", "MT-RNR2"))]

pheatmap(tmp,
         cluster_cols = F,
         breaks = seq(-20, 20, length = 100),
         color = viridisLite::viridis(100),
         main = "Correlation of CNV with gene-specific cleavage")

pheatmap(tmp[c("TBRG4", "FASTKD5"), ],
         cluster_cols = F,
         breaks = seq(-20, 20, length = 100),
         main = "Correlation of CNV with gene-specific cleavage")


# Create GWAS plot - global rate ####
cleave_eff <- rowMeans(
  apply(gene_rates, 2, function(x) (x - mean(x))/sd(x)), na.rm = T)
names(cleave_eff) <- sra$Sample.Name[match(names(cleave_eff),
                                           rownames(sra))]
sample_info$cleave <- cleave_eff[
  match(sample_info$CCLE_Name, names(cleave_eff))]

ok <- intersect(rownames(sample_info)[!is.na(sample_info$cleave)],
                rownames(cnv))

correl <- apply(cnv[ok,], 2, function(x)
  cor.test(x, sample_info$cleave[match(ok, sample_info$DepMap_ID)]))

res <- do.call(rbind, lapply(correl, function(x)
  unlist(x[c("estimate", "p.value")])))
res <- data.frame(res)
colnames(res) <- c("coef", "pval")

res$score <- -log10(res$pval)
res$score[which(res$coef < 0)] <- (-1)*res$score[which(res$coef < 0)]

res <- data.frame(
  res, genelist[match(rownames(res), genelist$hgnc_symbol), ])

res$chromosome_name <- factor(
  res$chromosome_name,
  levels = chrs) 
res <- res[order(res$chromosome_name, res$start_position), ]
res$rank <- 1:nrow(res)
res <- res[!is.na(res$chromosome_name), ]

res$label <- NA
res$label[rownames(res) %in% genes] <- rownames(res)[rownames(res) %in% genes]

ggplot(res,
       aes(rank, score, color = chromosome_name, label = label)) +
  geom_point() +
  labs(title = "CNV association with global cleavage",
       y = "Association score",
       x = "Genomic positions") +
  geom_label(size = 5, color = "black") +
  theme_classic()


# Create GWAS plot - MT-CO1 ####
res <- gene_res[["MT-CO1"]]
res <- data.frame(
  genelist[match(rownames(res), genelist$hgnc_symbol), ],
  res
)
res$score <- -log10(res$pval)
res$score[which(res$coef < 0)] <- (-1)*res$score[which(res$coef < 0)]
res$chromosome_name <- factor(
  res$chromosome_name,
  levels = chrs) 
res <- res[order(res$chromosome_name, res$start_position), ]
res$rank <- 1:nrow(res)
res <- res[!is.na(res$chromosome_name), ]

res$label <- NA
res$label[which(res$hgnc_symbol == "FASTKD5")] <- "FASTKD5"
res$label[which(res$hgnc_symbol == "TBRG4")] <- "FASTKD4"

ggplot(res,
       aes(rank, score, color = chromosome_name, label = label)) +
  geom_point() +
  labs(title = "CNV association with MT-CO1 cleavage",
       y = "Association score",
       x = "Genomic positions") +
  geom_label(size = 5, color = "black") +
  theme_classic()


# Plot FASTKD5 CNV vs MT-CO1 cleavage ####
cleave_eff <- gene_rates[, "MT-CO1"]
cleave_eff <- rowMeans(gene_rates)
names(cleave_eff) <- sra$Sample.Name[match(names(cleave_eff),
                                           rownames(sra))]
sample_info$cleave <- cleave_eff[
  match(sample_info$CCLE_Name, names(cleave_eff))]

aframe <- data.frame(
  sample_info,
  cnv = cnv[match(rownames(sample_info), rownames(cnv)), "FASTKD5"]
)
aframe$group <- "WT"
aframe$group[aframe$cnv > 1.25] <- "gain"
aframe$group[aframe$cnv < 0.75] <- "loss"
aframe$group <- factor(aframe$group, levels = c("loss", "WT", "gain"))

ggplot(aframe, aes(
  group, cleave)) +
  labs(
    title = "DepMap",
    subtitle = "FASTKD5 copy number vs MT-CO1 cleavage",
    y = "MT-CO1 cleavage rate",
    x = "FASTKD5 copy number") +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()


# Correlate RNA expression w cleavage in WT ####
ok <- aframe$DepMap_ID[which(aframe$group == "WT" & 
                               !is.na(aframe$cleave))]
ok <- intersect(ok, rownames(rnaseq))

res <- apply(rnaseq[match(ok, rownames(rnaseq)), ], 2, function(x){
  subm <- data.frame(
    gene = x,
    aframe[match(ok, rownames(aframe)), ])
  lm(x ~ lineage + cleave, data = subm)
})
res <- do.call(rbind, lapply(res, function(x)
  coefficients(summary(x))["cleave", c(1, 4)]
))
res <- data.frame(res)
colnames(res) <- c("coef", "pval")
res <- res[order(res$pval), ]

res$label <- NA
res$label[1:50] <- rownames(res)[1:50]

ggplot(res, aes(coef, -log10(pval), label = label)) +
  geom_point() +
  labs(
    title = "DepMap",
    subtitle = "RNA expression associated with cleavage",
    y = "-log10(p-value)",
    x = "Coefficient") +
  ggrepel::geom_text_repel(max.overlaps = 10) +
  theme_classic()


# Create heatmap ####
anno_col <- aframe[, c("lineage", "cleave", "group")]
good <- tail(names(sort(table(anno_col$lineage))), 10)
anno_col$lineage[!anno_col$lineage %in% good] <- "Other"

pheatmap(
  resids[match(sig, rownames(resids)),
         order(aframe$cleave[match(colnames(resids), rownames(aframe))])],
         breaks = seq(-2, 2, length = 100),
         scale = "row",
         cluster_cols = F,
         annotation_col = anno_col,
         show_colnames = F)


# Correlate with dependency ####
ok <- intersect(
  rownames(kronos),
  sample_info$DepMap_ID[!is.na(sample_info$cleave)])

correl <- cor(
  kronos[ok,], 
  sample_info$cleave[match(ok, rownames(sample_info))],
  use = "pairwise.complete.obs")[,1]
correl <- sort(correl)

library(enrichR)
up <- tail(names(correl), 500)
enrich <- enrichr(
  up,
  databases = c(
    "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
    "Disease_Signatures_from_GEO_down_2014",
    "Disease_Signatures_from_GEO_up_2014"))



