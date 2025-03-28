# Load R libs ####
library(biomaRt)
library(data.table)
library(ggpubr)
library(ggplot2)
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
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genelist <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                  filters = "hgnc_symbol",
                  values = colnames(cnv),
                  mart = mart)
chrs <- c(as.character(seq(1:22)), "X", "Y")
genelist <- genelist[genelist$chromosome_name %in% chrs,]
genelist$chromosome_name <- factor(
  genelist$chromosome_name, levels = chrs
)


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


# Calculate gene specific cleavage rates ####
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


# Make waterfall plot of cleavage rate of select gene ####
subm <- data.frame(
  sra,
  rate = gene_rates[, "MT-CYB"]
)
subm$rank <- rank(subm$rate)
ggplot(subm, aes(x = rank, y = rate)) +
  geom_point() +
  labs(
    title = "DepMap",
    x = "Cell lines ranked by cleavage rate",
    y = "MT-CYB cleavage rate") +
  theme_classic()

correl <- cor(gene_rates, use = 'pairwise.complete')
pheatmap(correl,
         breaks = seq(-0.75, 0.75, length = 100),
         main = "Correlation of processing rates across cell lines")


# Correlate CNV w cleavage of each gene ####
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


# Create Manhattan ####
subm <- gene_res[["MT-CYB"]]
subm$score <- -log10(subm[ , 'pval'])
subm$score[which(subm[, 'coef'] < 0)] <- (-1)*subm$score[which(subm[, 'coef'] < 0)]

subm <- data.frame(
  genelist[match(rownames(subm), genelist$hgnc_symbol),],
  subm
)
subm <- subm[order(subm$chromosome_name, subm$start_position), ]
subm$rank <- 1:nrow(subm)

ggplot(subm[!is.na(subm$chromosome_name), ], 
       aes(x = rank, y = score, color = chromosome_name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  labs(
    title = "DepMap",
    x = "Genes ranked by genomic location",
    y = "MT-CYB association score") +
  scale_color_manual(
    values = rep(c("black", "grey"), 22)) +
  theme_classic() +
  theme(legend.position = "none")

fastk_family <- c("FASTK", "FASTKD1", "FASTKD2", 
                  "FASTKD3", "FASTKD5", "TBRG4")

ggplot(subm[!is.na(subm$chromosome_name), ], 
       aes(x = rank, y = score, color = chromosome_name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  labs(
    title = "DepMap",
    x = "Genes ranked by genomic location",
    y = "MT-CYB association score") +
  scale_color_manual(
    values = rep(c("black", "grey"), 22)) +
  geom_text(
    data = subm[subm$hgnc_symbol %in% fastk_family, ],
    aes(label = hgnc_symbol),
    color = "red") +
  theme_classic() +
  theme(legend.position = "none")


# Plot TRBG4 CNV vs MT cleavage ####
sra$cleave <- gene_rates[, "MT-CYB"]

aframe <- data.frame(
  sample_info[match(sra$Sample.Name, sample_info$CCLE_Name),],
  sra
)
aframe$cnv <- cnv[match(aframe$DepMap_ID, rownames(cnv)), "TBRG4"]
aframe$group <- "WT"
aframe$group[aframe$cnv > 1.2] <- "gain"
aframe$group[aframe$cnv < 0.8] <- "loss"
aframe$group <- factor(aframe$group, levels = c("loss", "WT", "gain"))

ggplot(aframe, aes(
  group, cleave)) +
  labs(
    title = "DepMap",
    subtitle = "TBRG4 copy number vs MT cleavage",
    y = "MT-CYB cleavage rate",
    x = "TBRG4 copy number") +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()


aframe <- data.frame(
  sample_info[match(sra$Sample.Name, sample_info$CCLE_Name),],
  rate = gene_rates[, c("MT-CYB", "MT-CO1")]
)
aframe <- data.frame(
  aframe,
  cnv[match(aframe$DepMap_ID, rownames(cnv)), c("TBRG4", "FASTKD5")]
)

tmp <- reshape2::melt(aframe, measure.vars = c("TBRG4", "FASTKD5"))

p1 <- ggplot(tmp[!is.na(tmp$value), ], aes(
  value < 0.8, rate.MT.CYB)) +
  facet_wrap(~ variable) +
  labs(
    title = "DepMap",
    subtitle = "TBRG4 copy number vs MT cleavage",
    y = "MT-CYB cleavage rate",
    x = "Copy number loss") +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()

p2 <- ggplot(tmp[!is.na(tmp$value), ], aes(
  value < 0.8, rate.MT.CO1)) +
  facet_wrap(~ variable) +
  labs(
    title = "DepMap",
    subtitle = "TBRG4 copy number vs MT cleavage",
    y = "MT-CO1 cleavage rate",
    x = "Copy number loss") +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 1)


# Load FASTK knockout data ####
files <- list.files(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/mt_cleavage_GSE156260/", full.names = T)
dat <- lapply(files, read.delim)
names(dat) <- gsub("_1_preprocessed_cleavage.tsv", "", fixed = T, basename(files))

meta <- read.csv(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/GSE156260_info.txt")
meta <- meta[match(names(dat), meta$Run), ]
rownames(meta) <- meta$Run

start_rate <- do.call(cbind, lapply(dat, function(x)
  x$num_of_starts/x$total_num_reads))
end_rate <- do.call(cbind, lapply(dat, function(x)
  x$num_of_ends/x$total_num_reads))

gene_rates_ko <- do.call(cbind, lapply(1:nrow(anno), function(x){
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
colnames(gene_rates_ko) <- nom

aframe <- data.frame(
  meta,
  rate = gene_rates_ko[, c("MT-ATP6", "MT-ATP8")]
)

conditions <- c(
  "HAP1_KD4 -/-", 
  "HAP1_KD5 -/-", 
  "HAP1_WT +/+")

p1 <- ggplot(aframe[aframe$source_name %in% conditions,], aes(
  source_name, rate.MT.ATP8)) +
  labs(
    title = "DepMap",
    y = "MT.ATP8 cleavage rate",
    x = "Condition") +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(
    comparisons = list(c(1, 3), c(2, 3)),
    method = "t.test") +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()

p2 <- ggplot(aframe[aframe$source_name %in% conditions,], aes(
  source_name, rate.MT.ATP6)) +
  labs(
    title = "DepMap",
    y = "MT.ATP6 cleavage rate",
    x = "Condition") +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(
    comparisons = list(c(1, 3), c(2, 3)),
    method = "t.test") +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 1)


