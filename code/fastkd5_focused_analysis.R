# Load R libs ####
library(biomaRt)
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
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genelist <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                  filters = "hgnc_symbol",
                  values = colnames(cnv),
                  mart = mart)
chrs <- c(as.character(seq(1:22)), "X", "Y")
genelist <- genelist[genelist$chromosome_name %in% chrs,]
genelist$chromosome_name <- factor(genelist$chromosome_name, levels = chrs)


# Load mutation data ####
mut <- fread("/Users/lukas/Downloads/OmicsSomaticMutations.csv", sep = ",")
mut_orig <- mut


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


# Split by genotype ####
asplit <- split(meta$Run, meta$genotype)
means <- lapply(asplit, function(x)
  rowMeans(start_rate_ko[, x], na.rm = T))

subm <- data.frame(
  pos =  1:nrow(start_rate_ko),
  start_ko = means[["KD5 -/-"]],
  start_wt = means[["WT +/+"]]
)

means <- lapply(asplit, function(x)
  rowMeans(end_rate_ko[, x], na.rm = T))
subm$end_ko <- means[["KD5 -/-"]]
subm$end_wt <- means[["WT +/+"]]

means <- lapply(asplit, function(x)
  rowMeans(coverage_ko[, x], na.rm = T))
subm$cov_ko <- means[["KD5 -/-"]]
subm$cov_wt <- means[["WT +/+"]]

ggplot(subm[subm$cov_ko > 100 & subm$cov_wt > 100,], 
       aes(pos, start_ko - start_wt)) +
  labs(
    title = "Start rate difference",
    y = "KD5 - WT",
    x = "Position on chrM"
  ) +
  geom_point() +
  theme_classic() + xlim(9100, 9400)

ggplot(subm[subm$cov_ko > 100 & subm$cov_wt > 100,], 
       aes(pos, end_ko - end_wt)) +
  labs(
    title = "End rate difference",
    y = "KD5 - WT",
    x = "Position on chrM"
  ) +
  geom_point() +
  theme_classic() + xlim(9100, 9400)
  

# Junction based calculation ####
junctions <- unique(anno$V5)
tmp <- do.call(rbind, lapply(junctions, function(x)
  c(
    mean(subm$end_ko[(x - 10):x]) - mean(subm$end_wt[(x - 10):x]),
    mean(subm$start_ko[x:(x + 10)]) - mean(subm$start_wt[x:(x + 10)]))))
tmp <- data.frame(
  start = tmp[,2],
  end = tmp[,1],
  junction = junctions,
  id = paste(anno$genes[1:(nrow(anno))], 
             c(anno$genes[2:(nrow(anno))], "end"), 
             sep = "|"))

ggplot(tmp, aes(start, end)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = id)) +
  theme_classic()


tmp <- reshape2::melt(subm, measure.vars = c("start_ko", "start_wt"))
tmp <- reshape2::melt(subm, measure.vars = c("end_ko", "end_wt"))
ggplot(tmp[tmp$pos > 7400 & tmp$pos < 7500,],
       aes(pos, value, color = variable)) +
  geom_point() + theme_classic() +
  labs(title = "Start rate", y = "Rate", x = "Position on chrM")


# Single base calculation ####
res_start <- do.call(rbind, lapply(1:nrow(start_rate_ko), function(x){
  pval <- try(t.test(
    start_rate_ko[x, asplit[["KD5 -/-"]]],
    start_rate_ko[x, asplit[["WT +/+"]]])$p.value)
  if(class(pval) == "try-error") pval <- NA
  
  diff <- mean(start_rate_ko[x, asplit[["KD5 -/-"]]]) - 
    mean(start_rate_ko[x, asplit[["WT +/+"]]])
  
  c(diff, pval)
}))

res_end <- do.call(rbind, lapply(1:nrow(start_rate_ko), function(x){
  pval <- try(t.test(
    end_rate_ko[x, asplit[["KD5 -/-"]]],
    end_rate_ko[x, asplit[["WT +/+"]]])$p.value)
  if(class(pval) == "try-error") pval <- NA
  
  diff <- mean(end_rate_ko[x, asplit[["KD5 -/-"]]]) - 
    mean(end_rate_ko[x, asplit[["WT +/+"]]])
  
  c(diff, pval)
}))

res <- data.frame(
  pos = 1:nrow(start_rate_ko),
  start = res_start, end = res_end
)
res$start_score <- -log10(res$start.2)
res$start_score[which(res$start.1 < 0)] <- (-1)*res$start_score[which(res$start.1 < 0)]

res$end_score <- -log10(res$end.2)
res$end_score[which(res$end.1 < 0)] <- (-1)*res$end_score[which(res$end.1 < 0)]

ggplot(res, aes(pos, start_score)) +
  labs(title = "Start rate difference", y = "KD5 - WT", x = "Position on chrM") +
  geom_point() + 
  theme_classic()

ggplot(res, aes(pos, end_score)) +
  labs(title = "End rate", y = "KD5 - WT", x = "Position on chrM") +
  geom_point() + 
  theme_classic() + xlim(5800, 6100)

tmp <- reshape2::melt(res, measure.vars = c("start_score", "end_score"))

ggplot(tmp, aes(pos, value, color = value)) +
  labs(
    y = "KD5 vs WT",
    x = "Position on chrM"
  ) +
  facet_wrap(~variable, nrow = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(alpha = abs(value))) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
  theme_classic()

ggplot(tmp, aes(pos, value, color = value)) +
  labs(
    y = "KD5 vs WT",
    x = "Position on chrM"
  ) +
  facet_wrap(~variable, nrow = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(alpha = abs(value))) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
  theme_classic() +
  xlim(5850, 6050)

ggplot(tmp, aes(pos, value, color = value)) +
  labs(
    y = "KD5 vs WT",
    x = "Position on chrM"
  ) +
  facet_wrap(~variable, nrow = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(alpha = abs(value))) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
  theme_classic() +
  xlim(14600, 14900)

p_anno

junc <- 9207
subm <- data.frame(
  start = colMeans(start_rate_ko[junc:(junc + 10), ]),
  end = colMeans(end_rate_ko[(junc - 10):junc, ]))
meta$cleave <- subm$start + subm$end

ggplot(meta, aes(genotype, cleave)) +
  geom_boxplot() +
  theme_classic()


# Check DepMap ####
path <- "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/mt_cleavage_depmap/"
files <- list.files(path)
files <- files[grep("cleavage", files)]
files <- files[grep("SRR", files)]

dat <- lapply(files, function(x)
  read.delim(paste0(path, x)))
names(dat) <- gsub("_cleavage.tsv", "", fixed = T, basename(files))

sra <- read.delim(
  "/Users/lukas/OneDrive/Documents/GitHub/mt_cleavage/data/SraRunTable_v2.txt",
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

junc <- 9207
subm <- data.frame(
  start = colMeans(start_rate[junc:(junc + 10), ]),
  end = colMeans(end_rate[(junc - 10):junc, ]))
subm$id <- sra$DepMap_ID[match(rownames(subm), rownames(sra))]

ok <- intersect(subm$id, rownames(cnv))
correl <- cor(
  cnv[ok,], 
  subm[match(ok,subm$id), "start"] - subm[match(ok,subm$id), "end"], 
  use = 'pairwise.complete')[,1]

subm$cnv <- cnv[match(subm$id, rownames(cnv)), "FASTKD5"]
subm$score <- subm$start + subm$end

ggplot(subm, aes(cnv, start + end)) +
  labs(
    title = "DepMap",
    y = "ATP6-COX3 cleavage",
    x = "FASTKD5 CNV"
  ) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson") +
  geom_density_2d() +
  theme_classic()

ggplot(subm[!is.na(subm$cnv), ],
       aes(cnv < 0.8, start + end)) +
  labs(
    title = "DepMap",
    y = "ATP6-COX3 cleavage",
    x = "FASTKD5 CNV < 0.8"
  ) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.1) +
  ggpubr::stat_compare_means() +
  theme_classic()


# Check mutations ####
mut <- fread("/Users/lukas/Downloads/OmicsSomaticMutations.csv", sep = ",")
mut <- mut[mut$HugoSymbol == "FASTKD5",]

mut$score <- (subm$start + subm$end)[match(mut$ModelID, subm$id)]

mut$pos <- gsub("p.", "", fixed = T, mut$ProteinChange)
mut$pos <- unlist(lapply(mut$pos, function(x)
  substr(x, 2, nchar(x) - 1)))
mut$pos <- as.numeric(gsub("f$", "", mut$pos))

mut$ProteinChange <- factor(
  mut$ProteinChange, levels = unique(mut$ProteinChange[order(mut$pos)]))

ggplot(mut, aes(score, ProteinChange)) +
  labs(
    title = "DepMap",
    y = "FASTKD5 mutation",
    x = "ATP6-COX3 cleavage"
  ) +
  geom_boxplot() +
  geom_vline(xintercept = mean(subm$score), linetype = 2) +
  geom_vline(
    xintercept = mean(subm$score[which(subm$cnv < 0.8)]),
    linetype = 2, color = "green") +
  theme_classic()


subm$group <- subm$id %in% mut$ModelID
ggplot(subm, aes(interaction(group, cnv < 0.8), end + start)) +
  labs(
    title = "DepMap",
    y = "ATP6-COX3 cleavage"
  ) +
  geom_boxplot() +
  theme_classic()

subm$protein <- mut$ProteinChange[match(subm$id, mut$ModelID)]
subm$protein[subm$cnv < 0.8] <- "CNV loss"

ggplot(subm, aes(start + end, protein)) +
  labs(
    title = "DepMap",
    y = "ATP6-COX3 cleavage"
  ) +
  geom_boxplot() +
  theme_classic()
 

# Plot KO & DepMap ####
plt_both <- function(junc){
  subm <- data.frame(
    start = colMeans(start_rate[junc:(junc + 10), ]),
    end = colMeans(end_rate[(junc - 10):junc, ]))
  subm$id <- sra$DepMap_ID[match(rownames(subm), rownames(sra))]
  
  ok <- intersect(subm$id, rownames(cnv))
  correl <- cor(
    cnv[ok,], 
    subm[match(ok,subm$id), "start"] - subm[match(ok,subm$id), "end"], 
    use = 'pairwise.complete')[,1]
  
  subm$cnv <- cnv[match(subm$id, rownames(cnv)), "FASTKD5"]
  subm$score <- subm$start + subm$end
  
  p_depmap <- ggplot(subm[!is.na(subm$cnv), ],
                     aes(cnv < 0.8, score)) +
    labs(
      title = "DepMap",
      y = "Cleavage",
      x = "FASTKD5 CNV < 0.8"
    ) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.1) +
    ggpubr::stat_compare_means() +
    theme_classic()
  
  subm <- data.frame(
    start = colMeans(start_rate_ko[junc:(junc + 10), ]),
    end = colMeans(end_rate_ko[(junc - 10):junc, ]))
  meta$cleave <- subm$start + subm$end
  
  p_ko <- ggplot(meta, aes(genotype, cleave)) +
    labs(
      title = "KO data",
      y = "Cleavage",
      x = "Genotype"
    ) +
    geom_boxplot() +
    theme_classic()
  
  gridExtra::grid.arrange(p_depmap, p_ko, ncol = 2)
}
plt_both(5904)
plt_both(9207)
plt_both(14747)

plot_manhattan <- function(junc){
  subm <- data.frame(
    start = colMeans(start_rate[junc:(junc + 10), ]),
    end = colMeans(end_rate[(junc - 10):junc, ]))
  subm$id <- sra$DepMap_ID[match(rownames(subm), rownames(sra))]
  
  ok <- intersect(subm$id, rownames(cnv))
  correl <- cor(
    cnv[ok,], 
    subm[match(ok,subm$id), "start"] + subm[match(ok,subm$id), "end"], 
    use = 'pairwise.complete', method = "spearman")[,1]
  
  aframe <- data.frame(
    gene = names(correl),
    cor = correl,
    genelist[match(names(correl), genelist$hgnc_symbol),]
  )
  aframe <- aframe[order(aframe$chromosome_name, aframe$start_position),]
  aframe$rank <- 1:nrow(aframe)
  
  ggplot(aframe[-which(aframe$chromosome_name %in% c("X", "Y") |
                         is.na(aframe$chromosome_name)), ],
         aes(rank, cor, color = chromosome_name)) +
    labs(
      title = "DepMap",
      y = "Correlation",
      x = "Genes sorted by genomic location"
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    geom_text(
      data = aframe[aframe$gene == "FASTKD5", ],
      aes(label = gene), color = "red", size = 5
    ) +
    theme_classic() +
    scale_color_manual(values = rep(c("black", "grey"), 11)) +
    theme(legend.position = "none")
}
plot_manhattan(5904)
plot_manhattan(9207)
plot_manhattan(14747)
