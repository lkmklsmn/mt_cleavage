# Load R libs ####
library(data.table)
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

# Calculate per base deltas ####
res <- t(apply(start_rate_ko, 1, function(x){
  subm <- data.frame(rate = x, meta)
  subm <- subm[subm$genotype %in% c('KD4 -/-', 'WT +/+'), ]
  coefs <- try(coefficients(summary(lm(rate ~ genotype, data = subm))))
  if(class(coefs) == "try-error") return(c(NA, NA))
  coefs[2, c(1, 4)]
}))
res <- data.frame(pos = 1:nrow(res), res)
colnames(res) <- c("pos", "coef", "pval")
res$score <- -log10(res$pval)
res$score[res$coef < 0] <- (-1)*res$score[res$coef < 0]
res_kd4 <- res

cells <- sample_info$DepMap_ID[match(sra$Sample.Name, sample_info$CCLE_Name)]
correl <- cor(cnv[match(cells, rownames(cnv)), "TBRG4"],
              t(start_rate), use = 'pairwise.complete')[1,]
res_kd4$correl <- correl

tmp <- data.frame(score = zoo::rollmean(res_kd4$score, 5, na.pad = T),
                  correl = zoo::rollmean(res_kd4$correl, 5, na.pad = T))
tmp$pos <- 1:nrow(tmp)

ggplot(tmp, aes(score, correl)) +
  geom_point(alpha = .2) +
  geom_density_2d() +
  labs(x = "delta KD4-/- vs WT",
       y = "correlation w TBRG4 CNV") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic()

p1 <- ggplot(tmp, aes(pos, score, color = score, alpha = abs(score))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed') +
  labs(y = "delta KD4-/- vs WT",
       x = "") +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p2 <- ggplot(tmp, aes(pos, correl, color = correl, alpha = abs(correl))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed') +
  labs(y = "correlation w TBRG4 CNV",
       x = "") +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

gridExtra::grid.arrange(p1, p2, p_anno, ncol = 1,
                        heights = c(3, 3, 1))


res_kd4[which(res_kd4$correl > 0.1 & res_kd4$score > 2),]

aframe <- data.frame(rate = start_rate[8375, ],
                     cnv = cnv[match(cells, rownames(cnv)), "TBRG4"])
aframe$group <- "WT"
aframe$group[aframe$cnv < 0.8] <- "loss"
aframe$group[aframe$cnv > 1.2] <- "gain"
aframe$group <- factor(aframe$group,
                       levels = c("loss", "WT", "gain"))

ggplot(aframe[aframe$group != "WT", ],
       aes(group, rate)) +
  labs(title = "DepMap",
       x = "TBRG4 CNV",
       y = "Start rate at position 8375") +
  geom_boxplot() + geom_point() +
  theme_classic()

aframe <- data.frame(rate = start_rate_ko[8375, ],
                     meta)
ggplot(aframe[aframe$genotype %in% c("KD4 -/-", "WT +/+"), ],
       aes(genotype, rate)) +
  labs(title = "PLOS paper",
       x = "Condition",
       y = "Start rate at position 8375") +
  geom_boxplot() + geom_point() +
  theme_classic()

res <- t(apply(start_rate_ko, 1, function(x){
  subm <- data.frame(rate = x, meta)
  subm <- subm[subm$genotype %in% c('KD5 -/-', 'WT +/+'), ]
  coefs <- try(coefficients(summary(lm(rate ~ genotype, data = subm))))
  if(class(coefs) == "try-error") return(c(NA, NA))
  coefs[2, c(1, 4)]
}))
res <- data.frame(pos = 1:nrow(res), res)
colnames(res) <- c("pos", "coef", "pval")
res$score <- -log10(res$pval)
res$score[res$coef < 0] <- (-1)*res$score[res$coef < 0]
res_kd5 <- res

cells <- sample_info$DepMap_ID[match(sra$Sample.Name, sample_info$CCLE_Name)]
correl <- cor(cnv[match(cells, rownames(cnv)), "FASTKD5"],
              t(start_rate), use = 'pairwise.complete')[1,]
res_kd5$correl <- correl

tmp <- data.frame(score = zoo::rollmean(res_kd5$score, 5, na.pad = T),
                  correl = zoo::rollmean(res_kd5$correl, 5, na.pad = T))
tmp$pos <- 1:nrow(tmp)

p1 <- ggplot(tmp, aes(pos, score, color = score, alpha = abs(score))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed') +
  labs(y = "delta KD5-/- vs WT",
       x = "") +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p2 <- ggplot(tmp, aes(pos, correl, color = correl, alpha = abs(correl))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed') +
  labs(y = "correlation w FASTKD5 CNV",
       x = "") +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

gridExtra::grid.arrange(p1, p2, p_anno, ncol = 1,
                        heights = c(3, 3, 1))

ggplot(tmp, aes(score, correl)) +
  geom_point(alpha = .2) +
  geom_density_2d() +
  labs(x = "delta KD5-/- vs WT",
       y = "correlation w TBRG4 CNV") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic()


aframe <- data.frame(
  cbind(start_rate_ko[,1],
        end_rate_ko[,1]))
aframe$pos <- 1:nrow(aframe)

ggplot(aframe, aes(x = pos)) +
  geom_point(aes(y = X1, color = "start")) +
  geom_point(aes(y = -X2, color = "end")) +
  theme_classic()

tmp <- data.frame(start = zoo::rollmean(start_rate_ko[,1], 5, na.pad = T),
                  end = zoo::rollmean(end_rate_ko[,1], 5, na.pad = T))
tmp$pos <- 1:nrow(tmp)

ggplot(aframe[(3230-100):(3230 + 100), ],
       aes(x = pos)) +
  labs(title = "PLOS paper",
       subtitle = "WT",
       y = "Rate",
       x = "Position on MT") +
  geom_point(aes(y = X1, color = "start")) +
  geom_point(aes(y = -X2, color = "end")) +
  geom_vline(xintercept = c(3230, 3307)) +
  theme_classic()

r <- sample(ncol(start_rate_ko), 1)
tmp <- data.frame(start = zoo::rollmean(start_rate[,r], 5, na.pad = T),
                  end = zoo::rollmean(end_rate[,r], 5, na.pad = T))
tmp$pos <- 1:nrow(tmp)

ggplot(tmp[(3230-100):(3230 + 100), ],
       aes(x = pos)) +
  labs(title = "DepMap",
       y = "Rate",
       x = "Position on MT") +
  geom_point(aes(y = start, color = "start")) +
  geom_point(aes(y = -end, color = "end")) +
  geom_vline(xintercept = c(3230, 3307)) +
  theme_classic()

# Calculate cleavage rate by junction ####
junx_rates <- do.call(cbind, lapply(1:nrow(anno), function(x){
  pos <- anno$V4[x]
  num_starts <- apply(starts[pos:(pos+20), ], 2, sum)
  num_ends <- apply(ends[(pos-20):pos, ], 2, sum)
  coverage <- apply(coverage[(pos-20):(pos+20), ], 2, mean)
  rate <- (num_starts + num_ends)/coverage
}))
rownames(junx_rates) <- colnames(starts)

junx_rates_ko <- do.call(cbind, lapply(1:nrow(anno), function(x){
  pos <- anno$V4[x]
  num_starts <- apply(starts_ko[pos:(pos+20), ], 2, sum)
  num_ends <- apply(ends_ko[(pos-20):pos, ], 2, sum)
  coverage <- apply(coverage_ko[(pos-20):(pos+20), ], 2, mean)
  rate <- (num_starts + num_ends)/coverage
}))

res_j_kd4 <- t(apply(junx_rates_ko, 2, function(x){
  subm <- data.frame(rate = x, meta)
  subm <- subm[subm$genotype %in% c('KD4 -/-', 'WT +/+'), ]
  coefs <- try(coefficients(summary(lm(rate ~ genotype, data = subm))))
  if(class(coefs) == "try-error") return(c(NA, NA))
  coefs[2, c(1, 4)]
}))
res_j_kd4 <- data.frame(res_j_kd4)
colnames(res_j_kd4) <- c("coef", "pval")
res_j_kd4$score <- -log10(res_j_kd4$pval)
res_j_kd4$score[res_j_kd4$coef < 0] <- (-1)*res_j_kd4$score[res_j_kd4$coef < 0]

cells <- sample_info$DepMap_ID[match(sra$Sample.Name, sample_info$CCLE_Name)]
correl <- cor(cnv[match(cells, rownames(cnv)), "TBRG4"],
              junx_rates, use = 'pairwise.complete')[1,]
res_j_kd4$correl <- correl

res_j_kd4$pos <- anno$V4
res_j_kd4$gene <- anno$genes

res_j_kd5 <- t(apply(junx_rates_ko, 2, function(x){
  subm <- data.frame(rate = x, meta)
  subm <- subm[subm$genotype %in% c('KD5 -/-', 'WT +/+'), ]
  coefs <- try(coefficients(summary(lm(rate ~ genotype, data = subm))))
  if(class(coefs) == "try-error") return(c(NA, NA))
  coefs[2, c(1, 4)]
}))
res_j_kd5 <- data.frame(res_j_kd5)
colnames(res_j_kd5) <- c("coef", "pval")
res_j_kd5$score <- -log10(res_j_kd5$pval)
res_j_kd5$score[res_j_kd5$coef < 0] <- (-1)*res_j_kd5$score[res_j_kd5$coef < 0]

cells <- sample_info$DepMap_ID[match(sra$Sample.Name, sample_info$CCLE_Name)]
correl <- cor(cnv[match(cells, rownames(cnv)), "FASTKD5"],
              junx_rates, use = 'pairwise.complete')[1,]
res_j_kd5$correl <- correl

res_j_kd5$pos <- anno$V4
res_j_kd5$gene <- anno$genes

tmp <- data.frame(
  kd4 = res_j_kd4,
  kd5 = res_j_kd5
)
ggplot(tmp, aes(kd4.score, kd4.correl, label = kd4.gene)) +
  labs(title = "KD4",
       subtitle = "",
       x = "Delta KD4-/- vs WT",
       y = "Correlation w TBRG4 CNV") +
  ggpubr::stat_cor() +
  geom_label() +
  theme_classic()

ggplot(tmp, aes(kd5.score, kd5.correl, label = kd5.gene)) +
  labs(title = "KD5",
       subtitle = "",
       x = "Delta KD5-/- vs WT",
       y = "Correlation w FASTDK5 CNV") +
  ggpubr::stat_cor() +
  geom_label() +
  theme_classic()

aframe <- data.frame(junx_rates_ko, meta)
ggplot(aframe,
       aes(genotype, X24)) +
  labs(title = "PLOS paper",
       y = "Junction rate [start MT-TG]") +
  geom_boxplot() +
  theme_classic()

aframe <- data.frame(junx_rates,
                     cnv = cnv[match(cells, rownames(cnv)), "TBRG4"])
ggplot(aframe,
       aes(cnv, X24)) +
  labs(title = "DepMap",
       y = "Junction rate [start MT-TG]",
       x = "TBRG4 CNV") +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor() +
  theme_classic()

res_j_kd4$pos <- anno$V4

p1 <- ggplot(res_j_kd4, aes(pos, score, color = score)) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed', alpha = 0.1) +
  labs(y = "delta KD4-/- vs WT",
       x = "") +
  xlim(0, nrow(coverage)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p2 <- ggplot(res_j_kd4, aes(pos, correl, color = correl, alpha = abs(correl))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed', alpha = 0.1) +
  labs(y = "correlation w TBRG4 CNV",
       x = "") +
  xlim(0, nrow(coverage)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
  
gridExtra::grid.arrange(p1, p2, p_anno, ncol = 1,
                        heights = c(3, 3, 1))

p1 <- ggplot(res_j_kd5, aes(pos, score, color = score)) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed', alpha = 0.1) +
  labs(y = "delta KD5-/- vs WT",
       x = "") +
  xlim(0, nrow(coverage)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p2 <- ggplot(res_j_kd5, aes(pos, correl, color = correl, alpha = abs(correl))) +
  geom_vline(xintercept = anno$V4, linetype = 'dashed', alpha = 0.1) +
  labs(y = "correlation w FASTKD5 CNV",
       x = "") +
  xlim(0, nrow(coverage)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

gridExtra::grid.arrange(p1, p2, p_anno, ncol = 1,
                        heights = c(3, 3, 1))


# Correlate cleavage with TBRG4 mutations ####
mut <- mut_orig

colnames(junx_rates) <- gsub("-", ".", fixed = T, anno$gene)

clean_up <- function(mut){
  mut$name <- sample_info$CCLE_Name[match(mut$ModelID, sample_info$DepMap_ID)]
  mut$sra <- rownames(sra)[match(mut$name, sra$Sample.Name)]
  mut <- mut[!is.na(mut$sra), ]
  mut$prot_pos <- unlist(lapply(mut$ProteinChange, function(x)
    as.numeric(substr(x, 4, nchar(x) - 1))))
  
  aframe <- data.frame(
    mut,
    junx_rates[match(mut$sra, rownames(junx_rates)),]
  )
  
  aframe$label <- paste(
    aframe$ProteinChange,
    sample_info[aframe$ModelID, "stripped_cell_line_name"])
  
  aframe
}
aframe_4 <- clean_up(mut[mut$HugoSymbol == "TBRG4", ])
aframe_5 <- clean_up(mut[mut$HugoSymbol == "FASTKD5", ])

plot_heatmap <- function(aframe){
  tmp <- junx_rates[match(aframe$sra, rownames(junx_rates)),]
  
  anno_row <- aframe
  rownames(tmp) <- rownames(anno_row) <- anno_row$label
  
  tmp <- tmp[order(anno_row$prot_pos), ]
  
  pheatmap(tmp, 
           scale = "column", breaks = seq(-2, 2, length = 100),
           annotation_row = anno_row[, c("prot_pos", "VariantInfo")],
           labels_row = anno_row[rownames(tmp), "label"],
           labels_col = anno$gene,
           cluster_cols = F, cluster_rows = F)
}
plot_heatmap(aframe_4)
plot_heatmap(aframe_5)

plot_junction <- function(aframe, junx = "MT.TT"){
  colnames(aframe) <- gsub(junx, "junx", colnames(aframe))
  
  ggplot(aframe, aes(prot_pos, junx,
                     color = VariantInfo,
                     label = label)) +
    labs(x = "Position in protein",
         y = paste(junx, "cleavage")) +
    geom_point() +
    ggrepel::geom_text_repel(size = 2) +
    #geom_vline(xintercept = c(369, 437), color = 'green', linetype = 'dashed') +
    #geom_vline(xintercept = c(450, 536), color = 'red', linetype = 'dashed') +
    #geom_vline(xintercept = c(569, 620), color = 'blue', linetype = 'dashed') +
    theme_classic()  
}

lapply(colnames(junx_rates), function(x){
  plot_junction(aframe = aframe_4, x)
  ggsave(paste0("/Users/lukas/Downloads/tbrg4_", 
                gsub(".", "_", fixed = T, x), ".png"))
})
lapply(colnames(junx_rates), function(x){
  plot_junction(aframe = aframe_5, x)
  ggsave(paste0("/Users/lukas/Downloads/fastkd5_", 
                gsub(".", "_", fixed = T, x), ".png"))
})



# Correlate w gene dependencies ####
cells <- sample_info$DepMap_ID[match(sra$Sample.Name, sample_info$CCLE_Name)]
nas <- apply(demeter2, 2, function(x) sum(is.na(x)))
correl <- cor(demeter2[match(cells, rownames(demeter2)), which(nas < 10)],
              junx_rates, use = 'pairwise.complete')
correl <- cor(kronos[match(cells, rownames(kronos)), ],
              junx_rates, use = 'pairwise.complete')

maxs <- apply(correl, 1, function(x) max(abs(x)))
ok <- names(tail(sort(maxs), 50))

pheatmap(correl[ok, ], 
         cluster_cols = FALSE,
         main = "Correlation w Chronos dependencies")

aframe <- data.frame(
  kronos[match(cells, rownames(kronos)), "FASTKD5"]
)

