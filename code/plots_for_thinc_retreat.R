# Load R libs ####
library(biomaRt)
library(data.table)
library(fgsea)
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

global_cleave <- apply(gene_rates, 1, median)

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
genes <- c("FASTK", "FASTKD1", "FASTKD2", 
           "FASTKD3", "FASTKD5", "TBRG4",
           "TRMT10C", "MRPP1", 
           "POLRMT", "TFAM", "TFB2M", "TEFM")
res$label[rownames(res) %in% genes] <- rownames(res)[rownames(res) %in% genes]

ggplot(res,
       aes(rank, score, color = chromosome_name, label = label)) +
  geom_point() +
  labs(title = "CNV association with global cleavage",
       y = "Association score",
       x = "Genomic positions") +
  geom_label(size = 5, color = "black") +
  theme_classic()


# Plot TRBG4 CNV vs MT cleavage ####
aframe <- data.frame(
  sample_info,
  cnv = cnv[match(rownames(sample_info), rownames(cnv)), "TBRG4"]
)
aframe$group <- "WT"
aframe$group[aframe$cnv > 1.2] <- "gain"
aframe$group[aframe$cnv < 0.8] <- "loss"
aframe$group <- factor(aframe$group, levels = c("loss", "WT", "gain"))

ggplot(aframe, aes(
  group, cleave)) +
  labs(
    title = "DepMap",
    subtitle = "TBRG4 copy number vs MT cleavage",
    y = "MT cleavage rate",
    x = "TBRG4 copy number") +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_classic()


# Correlate RNA expression w cleavage in WT ####
ok <- aframe$DepMap_ID[which(aframe$group == "WT" & 
                               !is.na(aframe$cleave))]
ok <- intersect(ok, rownames(rnaseq))

res_rna <- apply(rnaseq[match(ok, rownames(rnaseq)), ], 2, function(x){
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


# Run enrichment ####
hallmark <- gmtPathways("/Users/lukas/OneDrive/Miko/UTHealth/projects/Marian/data/hallmark.gmt")
ranks <- -log10(res$pval)
ranks[which(res$coef < 0)] <- (-1)*ranks[which(res$coef < 0)]
names(ranks) <- res$hgnc_symbol
enrich <- fgsea(hallmark, stats = ranks)


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

correl_chr <- t(apply(kronos[ok,], 2, function(x)
  coefficients(summary(lm(
    x ~ sample_info$cleave[match(ok, rownames(sample_info))])))[2, c(1, 4)]))

ok <- intersect(
  rownames(demeter2),
  sample_info$DepMap_ID[!is.na(sample_info$cleave)])

correl_dem <- t(apply(demeter2[ok,], 2, function(x)
  coefficients(summary(lm(
    x ~ sample_info$cleave[match(ok, rownames(sample_info))])))[2, c(1, 4)]))

colnames(correl_chr) <- colnames(correl_dem) <- c("coef", "pval")


# Plot MT cleavage by cancer type ####
good_tissues <- tail(names(sort(table(sample_info$lineage))), 10)
sample_info$group <- "rest"
sample_info$group[sample_info$lineage %in% good_tissues] <- sample_info$lineage[sample_info$lineage %in% good_tissues]

ggplot(sample_info, aes(cleave, group)) +
  geom_boxplot() +
  theme_classic()


# Volcano plots for dependency ####
aframe <- data.frame(
  gene = rownames(correl_chr),
  correl_chr
)
aframe <- aframe[order(aframe$pval), ]
aframe$group <- aframe$gene %in% hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION
aframe <- aframe[order(aframe$group), ]
ggplot(aframe, aes(coef, -log10(pval),
                   color = group)) +
  labs(
    title = "DepMap",
    subtitle = "Chronos"
  ) +
  geom_point() +
  ggrepel::geom_text_repel(
    data = aframe[1:20,],
    aes(label = gene), color = "red") +
   geom_text(
     data = aframe[aframe$gene %in% c("COX5B"),],
     aes(label = gene), color = "blue") +
  theme_classic()

aframe <- data.frame(
  gene = rownames(correl_dem),
  correl_dem
)
aframe <- aframe[order(aframe$pval), ]
aframe$group <- aframe$gene %in% hallmark$HALLMARK_G2M_CHECKPOINT
aframe <- aframe[order(aframe$group), ]
ggplot(aframe, aes(coef, -log10(pval),
                   color = group)) +
  labs(
    title = "DepMap",
    subtitle = "D2"
  ) +
  geom_point() +
  ggrepel::geom_text_repel(
    data = aframe[1:20,],
    aes(label = gene), color = "red") +
  # ggrepel::geom_text_repel(
  #   data = aframe[aframe$gene %in% genes,],
  #   aes(label = gene), color = "blue") +
  theme_classic()

plot_gene <- function(gene, type = "d2"){
  if(type == "d2") matr <- demeter2
  if(type != "d2") matr <- kronos
  
  ok <- intersect(
    rownames(matr),
    rownames(sample_info))
  subm <- data.frame(
    chr = matr[ok, gene],
    cnv = cnv[match(ok, rownames(cnv)), gene],
    sample_info[match(ok, rownames(sample_info)), ])
  ggplot(subm, aes(cleave, chr)) +
    labs(
      y = paste(gene, "dependency"),
      x = "Global cleavage rate"
    ) +
    geom_point() +
    ggpubr::stat_cor() +
    theme_classic()
  
}
plot_gene("COX5B", "chr")

plot_gene_box <- function(gene, type = "d2"){
  if(type == "d2") matr <- demeter2
  if(type != "d2") matr <- kronos
  
  ok <- intersect(
    rownames(matr),
    rownames(sample_info))
  subm <- data.frame(
    chr = matr[ok, gene],
    cnv = cnv[match(ok, rownames(cnv)), gene],
    sample_info[match(ok, rownames(sample_info)), ])
  subm <- subm[!is.na(subm$cleave) & !is.na(subm$chr), ]
  
  subm$group <- cut(
    subm$cleave,
    breaks = quantile(subm$cleave, seq(0, 1, length = 4)),
    labels = c("low", "medium", "high"))
  
  ggplot(subm[subm$group %in% c("low", "high"),],
         aes(group, chr)) +
    labs(
      title = "DepMap",
      y = paste(gene, "dependency [Chronos]"),
      x = "Global cleavage rate"
    ) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.1) +
    ggpubr::stat_compare_means() +
    theme_classic()
  
}
plot_gene_box(gene = "COX5B", type = "chr")


# Run enrichment ####
hallmark <- gmtPathways("/Users/lukas/OneDrive/Miko/UTHealth/projects/Marian/data/hallmark.gmt")
ranks <- correl_chr[,1]
names(ranks) <- rownames(correl_chr)
enrich_chr <- fgsea(hallmark, stats = ranks)

tmp <- data.frame(enrich_chr)
tmp$pathway <- gsub("HALLMARK_", "", tmp$pathway)
tmp$pathway <- factor(
  tmp$pathway, levels = tmp$pathway[order(tmp$NES)])
ggplot(tmp) +
  labs(
    title = "Chronos",
    subtitle = "Hallmark pathways",
    x = "Enrichment score",
    y = ""
  ) +
  geom_bar(
    aes(y = pathway, x = NES, fill = NES),
    stat = "identity") +
  scale_fill_gradient2(
    low = "blue", mid = "grey", high = "red") +
  theme_classic()

# Run enrichment ####
ranks <- correl_dem[,1]
names(ranks) <- rownames(correl_dem)
enrich_d2 <- fgsea(hallmark, stats = ranks)


# Correlate with drug sensitivity ####
auc <- read.csv(
  "/Users/lukas/OneDrive/Documents/GitHub/HNRNPH1-in-MCL/data/depmap/24Q2/depmap_drug_auc.csv",
  row.names = 1)
drug <- read.csv(
  "/Users/lukas/OneDrive/Documents/GitHub/HNRNPH1-in-MCL/data/depmap/24Q2/depmap_drug_metadata.csv",
  row.names = 1)

ok <- intersect(rownames(auc), rownames(sample_info))
res_drug <- t(apply(auc[ok, ], 2, function(x)
  summary(lm(
    x ~ sample_info[ok, "cleave"]))$coefficients[2, c(1, 4)]))
res_drug <- data.frame(
  drug,
  res_drug)


# Load metabolism data ####
meta <- read.csv("/Users/lukas/Downloads/CCLE_metabolomics_20190502.csv")
metabolome <- data.frame(meta[, -c(1:2)])
ok <- which(!is.na(meta$DepMap_ID))
metabolome <- metabolome[ok, ]
rownames(metabolome) <- meta$DepMap_ID[ok]

ok <- intersect(
  rownames(metabolome),
  rownames(sample_info))

res_meta <- t(apply(metabolome[ok, ], 2, function(x)
  coefficients(summary(lm(
    #x ~ sample_info[ok, "lineage"] + sample_info[ok, "cleave"])))[3, c(1, 4)]))
    x ~ sample_info[ok, "cleave"])))[2, c(1, 4)]))
res_meta <- data.frame(res_meta)
colnames(res_meta) <- c("coef", "pval")
res_meta$metabolite <- rownames(res_meta)

res_meta <- res_meta[order(res_meta$pval), ]
ggplot(res_meta, aes(coef, -log10(pval))) +
  labs(
    title = "DepMap",
    subtitle = "Association with global MT cleavage"
  ) +
  geom_point() +
  ggrepel::geom_text_repel(
    data = res_meta[1:20,],
    aes(label = metabolite), color = "red") +
  theme_classic()

subm <- data.frame(
  xan = metabolome[ok, "oxalate"],
  sample_info[ok,]
)
summary(lm(
  xan ~ lineage + cleave,
  data = subm))

ggplot(subm, aes(xan, cleave)) +
  labs(
    title = "DepMap",
    x = "Xanthine levels",
    y = "MT cleavage"
  ) +
  stat_cor() +
  geom_point() +
  theme_classic()
