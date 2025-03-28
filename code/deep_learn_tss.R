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
  x$num_of_starts/(x$total_num_reads + 1)))
end_rate <- do.call(cbind, lapply(dat, function(x)
  x$num_of_ends/(x$total_num_reads + 1)))

# Build model ####
library(keras)

n_bases <- nrow(coverage)
n_width <- 100

y_anno <- rep(0, n_bases)
y_anno[anno$V4] <- 1

model <- keras_model_sequential()

model %>%
  layer_conv_1d(
    filter = 16, kernel_size = c(5), padding = "same",
    input_shape = c(n_width, 2)
  ) %>%
  layer_conv_1d(
    filter = 8, kernel_size = c(5), padding = "same",
    input_shape = c(n_width, 2)
  ) %>%
  layer_conv_1d(2, 1, activation = "sigmoid")

model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam',
  metrics = 'accuracy'
)

training <- lapply(1:10, 
                   function(x){
                     
  pairs <- data.frame(
    pos = sample(n_bases - n_width, 1000000, replace = T),
    sample = sample(ncol(coverage), 1000000, replace = T))
                     
  r <- lapply(pairs$pos, function(x) c(x:(x + (n_width - 1))))
  
  #n_sample <- 1
  
  X <- unlist(lapply(1:length(r), function(x)
    rbind(
      start_rate[r[[x]], pairs$sample[x]],
      end_rate[r[[x]], pairs$sample[x]])))
  X <- array_reshape(X, c(length(r), n_width, 2))
  
  y <- do.call(rbind, lapply(r, function(x) y_anno[x]))
  ok <- which(rowSums(y) > 0)
  
  y <- to_categorical(y[ok,], 2)
  X <- X[ok,,]
  
  model %>% fit(X, y, epochs = 10, batch_size = 128)
  
  model$history
})

plot(X[1,,1], col = "red", ylim = c(0, 0.4))
points(X[1,,2], col = "blue")
abline(v = which(y[3,] == 1))

pos <- 577
range <- (pos - 50):(pos + 49)
X_train <- rbind(
  start_rate[range, 1], end_rate[range, 1])
X_train <- array_reshape(X_train, c(1, 100, 2))
preds <- predict(model, X_train)
plot(preds[1,,2])
