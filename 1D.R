library(prim)
source("prim_edits.R")
setwd("Path_to_the_experimental_data")

qual <- rep(0,5)

for(k in 1:5){
  d <- read.table(paste0("data_aggr_dimensions=", k, "_multi_False"), sep=',')
  inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.6 < row & row < 0.9)))
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0.6, k),rep(0.9, k)), ncol = k, byrow = TRUE))
  mass.min = 0.3^k
  # if(length(inds)/nrow(d) < 0.9) mass.min = length(inds)/nrow(d)
  thr = 0.9*mean(d[inds, k + 1])
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05)
  qual[k] <- evaluate(res$box[res$y.fun > thr], gt)
}


qual <- rep(0,5)

for(k in 1:4){
  print(k)
  d <- read.table(paste0("data_density_dimensions=", k, "_multi_False"), sep=',')
  
  if (k == 1){
    kde_result = density(d[,1])
    d <- as.data.frame(matrix(runif(100000*k), ncol = k))
    d[,2] = approx(x = kde_result$x, y = kde_result$y, xout = d[,1])$y
  } else {
    library(ks)
    kde_result <- kde(d)
    d <- as.data.frame(matrix(runif(100000*k), ncol = k))
    d[, k + 1] <- predict(kde_result, x = d)
  }
  
  inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.6 < row & row < 0.9)))
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0.6, k),rep(0.9, k)), ncol = k, byrow = TRUE))
  mass.min = 0.3^k
  # if(length(inds)/nrow(d) < 0.9) mass.min = length(inds)/nrow(d)
  thr = 0.9*mean(d[inds, k + 1])
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05, verbose = TRUE)
  qual[k] <- evaluate(res$box[res$y.fun > thr], gt)
}

# 0.9750076 0.9080086 0.9342720 0.9345526 0.0000000

qual <- rep(0,5)

for(k in 1:5){
  print(k)
  d <- read.table(paste0("data_density_dimensions=", k, "_multi_False"), sep=',')

  d.tmp <- as.data.frame(matrix(runif(100000*k), ncol = k))
  d.tmp[, k + 1] <- -get.knnx(d, d.tmp, k = 1)$nn.dist[, 1]
  d <- d.tmp
  
  inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.6 < row & row < 0.9)))
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0.6, k),rep(0.9, k)), ncol = k, byrow = TRUE))
  mass.min = 0.3^k
  thr = 1.1*mean(d[inds, k + 1])
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05, verbose = TRUE)
  qual[k] <- evaluate(res$box[res$y.fun > thr], gt)
}

# 0.9781610 0.9118547 0.9086797 0.8636826 0.7006626 # knn = 5
# 0.9800453 0.9432002 0.9196315 0.9009851 0.8112448 # knn = 1


