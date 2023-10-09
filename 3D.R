library(prim)
source("prim_edits.R")
setwd("Path_to_the_experimental_data")

qual <- rep(0,5)

for(k in 1:5){
  d <- read.table(paste0("data_aggr_dimensions=", k, "_multi_true"), sep=',')

  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0, k),rep(0.2, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.3, k),rep(0.5, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.6, k),rep(0.8, k)), ncol = k, byrow = TRUE))
  mass.min = 0.2^k
  thr = 3
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05)
  qual[k] <- evaluate(res$box[res$y.fun > thr*0], gt)
}
# check evaluation for k=1


qual <- rep(0,5)

for(k in 1:4){
  print(k)
  d <- read.table(paste0("data_density_dimensions=", k, "_multi_True"), sep=',')
  
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
  
  avgs <- numeric()
  for(i in 0:2){
    inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.3*i < row & row < 0.3*i+0.2)))
    avgs <- c(avgs, mean(d[inds, k + 1]))
  }
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0, k),rep(0.2, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.3, k),rep(0.5, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.6, k),rep(0.8, k)), ncol = k, byrow = TRUE))
  mass.min = 0.2^k
  thr = 0.66*min(avgs)
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05, verbose = TRUE)
  qual[k] <- evaluate(res$box[res$y.fun > thr], gt)
}

# 0.3404283 0.8132877 0.6017588 0.6719577 0.0000000

qual <- rep(0,5)
library(FNN)

for(k in 1:5){
  print(k)
  d <- read.table(paste0("data_density_dimensions=", k, "_multi_True"), sep=',')

  d.tmp <- as.data.frame(matrix(runif(100000*k), ncol = k))
  d.tmp[, k + 1] <- -get.knnx(d, d.tmp, k = 1)$nn.dist[, 1]
  d <- d.tmp
  
  avgs <- numeric()
  for(i in 0:2){
    inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.3*i < row & row < 0.3*i+0.2)))
    avgs <- c(avgs, mean(d[inds, k + 1]))
  }
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0, k),rep(0.2, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.3, k),rep(0.5, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.6, k),rep(0.8, k)), ncol = k, byrow = TRUE))
  mass.min = 0.2^k
  thr = 1.5*min(avgs)
  res <- prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
                  threshold.type = 1, mass.min = mass.min, 
                  peel.alpha = 0.05, verbose = TRUE)
  qual[k] <- evaluate(res$box[res$y.fun > thr], gt)
}

# 0.9781610 0.9118547 0.9086797 0.8636826 0.7006626 # knn = 5
# 0.9800453 0.9432002 0.9196315 0.9009851 0.8112448 # knn = 1


