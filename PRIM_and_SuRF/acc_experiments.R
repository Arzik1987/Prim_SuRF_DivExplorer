
library(prim)
source("prim_edits.R")

#### Statistic 'Aggregate', 1 GT region

for(k in 1:5){
  print(k)
  d <- read.table(paste0("input\\data_aggr_dimensions=", k, "_multi_False"), sep=',')
  inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.6 < row & row < 0.9)))
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0.6, k),rep(0.9, k)), ncol = k, byrow = TRUE))
  mass.min = 0.3^k
  thr = 2
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_boxes\\data_aggr_dimensions=", k, "_multi_False.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_boxes\\data_aggr_dimensions=", k, "_multi_False.csv"))
  } 
}

#### Statistic 'Aggregate', 3 GT regions

for(k in 1:5){
  print(k)
  d <- read.table(paste0("input\\data_aggr_dimensions=", k, "_multi_True"), sep=',')
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0, k),rep(0.2, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.3, k),rep(0.5, k)), ncol = k, byrow = TRUE),
             matrix(c(rep(0.6, k),rep(0.8, k)), ncol = k, byrow = TRUE))
  mass.min = 0.2^k
  thr = 2
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_boxes\\data_aggr_dimensions=", k, "_multi_True.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_boxes\\data_aggr_dimensions=", k, "_multi_True.csv"))
  } 
}


#### Statistic 'Density', 1 GT region

library(FNN)

for(k in 1:5){
  print(k)
  d <- read.table(paste0("input\\data_density_dimensions=", k, "_multi_False"), sep=',')
  
  d.tmp <- as.data.frame(matrix(runif(100000*k), ncol = k))
  d.tmp[, k + 1] <- -get.knnx(d, d.tmp, k = 1)$nn.dist[, 1]
  d <- d.tmp
  
  inds <- which(apply(d[, 1:k, drop=FALSE], 1, function(row) all(0.6 < row & row < 0.9)))
  
  box <- matrix(c(rep(0, k),rep(1, k)), ncol = k, byrow = TRUE)
  gt <- list(matrix(c(rep(0.6, k),rep(0.9, k)), ncol = k, byrow = TRUE))
  mass.min = 0.3^k
  thr = 1.1*mean(d[inds, k + 1])
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_boxes\\data_density_dimensions=", k, "_multi_False.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_boxes\\data_density_dimensions=", k, "_multi_False.csv"))
  }
}

#### Statistic 'Density', 3 GT regions

for(k in 1:5){
  print(k)
  d <- read.table(paste0("input\\data_density_dimensions=", k, "_multi_True"), sep=',')
  
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
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_boxes\\data_density_dimensions=", k, "_multi_True.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_boxes\\data_density_dimensions=", k, "_multi_True.csv"))
  }
}


#############
#### KDE ####
#############

#### NOTE: For 5 dimensions it takes some time to estimate kde. So we just store 
#### the initial box, as we use kde for exclusively for demonstration

#### Statistic 'Density', 1 GT region

for(k in 1:4){
  print(k)
  d <- read.table(paste0("input\\data_density_dimensions=", k, "_multi_False"), sep=',')

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
  thr = 0.9*mean(d[inds, k + 1])
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_kde\\data_density_dimensions=", k, "_multi_False.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_kde\\data_density_dimensions=", k, "_multi_False.csv"))
  }
}

#### Statistic 'Density', 3 GT regions

library(ks)

for(k in 1:4){
  print(k)
  d <- read.table(paste0("input\\data_density_dimensions=", k, "_multi_True"), sep=',')
  
  if (k == 1){
    kde_result = density(d[,1])
    d <- as.data.frame(matrix(runif(100000*k), ncol = k))
    d[,2] = approx(x = kde_result$x, y = kde_result$y, xout = d[,1])$y
  } else {
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
  
  res <- tryCatch({
    prim.box(d[, 1:k, drop = FALSE], d[,k + 1], box.init = box, threshold = thr, 
             threshold.type = 1, mass.min = mass.min, 
             peel.alpha = 0.05)
  }, error=function(e) {
    NA
  })
  
  if(!is.list(res)){
    save_res_to_file(list(box), paste0("prim_kde\\data_density_dimensions=", k, "_multi_True.csv"))
  } else {
    save_res_to_file(res$box, paste0("prim_kde\\data_density_dimensions=", k, "_multi_True.csv"))
  }
}

#######################################
#### PRIM modification for density ####
#######################################

for(k in 1:5){
  print(k)
  d <- read.table(paste0("input\\data_density_dimensions=", k, "_multi_False"), sep=',')
  res <- PRIMdens(d)
  mass.min = 0.3^k
  best_box <- res$boxes[max(which(sapply(res$boxes, function(x) prod(x[2,]-x[1,])) >= mass.min))]
  save_res_to_file(best_box, paste0("prim_dens\\data_density_dimensions=", k, "_multi_False.csv"))
}


