

dir.create(file.path(getwd(), "prim_boxes"), showWarnings = FALSE)
dir.create(file.path(getwd(), "prim_kde"), showWarnings = FALSE)
dir.create(file.path(getwd(), "prim_dens"), showWarnings = FALSE)


save_res_to_file <- function(boxes, filename) {
  for(i in 1:length(boxes)){
    boxes[[i]] <- c(boxes[[i]][1,], boxes[[i]][2,]-boxes[[i]][1,])
  }
  mat <- do.call(rbind, boxes)
  write.table(mat, filename, row.names = FALSE, col.names = FALSE, sep = ",")
}

### to work with 1D data, we do a small change in the internal function of 
### peel package: We replace "x.new <- x[x.index, ]" with 
### "x.new <- x[x.index, , drop = FALSE]"

assignInNamespace("peel.one",
                  value = function (x, y, box, peel.alpha, mass.min, threshold, d, n, type = 8, 
                             y.fun = mean, verbose){
  box.new <- box
  mass <- length(y)/n
  if(is.vector(x)) 
    return(NULL)
  
  y.fun.val <- do.call(y.fun, list(x = y))
  y.fun.peel <- matrix(0, nrow = 2, ncol = d)
  box.vol.peel <- matrix(0, nrow = 2, ncol = d)
  for(j in 1:d){
    box.min.new <- quantile(x[, j], peel.alpha, type = type)
    box.max.new <- quantile(x[, j], 1 - peel.alpha, type = type)
    y.fun.peel[1, j] <- do.call(y.fun, list(x = y[x[, j] >= box.min.new]))
    y.fun.peel[2, j] <- do.call(y.fun, list(x = y[x[, j] <= box.max.new]))
    box.temp1 <- box
    box.temp2 <- box
    box.temp1[1, j] <- box.min.new
    box.temp2[2, j] <- box.max.new
    box.vol.peel[1, j] <- prim:::vol.box(box.temp1)
    box.vol.peel[2, j] <- prim:::vol.box(box.temp2)
  }
  y.fun.peel.max.ind <- which(y.fun.peel == max(y.fun.peel, na.rm = TRUE), arr.ind = TRUE)
  nrr <- nrow(y.fun.peel.max.ind)
  if(nrr > 1){
    box.vol.peel2 <- rep(0, nrr)
    for (j in 1:nrr) box.vol.peel2[j] <- box.vol.peel[y.fun.peel.max.ind[j,1], y.fun.peel.max.ind[j, 2]]
    row.ind <- which(max(box.vol.peel2) == box.vol.peel2)
  } else row.ind <- 1
  
  y.fun.peel.max.ind <- y.fun.peel.max.ind[row.ind, ]
  j.max <- y.fun.peel.max.ind[2]
  
  if (y.fun.peel.max.ind[1] == 1) {
    box.new[1, j.max] <- quantile(x[, j.max], peel.alpha, type = type)
    x.index <- x[, j.max] >= box.new[1, j.max]
  }
  else if(y.fun.peel.max.ind[1] == 2){
    box.new[2, j.max] <- quantile(x[, j.max], 1 - peel.alpha, type = type)
    x.index <- x[, j.max] <= box.new[2, j.max]
  }
  if(verbose){
    cat("Peeled in dimension", j.max, ": new limits are", 
        paste("[", signif(box.new[1, j.max], 4), ",", 
              signif(box.new[2, j.max], 4), "]\n", sep = ""))
  }
    
  x.new <- x[x.index, , drop = FALSE]
  y.new <- y[x.index]
  mass.new <- length(y.new)/n
  y.fun.new <- do.call(y.fun, list(x = y.new))
  if((y.fun.new >= threshold) & (mass.new >= mass.min) & (mass.new < mass)){
    return(list(x = x.new, y = y.new, y.fun = y.fun.new, 
                box = box.new, mass = mass.new))
  }
}, ns = "prim")


### there is no need to discover thousands of boxes when we set the minimal
### size of a box to be small. So we edited the "prim.one" function by
### restricting their number as "num.boxes <- min(k.max, 5)"

assignInNamespace("prim.one",
                  value = function (x, y, box.init = NULL, peel.alpha = 0.05, paste.alpha = 0.01, 
          mass.min = 0.05, threshold, pasting = FALSE, threshold.type = 1, 
          verbose = FALSE, y.fun = mean) 
{
  d <- ncol(x)
  n <- nrow(x)
  k.max <- ceiling(1/mass.min)
  num.boxes <- min(k.max, 5)
  y.fun.val <- do.call(y.fun, list(x = y))
  mass.init <- length(y)/n
  if (is.null(box.init)) {
    box.init <- apply(x, 2, range)
    box.diff <- box.init[2, ] - box.init[1, ]
    box.init[1, ] <- box.init[1, ] - 10 * paste.alpha * box.diff
    box.init[2, ] <- box.init[2, ] + 10 * paste.alpha * box.diff
  }
  k <- 1
  boxk <- prim:::find.box(x = x, y = y, box = box.init, peel.alpha = peel.alpha, 
                   paste.alpha = paste.alpha, mass.min = mass.min, threshold = min(y) - 
                     0.1 * abs(min(y)), d = d, n = n, pasting = pasting, 
                   verbose = verbose, y.fun = y.fun)
  if (is.null(boxk)) {
    if (verbose) 
      warning(paste("Unable to find box", k, "\n"))
    x.prim <- list(x = list(x), y = list(threshold.type * 
                                           y), y.fun = threshold.type * y.fun.val, box = list(box.init), 
                   box.mass = mass.init, num.class = 1, num.hdr.class = 1, 
                   threshold = do.call(y.fun, list(x = y)))
    class(x.prim) <- "prim"
    return(x.prim)
  } else {
    if (verbose) 
      cat(paste("Found box ", k, ": y.fun=", 
                signif(threshold.type * boxk$y.fun, 4), ", mass=", 
                signif(boxk$mass, 4), "\n\n", sep = ""))
    boxes <- list(x = list(boxk$x), y = list(boxk$y), y.fun = list(boxk$y.fun), 
                  box = list(boxk$box), mass = list(boxk$mass))
  }
  if (num.boxes > 1) {
    boxk <- list(x = boxes$x[[k]], y = boxes$y[[k]], y.fun = boxes$y.fun[[k]], 
                 box = boxes$box[[k]], mass = boxes$mass[[k]])
    x.out.ind.mat <- matrix(TRUE, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:d) x.out.ind.mat[, j] <- (x[, j] < boxk$box[1, 
                                                            j]) | (x[, j] > boxk$box[2, j])
    x.out.ind <- apply(x.out.ind.mat, 1, sum) != 0
    x.out <- x[x.out.ind, ]
    if (is.vector(x.out)) 
      x.out <- as.matrix(t(x.out))
    y.out <- y[x.out.ind]
    while ((length(y.out) > 0) & (k < num.boxes) & (!is.null(boxk))) {
      k <- k + 1
      boxk <- prim:::find.box(x = x.out, y = y.out, box = box.init, 
                       peel.alpha = peel.alpha, paste.alpha = paste.alpha, 
                       mass.min = mass.min, threshold = min(y) - 0.1 * 
                         abs(min(y)), d = d, n = n, pasting = pasting, 
                       verbose = verbose, y.fun = y.fun)
      if (is.null(boxk)) {
        if (verbose) 
          cat(paste("Bump", k, "includes all remaining data\n\n"))
        boxes$x[[k]] <- x.out
        boxes$y[[k]] <- y.out
        boxes$y.fun[[k]] <- do.call(y.fun, list(x = y.out))
        boxes$box[[k]] <- box.init
        boxes$mass[[k]] <- length(y.out)/n
      } else {
        if (verbose) 
          cat(paste("Found box ", k, ": y.fun=", 
                    signif(threshold.type * boxk$y.fun, 4), ", mass=", 
                    signif(boxk$mass, 4), "\n\n", sep = ""))
        x.out.ind.mat <- matrix(TRUE, nrow = nrow(x), 
                                ncol = ncol(x))
        for (j in 1:d) x.out.ind.mat[, j] <- (x[, j] < 
                                                boxk$box[1, j]) | (x[, j] > boxk$box[2, j])
        x.out.ind <- x.out.ind & (apply(x.out.ind.mat, 
                                        1, sum) != 0)
        x.out <- x[x.out.ind, ]
        if (is.vector(x.out)) 
          x.out <- as.matrix(t(x.out))
        y.out <- y[x.out.ind]
        boxes$x[[k]] <- boxk$x
        boxes$y[[k]] <- boxk$y
        boxes$y.fun[[k]] <- boxk$y.fun
        boxes$box[[k]] <- boxk$box
        boxes$mass[[k]] <- boxk$mass
      }
    }
  }
  for (k in 1:length(boxes$y.fun)) {
    boxes$y[[k]] <- threshold.type * boxes$y[[k]]
    boxes$y.fun[[k]] <- threshold.type * boxes$y.fun[[k]]
  }
  prim.res <- prim:::prim.hdr(prim = boxes, threshold = threshold, 
                       threshold.type = threshold.type, y.fun = y.fun)
  return(prim.res)
}, ns = "prim")



PRIMdens <- function(X, alpha=0.05) {
  # Initialize
  boxes <- list()
  densities <- numeric()
  
  n_points <- nrow(X)
  
  # Initial box is the unit box
  box <- rbind(rep(0, ncol(X)), rep(1, ncol(X)))
  
  for(iteration in 1:120) {
    max_density <- -Inf
    best_cut <- NULL
    best_box <- NULL
    n_in_best_box <- n_points
    ind_in_best_box <- NULL
    
    for(dim in 1:ncol(X)) {
      for(direction in c(0, 1)) { # 0 for lower bound, 1 for upper bound
        trial_box <- box
        current_dim_length <- trial_box[2, dim] - trial_box[1, dim]
        
        # Adjust the boundary in the dimension
        adjustment <- current_dim_length * alpha
        if(direction == 0) {
          trial_box[1, dim] <- trial_box[1, dim] + adjustment
          in_box <- X[, dim] >= trial_box[1, dim]
        } else {
          trial_box[2, dim] <- trial_box[2, dim] - adjustment
          in_box <- X[, dim] <= trial_box[2, dim]
        }
        
        # Count points in the trial box
        n_in_box <- sum(in_box)
        
        # Calculate volume of the trial box
        volume <- prod(trial_box[2, ] - trial_box[1, ])
        
        # Calculate density
        density <- ifelse(volume > 0, n_in_box / volume, 0)
        
        if(density > max_density) {
          max_density <- density
          best_cut <- list(dim, direction)
          best_box <- trial_box
          n_in_best_box <- n_in_box
          ind_in_best_box <- in_box
        }
      }
    }
    
    # Break if no cut is found or too few points remain
    if(is.null(best_cut) || n_in_best_box < 1/alpha) break
    
    # Update the box and record
    box <- best_box
    boxes <- c(boxes, list(box))
    densities <- c(densities, max_density)
    X <- X[ind_in_best_box, , drop=FALSE]
  }
  
  return(list(boxes=boxes, densities=densities))
}



### Evaluation boxes

# iou <- function(box1, box2){
#   
#   inter <- box1
#   inter[1, ] <- apply(rbind(box1[1, ], box2[1, ]), 2, max)
#   inter[2, ] <- apply(rbind(box1[2, ], box2[2, ]), 2, min)
#   
#   sides1 <- box1[2, ] - box1[1, ]
#   sides2 <- box2[2, ] - box2[1, ]
#   sides.inter <- inter[2, ] - inter[1, ]
#   
#   if(sum(sides.inter<=0) > 0){
#     iou = 0
#   } else {
#     iou = 1/(prod(sides1/sides.inter) + prod(sides2/sides.inter) - 1)
#   }
#   iou
# }

# box1 <- matrix(c(0,0,0.3,0.3), ncol = 2, byrow = TRUE)
# box2 <- matrix(c(0.2,0.2,0.5,0.5), ncol = 2, byrow = TRUE)  
# box3 <- matrix(c(0.7,0.7,1,1), ncol = 2, byrow = TRUE)
# box4 <- matrix(c(0,0,1,1), ncol = 2, byrow = TRUE)

# iou(box1,box2) # 0.05882353
# iou(box2,box1) # 0.05882353
# iou(box2,box3) # 0
# iou(box4,box1) # 0.09


# evaluate <- function(boxes, gt){
#   if(length(boxes) != length(gt)){
#     warning("The number ob boxes discovered is not equal to the true number of boxes")
#   }
#   
#   k <- length(boxes)
#   l <- length(gt)
#   res <- matrix(rep(NA, k*l), ncol = l)
#   for(i in 1:k){
#     for(j in 1:l){
#       res[i, j] <- iou(boxes[[i]], gt[[j]])
#     }
#   }
#   mean(apply(res, 2, max))
# }








