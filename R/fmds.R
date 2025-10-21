#' Distance between vectors
#'
#' @param z Matrix or vector of observations
#'
#' @return Distance matrix
#' @importFrom stats dist
#' @export
#' @examples
#' set.seed(100)
#' z <- rbind(matrix(rnorm(100), ncol=4),
#' matrix(rnorm(100, 2), ncol=4))
#' getDistMat(z)
getDistMat <- function(z){
    z_dist <- as.matrix(dist(z)) # L2 distance
    return(z_dist)
}


#' Get p-value matrix
#'
#' @param z Lower dimension representation
#' @param D Original distance matrix
#' @param y Treatment vector
#'
#' @return pseudo-F values matrix
#' 1st col of original data, 2nd col of reduced dim
#' @importFrom stats cmdscale
#' @importFrom phyloseq distance sample_data
#' @export
#' @examples
#' require(phyloseq)
#' data(microbiome)
#' D <- distance(microbiome, method = 'wunifrac') # requires phyloseq package
#' y <- sample_data(microbiome)$Treatment
#' z0 <- cmdscale(d = D)
#' pairByRank(z = z0, D = D, y = y)
pairByRank <- function(z, D, y){
    f0_sorted <- getP(D = D, y = y)$ratio_all
    fz_sorted <- getP(z = z, y = y)$ratio_all
    mat_pair <- cbind(f0_sorted, fz_sorted)
    return(mat_pair)
}


#' Objective term of MDS
#'
#' @param D Original distance matrix
#' @param z Lower dimension representation
#' @param N Number of observations--scaling factors
#'
#' @return Scalar of objective function value of MDS
#' @importFrom stats dist cmdscale
#' @importFrom phyloseq distance
#' @export
#' @examples
#' require(phyloseq)
#' data(microbiome)
#' D <- distance(microbiome, method = 'wunifrac') # requires phyloseq package
#' z0 <- cmdscale(d = D)
#' N <- dim(z0)[1]
#' mdsObj(D = D, z = z0, N = N)
mdsObj <- function(D, z, N){
    z_distmat <- getDistMat(z)
    return(sum((D - z_distmat)^2) / N)
}



#' FMDS calculation using MM algorithm
#'
#' @param D Square matrix of pairwise distance, size of N by N
#' @param y Vector of label or group set, size of N
#' @param X Object matrix; used to build distance matrix D; D is prioritized
#' @param nit Number of iterations; 100 by default
#' @param lambda Hyperparameter; 0.5 by default
#' @param threshold_p Lower limit of p-value difference that allows iteration
#' @param z0 Initialization of configuration; NULL by default
#'
#' @return 2D representation vector, size of N by 2
#' @importFrom stats dist
#' @importFrom stats cmdscale
#' @importFrom phyloseq distance sample_data
#' @export
#' @examples
#' set.seed(100)
#' require(phyloseq)
#' data(microbiome)
#' D <- distance(microbiome, method = 'wunifrac') # requires phyloseq package
#' y <- sample_data(microbiome)$Treatment
#' z0 <- cmdscale(d = D)
#' fmds(z0 = z0, D = D, y = y)
fmds <- function(D = NULL, y = NULL, X, nit = 100, lambda = 0.5, threshold_p = 0.05, z0 = NULL){
    if(is.null(y)){
        stop("The label input is not provided.")
    }
    if(is.null(D)){
        D <- getDistMat(X)
    } else {
        D <- as.matrix(D)
    }
    if(is.null(z0)){
        z0 <- cmdscale(d = D)
    }
    N <- dim(z0)[1]
    S <- dim(z0)[2]
    a <- length(unique(y))
    if(a == 1){
        stop("The label input is inappropriate: expected two or more groups.")
    }
    y_indmat <- getIndMat(y)
    f_ratio <- pseudoF(z = X, D = D, y = y)
    z_temp <- z_up <- z0
    p0 <- getP(D = D, y = y)$p
    if (p0 > 0.1) {
        message("p-value >= .1: FMDS is not necessary, returning z0 or MDS output.")
        return(z0)     # <-- This exits the function immediately
    }
    log_iter_mat <- matrix(0, nrow=0, ncol=6)
    colnames(log_iter_mat) <-
        c('epoch', 'obj', 'obj_mds', 'obj_confr', 'p_z', 'p_0')
    p_prev <- 1
    for(t in seq(0, nit)){
        p_up <- getP(z = z_up, y = y)$p

        if((abs(p_up-p0) >= abs(p_prev-p0)) & (abs(p_prev-p0) <= threshold_p)){
            message(sprintf('Lambda %.2f ...halt iteration', lambda))
            z_up <- z_prev # revert to prev
            break
        }

        list_pair <- pairByRank(D = D, z = z_up, y = y) # _0, _z
        ind_f_ratio <- which.min(abs(f_ratio - list_pair[,1]))[1]
        f_ratio_pred <- list_pair[,2][ind_f_ratio]

        z_distmat <- as.matrix(dist(z_up))
        f_diff_nominator <- sum((1 - a * y_indmat * (1+f_ratio_pred*(a-1)/(N-a))) *
                                    z_distmat^2)
        delta <- sign(f_diff_nominator)
        obj_conf <- abs(f_diff_nominator)
        obj_mds <- mdsObj(D, z_up, N)
        obj <- lambda*obj_conf + obj_mds

        message(paste('epoch', t,
                    '  lambda', lambda,
                    '  total', sprintf(obj, fmt = '%#.2f'),
                    '  mds', sprintf(obj_mds, fmt = '%#.2f'),
                    '  conf', sprintf(obj_conf, fmt = '%#.2f'),
                    '  p_z', sprintf(p_up, fmt = '%#.3f'),
                    '  p_0', sprintf(p0, fmt = '%#.3f')
        ))
        log_iter_mat <- rbind(log_iter_mat,
                              c(t, obj, obj_mds, obj_conf, p_up, p0))


        for(i in seq(1, N)){
            z_distmat <- as.matrix(dist(z_up))  # (N,N)
            coeff <- D/z_distmat  # final term in the update
            coeff[is.nan(coeff)] <- 0
            z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")

            z_temp[i,] <- (1+lambda*delta) * (apply(z_up[y!=y[i],], 2, sum)) +
                (1 - lambda * delta * (a-1 + f_ratio_pred*(a-1)/(N-2))) *
                (apply(z_up[y==y[i],], 2, sum)) +
                apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)

            z_temp[i,] <- z_temp[i,] / (N - N*lambda*delta*f_ratio_pred/(N-a))
        } # end z_temp

        z_prev <- z_up
        p_prev <- p_up
        z_up <- z_temp
    } # end iteration

    return(z_up)
}
