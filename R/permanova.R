#' Get index matrix
#'
#' @param y a vector of treatments of length N
#'
#' @return A N by N matrix of indicators of equal treatment
#' @export
get_ind_mat <- function(y){
    N <- length(y)
    mat <- matrix(0, nrow = N, ncol = N)
    for(i in 1:N){
        for(j in 1:(i-1)){
            mat[i,j] <- as.numeric(y[i]==y[j])
        }
    }
    return(mat + t(mat))
}


#' Title Get pseudo-F value for PERMANOVA
#'
#' @param mat Object matrix; used to build distance matrix d; d is prioritized
#' @param trt Vector of treatments
#' @param d Distance matrix; if NULL, obtain from mat using Euclidean distance
#'
#' @return pseudo-F value
#' @export
pseudo_F <- function(mat=NULL, trt, d = NULL){
    if(is.null(d)){
        d <- as.matrix(dist(mat))
    } else {
        d <- as.matrix(d)
    }
    N <- nrow(d)
    a <- length(unique(trt))
    n <- N/a

    SST <- SSW <- 0
    SST <- sum(d * d) / 2
    SSW <- sum(d * d * get_ind_mat(trt)) / 2
    SST <- SST / N
    SSW <- SSW / n
    SSA <- SST - SSW
    return((SSA * (N-a))/(SSW * (a-1)))
}


#' Get p-value of PERMANOVA
#'
#' @param mat Object matrix; used to build distance matrix d; d is prioritized
#' @param d Distance matrix; if NULL, obtain from mat using Euclidean distance
#' @param trt Vector of treatments
#' @param n_iter Number of iterations; defaults to 999
#'
#' @return list of ratio_all: vector of obtained pseudo-F values from permutations,
#' ratio: pseudo-F value, p: p-value from PERMANOVA
#' @export
get_p <- function(mat=NULL, d=NULL, trt, n_iter=999){
    # initialize
    trt <- as.matrix(trt)
    f_permuted <- matrix(0, nrow=n_iter, ncol=1)  # pseudo-F only
    # iterate to get pseudo F
    N <- length(trt)
    a <- length(unique(trt))
    tbl <- table(trt)
    for (iter in 1:n_iter){
        y_rand <- rep(1,N)
        pool_v <- 1:N
        for(cl in 2:a){
            ind_rand <- sample(pool_v, tbl[cl], replace=F)
            y_rand[ind_rand] <- cl
            pool_v <- setdiff(pool_v, ind_rand)
        }
        f_permuted[iter,1] <- pseudo_F(mat=mat, d = d, trt = y_rand)
    }
    # compute p value
    f_sorted <- sort(f_permuted[,1], decreasing = TRUE)
    f_val <- pseudo_F(mat=mat, d = d, trt = trt)
    p_val <- which(f_val > f_sorted)[1]
    p_val <- (p_val-1)/(n_iter+1)

    return(list(ratio_all = f_sorted, ratio = f_val, p = p_val))
}
