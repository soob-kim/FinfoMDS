#' Get index matrix
#'
#' @param y a vector of treatments of length N
#'
#' @return A N by N matrix of indicators of equal treatment
#' @export
#' @examples
#' getIndMat(microbiome$host)
getIndMat <- function(y){
    N <- length(y)
    mat <- matrix(0, nrow = N, ncol = N)
    for(i in seq(1, N, by = 1)){
        for(j in seq(1, N, by = 1)){
            mat[i,j] <- as.numeric(y[i]==y[j])
        }
    }
    return(mat)
}


#' Get pseudo-F value for PERMANOVA
#'
#' @param z Object matrix; used to build distance matrix d; d is prioritized
#' @param D Distance matrix; if NULL, obtain from mat using Euclidean distance
#' @param y Vector of treatments
#'
#' @return pseudo-F value
#' @export
#' @examples
#' pseudoF(D = microbiome$dist, y = microbiome$host)
pseudoF <- function(z=NULL, D = NULL, y){
    if(is.null(D)){
        D <- as.matrix(dist(z))
    } else {
        D <- as.matrix(D)
    }
    N <- nrow(D)
    a <- length(unique(y))
    n <- N/a

    SST <- SSW <- 0
    SST <- sum(D * D) / 2
    SSW <- sum(D * D * getIndMat(y)) / 2
    SST <- SST / N
    SSW <- SSW / n
    SSA <- SST - SSW
    return((SSA * (N-a))/(SSW * (a-1)))
}


#' Get p-value of PERMANOVA
#'
#' @param z Object matrix; used to build distance matrix d; d is prioritized
#' @param D Distance matrix; if NULL, obtain from mat using Euclidean distance
#' @param y Vector of treatments
#' @param n_iter Number of iterations; defaults to 999
#'
#' @return list of ratio_all: vector of obtained pseudo-F values from permutations,
#' ratio: pseudo-F value, p: p-value from PERMANOVA
#' @export
#' @examples
#' getP(D = microbiome$dist, y = microbiome$host)
getP <- function(z = NULL, D = NULL, y, n_iter = 999){
    # initialize
    f_permuted <- matrix(0, nrow = n_iter, ncol = 1)  # pseudo-F only
    # iterate to get pseudo F
    N <- length(y)
    a <- length(unique(y))
    tbl <- table(y)
    for (iter in 1:n_iter){
        y_rand <- rep(1, N)
        pool_v <- seq(1, N)
        for(cl in seq(2, a)){
            ind_rand <- sample(pool_v, tbl[cl], replace=FALSE)
            y_rand[ind_rand] <- cl
            pool_v <- setdiff(pool_v, ind_rand)
        }
        f_permuted[iter, 1] <- pseudoF(z = z, D = D, y = y_rand)
    }
    # compute p value
    f_sorted <- sort(f_permuted[, 1], decreasing = TRUE)
    f_val <- pseudoF(z = z, D = D, y = y)
    p_val <- which(f_val > f_sorted)[1]
    p_val <- (p_val-1)/(n_iter+1)

    return(list(ratio_all = f_sorted, ratio = f_val, p = p_val))
}
