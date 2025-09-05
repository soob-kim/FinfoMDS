test_that("fmds returns correct structure", {
    y = c(rep(0, 10), rep(1, 10))
    X = rbind(matrix(rnorm(30), ncol = 3),
              matrix(rnorm(30, mean = 1), ncol = 3))
    res <- fmds(y = y, X = X)

    expect_true(is.matrix(res))
    expect_equal(dim(res), c(20, 2))
    expect_false(anyNA(res))
})

test_that("fmds is reproducible with fixed seed", {
    y = c(rep(0, 10), rep(1, 10))
    X = rbind(matrix(rnorm(30), ncol = 3),
              matrix(rnorm(30, mean = 1), ncol = 3))
    set.seed(123)
    res1 <- fmds(y = y, X = X)
    set.seed(123)
    res2 <- fmds(y = y, X = X)

    expect_equal(res1, res2)
})

test_that("fmds preserves distances at least somewhat", {
    y = c(rep(0, 10), rep(1, 10))
    X = rbind(matrix(rnorm(30), ncol = 3),
              matrix(rnorm(30, mean = 1), ncol = 3))
    res <- fmds(y = y, X = X)

    d_orig <- dist(X)
    d_emb  <- dist(res)

    cor_val <- cor(as.vector(d_orig), as.vector(d_emb))
    expect_gt(cor_val, 0.5)   # loose threshold to catch total failures
})

test_that("fmds rejects single-label input", {
    y = rep(0, 10)
    X = matrix(rnorm(30), ncol = 3)
    expect_error(
        fmds(y = y, X = X),
        regexp = "inappropriate"
    )
})

test_that("fmds handles degenerate inputs", {
    y = c(rep(0, 10), rep(1, 10))
    X = rbind(matrix(0, ncol = 3, nrow = 10),
              matrix(1, ncol = 3, nrow = 10))
    res <- fmds(y = y, X = X)
    expect_equal(nrow(res), 20)
})
