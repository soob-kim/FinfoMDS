test_that("PERMANOVA p-value works", {
    z <- c(rnorm(10), rnorm(10, mean = 3))
    y <- c(rep(0, 10), rep(1, 10))
    p_val <- getP(z = z, y = y)$p
    if(p_val <= 0.05){
        flag = 1
    } else {
        flag = 0
    }
    expect_equal(flag, 1)
})
