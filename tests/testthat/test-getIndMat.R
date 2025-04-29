test_that("Obtaining index matrix from treatment vector", {
    mat <- getIndMat(c(1, 1, 0, 0))
    expect_equal(mat, cbind(c(1,1,0,0), c(1,1,0,0), c(0,0,1,1), c(0,0,1,1)))
})
