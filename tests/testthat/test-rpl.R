##-------------------------------------------------------------------------------
context("rpl simulation tests")

# set.seed(12345)
# sim_6a <- mv_simulate(type = "D6")
# sim_6b <- mv_simulate(type = "D6", beta=7, n=200, K=5, sigma=0.5)
#
# test_that("delta and sigma are nonnnegative", {
#     expect_error(mv_simulate(type = "D6", beta=-1))
#     expect_error(mv_simulate(type = "D6", sigma=-1))
# })
#
# test_that("n and K are positive scalars", {
#     expect_error(mv_simulate(type = "D6", K=-1))
#     expect_error(mv_simulate(type = "D6", K=0))
#     expect_error(mv_simulate(type = "D6", n=-1))
#     expect_error(mv_simulate(type = "D6", n=0))
# })