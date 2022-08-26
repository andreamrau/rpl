
## 1) Simulate exchangeable correlation matrix ---------------------------------
d <- 10
sim_exchangeable <-
    simulate_mvt_poisson(n = 1000, d = d,
                         param = list(margins = rep(1, d), disp = 0.5),
                         type = "exchangeable")

init_exchangeable <- init(sim_exchangeable$data, type = "exchangeable")
run_exchangeable <- rpl_optim(init_exchangeable, sim_exchangeable$data, pi = 0.1)

## 2) Simulate unstructured correlation matrix ---------------------------------
d <- 4
sim_unstructured <-
    simulate_mvt_poisson(n = 1000, d = d,
                         param = list(margins = rep(1, d),
                                      disp = c(0.4,0.5,0.6,0.7,0.8,0.9)),
                         type = "unstructured")

init_unstructured <- init(sim_unstructured$data, type = "unstructured")
run_unstructured <- rpl_optim(init_unstructured, sim_unstructured$data, pi = 0.1)

## Calculate SE of estimated parameters
se_unstructured <- rpl_se(mydata=sim_unstructured$data,
                          parameters=run_unstructured$par,
                          pi = 0.1,
                          offsets=NULL, verbose=FALSE)


## 3) Simulate one-factor correlation matrix -----------------------------------
d <- 4
sim_factor <-
    simulate_mvt_poisson(n = 1000, d = d,
                         param = list(margins = rep(1, d),
                                      disp = seq(0, 1, length = d)),
                         type = "factor")

init_factor <- init(sim_factor$data, type = "factor")
run_factor <- rpl_optim(init_factor, sim_factor$data, pi = 0.1)

## Calculate SE of estimated parameters
se_factor <- rpl_se(mydata=sim_factor$data,
                          parameters=run_factor$par,
                          pi = 0.1,
                          offsets=NULL, verbose=FALSE)

## 4) Simulate block-exchangeable correlation matrix ---------------------------
d <- 10
block_indices <- list(c(1:5), c(6:10))
sim_block <-
    simulate_mvt_poisson(n = 1000, d = d,
                         block_indices = block_indices,
                         param = list(margins = rep(1, d),
                                      disp = c(0.5,0.1,0.8)),
                         type = "block_exchangeable")

init_block<- init(sim_block$data, type = "block_exchangeable",
                  block_indices = block_indices)
run_block <- rpl_optim(init_block, sim_block$data, pi = 0.1,
                       block_indices = block_indices)
