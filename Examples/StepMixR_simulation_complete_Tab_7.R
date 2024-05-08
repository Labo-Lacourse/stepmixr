library(stepmixr)

simul  <- function(n = 2000, n_steps = 1, sep_level = 0.8, nan_ratio=0.00, correction = NA, rs = 1){
  # Simulate data
  datasim <- data_bakk_complete(n_samples=as.integer(n), 
                                sep_level = 0.8, 
                                nan_ratio=nan_ratio, 
                                random_state = as.integer(rs))
  
  cov_params = list(model = 'covariate',
                    method = 'newton-raphson',
                    n_columns = as.integer(1),
                    max_iter  = as.integer(1))
  
  dist_param = list(model = 'gaussian_unit_nan',
                    n_columns = as.integer(1))
  
  structural_descriptor = list(covariate = cov_params,
                               response = dist_param)
  # Fit StepMix Estimator
  if(!is.na(correction))
    model = stepmix(n_components=3, 
                    measurement='bernoulli_nan', 
                    structural= structural_descriptor, 
                    max_iter=500,  
                    abs_tol=1e-8,
                    n_steps=3, 
                    progress_bar=0,
                    correction = correction, 
                    random_state = as.integer(rs)) 
  else 
    model = stepmix(n_components=3, 
                    measurement='bernoulli_nan', 
                    structural= structural_descriptor, 
                    max_iter=500,  
                    abs_tol=1e-8,
                    n_steps=n_steps, 
                    progress_bar=0,
                    correction = NULL,
                    random_state = as.integer(rs))
  fit1 <- fit(model, datasim[[1]], datasim[[2]])
  
  means = fit1$get_parameters()[["structural"]][["response"]][["means"]]
  dist = (max(means))
  
  # Return the maximum value.
  return(dist)
  
}

### Paramater grid. 
test.par <- expand.grid(n         = c(500, 1000, 2000), # Number of participants
                        nan_ratio = c(0.00, 0.25, 0.50),
                        n_steps   = c(1, 2, 3),          # Stepwise approach
                        correction = c(NA, "BCH", "ML"))  # Correction. 

test.par <- test.par[test.par$n_steps %in%  1:2 & is.na(test.par$correction) |
                       test.par$n_steps == 3,]



### Run the 500 simulations.
test = NULL
for(i in 1:500){
  print(i)
  rs <- 12345 * (i - 1)
  test <- cbind(test, mapply(simul,
                             n          = test.par$n,
                             nan_ratio  = test.par$nan_ratio,
                             n_steps    = test.par$n_steps,
                             correction = test.par$correction,
                             rs = rs))
  
}

test[test>=1000] <- NA 


### Compute Mean Bias and RMSE
test.par$Dist_resMean_bias = apply(1-test, 1, mean, na.rm = TRUE)
test.par$Dist_rmse = sqrt(apply((1-test)^2, 1, mean, na.rm = TRUE))

### Combine n_steps, correction into method
test.par$Dist_method = with(test.par, sprintf("%s (%s)", n_steps, correction))

### Format results
test.par$Dist_method <- factor(test.par$Dist_method,
                               levels = c("1 (NA)",
                                          "2 (NA)",
                                          "3 (NA)",
                                          "3 (BCH)",
                                          "3 (ML)"))

### Table 7
(ft1 <- round(ftable(xtabs(Dist_resMean_bias ~ nan_ratio + n + Dist_method, test.par)), 2))
(ft2 <- round(ftable(xtabs(Dist_rmse ~ nan_ratio + n + Dist_method, test.par)), 2))

### Pour tableau en LaTeX. 
xtable::xtableFtable(ft1)
xtable::xtableFtable(ft2)

save(test, test.par, file = "complete.Rdata")
