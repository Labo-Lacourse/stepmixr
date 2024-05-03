library(stepmixr)

simul  <- function(n = 2000, n_steps = 1, sep_level = 0.9, correction = NA){
  # Simulate data
  datasim <- data_bakk_covariate(n_samples=as.integer(n), sep_level=sep_level)
  opt_params = list(
    method = 'newton-raphson',
    max_iter  = as.integer(1),
    intercept = TRUE)
  
  # Fit StepMix Estimator
  if(!is.na(correction))
    
    model = stepmix(n_components=3, 
                       measurement='binary', 
                       structural='covariate', 
                       max_iter=250,  
                       abs_tol=1e-8,
                       n_steps=n_steps, 
                       progress_bar=0,
                       structural_params = opt_params,
                      correction = correction)
  else 
    model = stepmix(n_components=3, 
                  measurement='binary', 
                  structural='covariate', 
                  max_iter=250,  
                  abs_tol=1e-8,
                  n_steps=n_steps, 
                  progress_bar=0,
                  structural_params = NULL,
                  correction = NULL)
  
  fit1 <- fit(model, datasim[[1]], datasim[[2]])
  
  # Retrieve mean parameters
  coeff = fit1$get_parameters()[["structural"]][["beta"]]
  coeff2 = identify_coef(coeff)
  mus = max(coeff2[, 2])
  
  # Return the maximum value.
  return(max(mus))
}


### Paramater grid. 
test.par <- expand.grid(n          = c(500, 1000, 2000), # Number of participants
                        sep_level  = c(0.7, .8, .9),     # Separation level
                        n_steps    = c(1, 2, 3),   # Stepwise approach
                        correction = c(NA, "BCH", "ML"))   

test.par <- test.par[test.par$n_steps %in%  1:2 & is.na(test.par$correction) |
                       test.par$n_steps == 3,]

### Run the 500 simulations.
test <- replicate(500, expr = mapply(simul,
                                     n          = test.par$n,
                                     sep_level  = test.par$sep_level,
                                     n_steps    = test.par$n_steps,
                                     correction = test.par$correction))

test[test>=1000] <- NA 

save(test, test.par, file = "Distal outcome/")

### Compute Mean Bias and RMSE

test.par$resMean_bias = apply(1-test, 1, mean, na.rm = TRUE)
test.par$rmse = sqrt(apply((1-test)^2, 1, mean, na.rm = TRUE))

### Combine n_steps, correction into method
test.par$method = with(test.par, sprintf("%s (%s)", n_steps, correction))

### Format results. 
(ft1 <- round(ftable(xtabs(resMean_bias ~ sep_level + n + method, test.par)), 2))
(ft2 <- round(ftable(xtabs(rmse ~ sep_level + n + method, test.par)), 2))


### Pour tableau en LaTeX. 
xtable::xtableFtable(ft1)

write.csv(test, "StepMixR simulation results data/test_df_cov.csv")
write.csv(test.par, "StepMixR simulation results data/testpar_df_cov.csv")
