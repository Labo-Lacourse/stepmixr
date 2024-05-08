### Stepmixr.
library(stepmixr)

### To be corrected in Stepmixr. 
data_bakk_response <- function (n_samples, sep_level, n_classes = 3, n_mm = 6,
                                random_state = NULL) {
  sm <- try(reticulate::import("stepmix"), silent = TRUE)
  if (!is.null(random_state))
    random_state = as.integer(random_state)
  if (inherits(sm, "try-error"))
    stop(paste("Unable to find stepmix library in your python repos\n",
               "Install it using pip install stepmix", collapse = ""))
  sm$datasets$data_bakk_response(as.integer(n_samples), sep_level,
                                 as.integer(n_classes), as.integer(n_mm), random_state)
}

### Simulation taking different parameters into account. 
simul  <- function(n = 2000, n_steps = 1, sep_level = 0.9, 
                   correction = NA, rs = 1){
  # Simulate data
  datasim <- data_bakk_response(n_samples=as.integer(n), 
                                sep_level=sep_level, 
                                random_state = as.integer(rs))

  # Fit StepMix Estimator
  if(!is.na(correction))
    model = stepmix(n_components=3,
                    measurement='binary',
                    structural='gaussian_unit',
                    n_steps=n_steps,
                    max_iter=250,  # Latent Gold default : 250 EM iterations + 50 NR iterations
                    abs_tol=1e-8,
                    progress_bar=0,
                    correction = correction,
                    random_state = as.integer(rs))
  else
    model = stepmix(n_components=3,
                    measurement='binary',
                    structural='gaussian_unit',
                    max_iter=250,
                    abs_tol=1e-8,
                    n_steps=n_steps,
                    progress_bar=0,
                    correction = NULL,
                    random_state = as.integer(rs))
  fit1 <- fit(model, datasim[[1]], datasim[[2]])

  # Retrieve mean parameters
  mus = fit1$get_parameters()[["structural"]][["means"]]

  # Return the maximum value.
  return(max(mus))
}

### Paramater grid.
test.par <- expand.grid(n         = c(500, 1000, 2000), # Number of participants
                        sep_level = c(0.7, .8, .9),     # Separation level
                        n_steps   = c(1, 2, 3),         # Stepwise approach
                        correction = c(NA, "BCH", "ML"))# Correction
test.par <- test.par[test.par$n_steps %in%  1:2 & is.na(test.par$correction) |
           test.par$n_steps == 3,]

# ### Run the 500 simulations.
test = NULL
for( i in 1:500){
  print(i)
  rs = 12345 * (i - 1)
  test <- cbind(test, mapply(simul,
                             n          = test.par$n,
                             sep_level  = test.par$sep_level,
                             n_steps    = test.par$n_steps,
                             correction = test.par$correction,
                             rs = rs))
}

test[test>=1000] <- NA #max= 1.806421, min -0.004786873

### Compute Mean Bias and RMSE
test.par$resMean_bias = apply(1-test, 1, mean, na.rm = TRUE)
test.par$rmse = sqrt(apply((1-test)^2, 1, mean, na.rm = TRUE))

### Combine n_steps, correction into method
test.par$method = factor(with(test.par, sprintf("%s (%s)", n_steps, correction)),
                         levels = c("1 (NA)",
                                    "2 (NA)",
                                    "3 (NA)",
                                    "3 (BCH)",
                                    "3 (ML)"))

save(test, test.par, file = "./Distal_outcome.Rdata")

### Format results.
(ft1 <- round(ftable(xtabs(-resMean_bias ~ sep_level + n + method, test.par)), 2))
(ft2 <- round(ftable(xtabs(rmse ~ sep_level + n + method, test.par)), 2))

