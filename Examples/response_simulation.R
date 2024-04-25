library(stepmixr)


n = 2000

### À corriger dans stepmix. À utiliser à la place de celle de Stepmix.
data_bakk_response <- function (n_samples, sep_level, n_classes = 3, n_mm = 6, 
                                random_state = NULL) 
{
  sm <- try(reticulate::import("stepmix"), silent = TRUE)
  if (!is.null(random_state)) 
    random_state = as.integer(random_state)
  if (inherits(sm, "try-error")) 
    stop(paste("Unable to find stepmix library in your python repos\n", 
               "Install it using pip install stepmix", collapse = ""))
  sm$datasets$data_bakk_response(as.integer(n_samples), sep_level, 
                                 as.integer(n_classes), as.integer(n_mm), random_state)
}


### Simulation qui tient compte des différents paramètre. 
### R NA   -> Python NaN
### R NULL -> Python None

simul  <- function(n = 2000, n_steps = 1, sep_level = 0.9, correction = NA){
  # Simulate data
  datasim <- data_bakk_response(n_samples=as.integer(n), sep_level=sep_level)
  
  # Fit StepMix Estimator
  if(!is.na(correction))
    model = stepmix(n_components=3, measurement='binary', structural='gaussian_unit', 
                    n_steps=n_steps, correction = correction)
  else 
    model = stepmix(n_components=3, measurement='binary', structural='gaussian_unit', 
                    n_steps=n_steps, correction = NULL)
  fit1 <- fit(model, datasim[[1]], datasim[[2]])
  
  # Retrieve mean parameters
  mus = fit1$get_parameters()[["structural"]][["means"]]  
  
  # Return the maximum value.
  return(max(mus))
}

### Paramater grid. 
test.par <- expand.grid(n         = c(500, 1000, 2000), # Number of participants
                        sep_level = c(0.7, .8, .9),     # Separation level
                        n_steps   = c(1, 2, 3, 3, 3))   # Stepwise approach
test.par$correction <- c(rep(NA, 18), rep(c(NA, "BCH", "ML"), 9)) # Correction.


### Run the 500 simulations.
test <- replicate(2, expr = mapply(simul,
                                   n          = test.par$n,
                                   sep_level  = test.par$sep_level,
                                   n_steps    = test.par$n_steps,
                                   correction = test.par$correction))


test[test>=1000] <- NA

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



#### Exemple complémentaire de R.

### Exemple de loop.
for(i in 1:10){
  for(j in 1:5)
    print(sprintf("i:%d j:%d", i, j))
}

for(let in letters)
  print(let)


### python
### for let in letters:
###     print(let)


replicate(10000, sum(sample(1:1000)))

for(i in 1:10000)
  sum(sample(1:1000))


### Exemple avec boucle for.
for(i in 1:45)
  simul(n = test.par$n[i], 
        sep_level = test.par$n[i])


simul2  <- function(n = 2000, n_steps = 1, sep_level = 0.9, correction = NA){
  # Simulate data
  datasim <- data_bakk_response(n_samples=as.integer(n), sep_level=sep_level)
  
  # Fit StepMix Estimator
  if(!is.na(correction))
    model = stepmix(n_components=3, measurement='binary', structural='gaussian_unit', 
                    n_steps=n_steps, correction = correction, verbose = 0, progress_bar = FALSE)
  else 
    model = stepmix(n_components=3, measurement='binary', structural='gaussian_unit', 
                    n_steps=n_steps, correction = NULL, verbose = 0, progress_bar = FALSE)
  X <- datasim[[1]]
  Y <- datasim[[2]]
  fit1 <- fit(model, X, Y)
  ll = fit1$score(X, Y) * n
  ll
}

ll10 <- replicate(10, simul2())
hist(ll10)

