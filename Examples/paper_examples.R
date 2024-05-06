library(stepmixr)

#Response example
datasim = data_bakk_response(
  n_samples = 2000, 
  sep_level = 0.9, 
  random_state = 42
  )

Y = datasim[[1]]
Z_o = datasim[[2]]

model = stepmix(
  n_components = 3, 
  measurement = 'binary', 
  structural = 'gaussian_unit', 
  n_steps = 1, 
  random_state = 42, 
  verbose = 1
  ) 

fit1 = fit(model, Y, Z_o)

mus = fit1$get_sm_df()
mus

#Covariate example
datasim = data_bakk_covariate(
  n_samples = 2000, 
  sep_level = 0.9, 
  random_state = 42
)

Y = datasim[[1]]
Z_p = datasim[[2]]

covariate_params = list(
  method = 'newton-raphson', 
  max_iter = as.integer(1), 
  intercept = TRUE
)

model = stepmix(
  n_components = 3, 
  measurement = 'binary', 
  structural = 'covariate', 
  n_steps = 1, 
  random_state = 42, 
  structural_params = covariate_params
)

fit1 = fit(model, Y, Z_p)

betas = fit1$get_parameters()

betas = betas$structural$beta

betas = betas - c(1, 1, 1) %o% betas[2, ]
betas

#Complete example
datasim = data_bakk_complete(
  n_samples = 2000, 
  sep_level = 0.8, 
  nan_ratio = 0.25, 
  random_state = 42
  )

Y = datasim[[1]]
Z = datasim[[2]]

structural_descriptor = list(
  covariate = list( 
    model = "covariate", 
    n_columns = as.integer(1), 
    method = "newton-raphson", 
    max_iter = as.integer(1)
  ),
  response = list(
    model = "gaussian_unit_nan", 
    n_columns = as.integer(1)
  )
)

model = stepmix(
  n_components = 3, 
  measurement = "binary_nan", 
  structural = structural_descriptor, 
  n_steps = 1, 
  random_state = 42
  )

fit1 = fit(model, Y, Z)

mus = fit1$get_parameters()
mus = mus$structural$response
mus