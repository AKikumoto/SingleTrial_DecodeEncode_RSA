library(here)
source(here("utils.R"))
source(here("params.R"))

## simulate ----


params <- new_params(
  l = get_conditions_axes("congruency_3by3"),
  feature_weights = c(
    "(Intercept)" = 0, "f1a-c" = 1, "f1b-c" = 1, "f2a-c" = 1, "f2b-c" = 1,
    "f1a-c:f2a-c" = 0, "f1b-c:f2a-c" = 0, "f1a-c:f2b-c" = 0, "f1b-c:f2b-c" = 0, "congruency" = 0
    ),
  sigma_type = "uniform", sigma_args = 1,
  n_trial = 10, n_channel = 20, n_sim = 100,
  rsa_model_names = c("f1", "f2", "congruency", "id"),
  encoding_contrast_names = list(
    f1 = c("f1a-c", "f1b-c"),
    f2 = c("f2a-c", "f2b-c"),
    id = c("f1a-c:f2a-c", "f1b-c:f2a-c", "f1a-c:f2b-c", "f1b-c:f2b-c")
  )
)

params  ## holds all parameters for simulation


## initial encoding models (for both decoding+RSA and encoding)

Y <- list(train = generate_data(params), test = generate_data(params))
fit <- lapply(Y, \(y, x = params$X) fit_encoder(X = x, Y = y))
trial_condition <- dummy2factor(params$X)

## decoding + rsa

rsa_coefs <- fit_rsa(Y$test, fit$train, params, trial_condition)

## encoding

encoding_projs <- fit_encoding(Y$test, fit$train, params, trial_condition)

## summarize

rsa_sum <- rsa_coefs %>% 
  summarize(across(where(is.numeric), \(x) mean(x) / sd(x))) %>%
  mutate(method = "rsa")
encoding_sum <- encoding_projs %>% 
  summarize(across(where(is.numeric), \(x) mean(x) / sd(x))) %>%
  mutate(method = "encoding")
sim_sum <- bind_rows(rsa_sum, encoding_sum)
