library(here)
library(matrixStats)
library(pracma)
library(data.table)
library(future)
library(doFuture)
library(dplyr)
library(ggplot2)

source(here("utils.R"))
source(here("params.R"))

plan(multicore(workers = length(future::availableWorkers())))

## simulate ----


params <- new_params(
  l = get_conditions_axes("congruency_3by3"),
  feature_weights = c(
    "(Intercept)" = 0, "f1a-c" = 1, "f1b-c" = 1, "f2a-c" = 1, "f2b-c" = 1,
    "f1a-c:f2a-c" = 1, "f1b-c:f2a-c" = 1, "f1a-c:f2b-c" = 1, "f1b-c:f2b-c" = 1,
    "congruency" = 0
    ),
  sigma_type = "uniform", sigma_args = 1000,
  n_trial = 50, n_channel = 100, n_sim = 1E3,
  rsa_model_names = c("f1", "f2", "congruency", "id"),
  encoding_contrast_names = list(
    f1 = c("f1a-c", "f1b-c"),
    f2 = c("f2a-c", "f2b-c"),
    id = c("f1a-c:f2a-c", "f1b-c:f2a-c", "f1a-c:f2b-c", "f1b-c:f2b-c")
  )
)

str(params)  ## holds all parameters for simulation

res <- simulate_singletrial(params)

res %>%
  filter(term != "intercept") %>%
  ggplot(aes(x = term, effsize, color = method)) +
  geom_boxplot(width = 0.25)
