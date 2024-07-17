new_params <- function(
  l,
  n_trial, n_channel, n_sim,
  feature_weights = NULL,
  sigma_type = "spherical",
  sigma_args = 1,
  scale_args = NULL,
  shift_args = NULL,
  rotation_args = NULL,
  rsa_model_names = NULL,
  encoding_contrast_names = NULL
) {

  conditions <- l$conditions
  A <- l$A

  ## numbers
  n <- c(trial = n_trial, channel = n_channel, sim = n_sim, condition = nrow(conditions))

  ## design matrix: maps each condition to a trial
  X <- model.matrix(~ 0 + factor(rep(1:n["condition"], each = n["trial"])))
  colnames(X) <- conditions$id

  ## scale axes non-selectively (i.e., in all conditions)
  if (!is.null(feature_weights)) {
    if (!is.null(names(feature_weights))) stopifnot(all(names(feature_weights) == colnames(A)))
    feature_weights <- matrix(rep(feature_weights, nrow(A)), nrow = nrow(A), byrow = TRUE)
  } else {
     feature_weights <- 1
  }

  ## selectively scale axes in specific conditions (via hadamard product)
  if (is.null(scale_args)) {
    B_scale <- 1
  } else {
    B_scale <- do.call(matrix_mask, c(list(conditions = conditions, axes = colnames(A), rest_ones = TRUE), scale_args))
  }
  
  ## selectively shift axes in specific conditions (via element-wise addition)
  if (is.null(shift_args)) {
    B_shift <- 0
  } else {
    B_shift <- do.call(matrix_mask, c(list(conditions = conditions, axes = colnames(A), rest_ones = FALSE), shift_args))
  }

  ## rotate one axis in the direction of another (via interpolation)
  if (is.null(rotation_args)) {
    R <- diag(ncol(A))
  } else {
    R <- do.call(interpolation, rotation_args)
  }
  dimnames(R) <- list(colnames(A), colnames(A))

  ## noise covariance matrix
  if (sigma_type == "uniform") {
    Sigma <- diag(n["channel"]) * sigma_args
  }
  sigma <- pracma::sqrtm(Sigma)$B

  ## compute means mu
  D <- basis(n["channel"], n["condition"] + 1)  ## channel basis matrix
  A_trans <- (A * feature_weights * B_scale + B_shift) %*% R
  mu <- tcrossprod(A_trans, D) # condition by channel means

  if (!is.null(rsa_model_names)) {
    rsa_models <- create_rsa_models(conditions, rsa_model_names)
  }
  
  if (!is.null(encoding_contrast_names)) {
    encoding_contrasts <- create_encoding_contrasts(A[, 1:nrow(A)], encoding_contrast_names)
  }
  
  ## return
  list(
    ## design:
    conditions = conditions, n = n, A = A, X = X, 
    ## data:
    A_trans = A_trans, mu = mu, sigma = sigma,
    ## rsa, encoding:
    rsa_models = rsa_models, encoding_contrasts = encoding_contrasts
    )

}

