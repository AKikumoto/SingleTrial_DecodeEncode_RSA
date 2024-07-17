matrix_mask <- function(
  conditions, condition_name, condition_val, axes, axis_name, mask_val,
  rest_ones = FALSE
  ) {
  is_condition <- conditions[[condition_name]] == condition_val
  is_axis <- axes == axis_name
  mask_01 <- is_condition %*% t(is_axis)
  mask <- mask_01 * mask_val
  if (rest_ones) mask <- mask + (1 - mask_01)
  mask
}


interpolation <- function(size, name_rotated, name_target, prop_alignment, nms = NULL) {
  stopifnot(prop_alignment > 0 && prop_alignment <= 1)
  R <- diag(size)
  if (!is.null(nms)) dimnames(R) <- list(nms, nms)
  R[name_rotated, name_rotated] <- 1 - prop_alignment
  R[name_target, name_rotated] <- prop_alignment
}


basis <- function(n_row, n_col, centered = FALSE) {
  x <- matrix(rnorm(n_row * n_col), nrow = n_row)
  if (centered) x <- sweep(x, 2, colMeans(x), "-")
  qr.Q(qr(x))
}



get_conditions_axes <- function(design_type) {

  # conditions grid: holds constellations of unique conditions
  if (design_type == "congruency_3by3") {

    conditions <- expand.grid(f1 = c("a", "b", "c"), f2 = c("a", "b", "c"))
    conditions$id <- paste0(conditions$f1, "_", conditions$f2)
    conditions$congruency <- ifelse(conditions$f1 != conditions$f2, "incon", "congr")

    ## feature matrix: maps each condition to coordinates in representational space
    contrast_matrix <- contr.sum(3)
    dimnames(contrast_matrix) <- list(c("a", "b", "c"), c("a-c", "b-c"))
    contrasts(conditions$f1) <- contrast_matrix
    contrasts(conditions$f2) <- contrast_matrix
    congruency <- sign((conditions$congruency == "incon") - 0.5)  ## contrast code
    A <- cbind(model.matrix(~ f1 * f2, conditions), congruency = congruency)
    rownames(A) <- conditions$id

  }

  list(A = A, conditions = conditions)

}



generate_data <- function(params) {
  eps <- basis(params$n["trial"] * params$n["condition"], params$n["channel"], centered = TRUE) %*% params$sigma
  Y <- params$X %*% params$mu + eps
  Y
}


fit_encoder <- function(X, Y) {
  fit <- .lm.fit(X, Y)
  rownames(fit$coefficients) <- colnames(X)
  fit
}


pdist <- function(A, B) {
  
  ## squared euclidean distances between each row of A and each row of B.
  ## the output is therefore a matrix of size nrow(A) by nrow(B)
  
  ## see links for more information:
  ## https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
  ## https://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
  
  
  an <- apply(A, 1, function(rvec) crossprod(rvec, rvec))
  bn <- apply(B, 1, function(rvec) crossprod(rvec, rvec))
  m <- nrow(A)
  n <- nrow(B)
  tmp <- matrix(rep(an, n), nrow = m) 
  tmp <- tmp + matrix(rep(bn, m), nrow = m, byrow = TRUE)
  
  tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
}



factor2simil <- function(x, nms, distance = FALSE) {
  res <- tcrossprod(model.matrix(~ 0 + x))
  dimnames(res) <- list(nms, nms)
  if (distance) res <- 1 - res
  res
}


dummy2factor <- function(x) {
  nms <- colnames(x)
  as.factor(apply(x, 1, function(x_i) nms[which.max(x_i)]))
}


create_rsa_models <- function(conditions, rsa_model_names) {
  
  n_condition <- nrow(conditions)
  n_models <- length(rsa_model_names)
  condition_names <- conditions$id

  ## get rsa model structures
  rsa_models <- lapply(
    rsa_model_names,
    \(x, nms = condition_names, distance = TRUE) factor2simil(x = conditions[[x]], nms, distance)
  )

  ## put in array
  size <- c(n_condition, n_condition, n_models)
  rsa_models <- array(unlist(rsa_models), dim = size)
  dimnames(rsa_models) <- list(condition_names, condition_names, rsa_model_names)

  rsa_models

}


fit_rsa <- function(Y_test, fit_train, params) {

  b_train <- coef(fit_train)  ## class means
  d2 <- t(pdist(Y_test, b_train))  ## single-trial squared euclidean distances

  .fit_singletrial_model(
    models = params$rsa_models,
    y = t(d2),
    trial_condition = params$trial_condition,
    fit_fun = \(model_i, y_i) {
      t(coef(.lm.fit(cbind(1, model_i), t(y_i))))
    }
  )

}


create_encoding_contrasts <- function(A, encoding_contrasts) {
  contrast_mats <- lapply(
    encoding_contrasts,
    \(x, a = A) {
      m <- a * matrix(rep(colnames(a) %in% x), nrow(a), nrow = nrow(a), byrow = TRUE)
      t(apply(m, 1, signed_scale))
    }
  )
  ## convert to array:
  size <- c(nrow(A), ncol(A), length(encoding_contrasts))
  nms <- list(rownames(A), colnames(A), names(encoding_contrasts))
  contrasts <- array(unlist(contrast_mats), dim = size, dimnames = nms)
  contrasts
}


signed_scale <- function(x) {
    is_pos <- x > 0
    is_neg <- x < 0
    x[is_pos] <- x[is_pos] / sum(x[is_pos])
    x[is_neg] <- x[is_neg] / -sum(x[is_neg])
    x
}


fit_encoding <- function(Y_test, fit_train, params) {

  contrast <- t(solve(params$A[, 1:nrow(params$A)]))
  b_axes_train <- crossprod(coef(fit_train), contrast)
  encoding_projs <- Y_test %*% b_axes_train

  .fit_singletrial_model(
    models = params$encoding_contrasts,
    y = encoding_projs,
    trial_condition = params$trial_condition,
    fit_fun = \(model_i, y_i) y_i %*% model_i
  )

}


.fit_singletrial_model <- function(models, y, trial_condition, fit_fun) {

  cond_names <- rownames(models)  ## unique trial conditions
  model_names <- dimnames(models)[[3]]
  coefs <- vector("list", length(cond_names))
  for (i in seq_along(cond_names)) {
    cond <- cond_names[i]
    y_i <- y[trial_condition == cond, ]  ## select trials belonging to condition
    model_i <- models[cond, , ]
    coefs[[i]] <- fit_fun(model_i, y_i)  ## returns coefs
  }

  coefs <- as.data.frame(do.call(rbind, coefs))
  nms <- model_names
  if (ncol(coefs) > length(model_names)) nms <- c("intercept", nms)
  names(coefs) <- nms
  coefs$condition <- trial_condition
  
  coefs

}

summarize_df <- function(df, fun = \(x) colMeans2(x) / colSds(x)) {
  m <- as.matrix(df[vapply(df, is.numeric, logical(1))])
  c(fun(m))
}


.simulate_singletrial <- function(params) {

    ## initial encoding models (for both decoding+RSA and encoding)

    Y <- list(train = generate_data(params), test = generate_data(params))
    fit <- lapply(Y, \(y, x = params$X) fit_encoder(X = x, Y = y))

    ## decoding + rsa

    rsa_coefs <- fit_rsa(Y$test, fit$train, params)

    ## encoding

    encoding_projs <- fit_encoding(Y$test, fit$train, params)

    ## summarize

    rsa_coefs_sum <- summarize_df(rsa_coefs)
    encoding_projs_sum <- summarize_df(encoding_projs)
    res <- c(rsa_coefs_sum, encoding_projs_sum)
    method <- c(rep("rsa", length(rsa_coefs_sum)), rep("encoding", length(encoding_projs_sum)))
    data.table(effsize = res, term = names(res), method = method)

}


simulate_singletrial <- function(params) {
  
  sim_idx <- seq_len(params$n["sim"])
  res <- foreach(
    x = sim_idx, 
    .options.future = list(seed = TRUE)
    ) %dofuture% .simulate_singletrial(params)
  
  rbindlist(res, idcol = "sim_idx")

}
