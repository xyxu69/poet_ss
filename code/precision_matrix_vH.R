# Precision-matrix simulation.
for (.pkg in c("mvtnorm", "SpatialNP", "ICSNP", "glassoFast", "pcaPP")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) stop("Missing R package: ", .pkg)
}
rm(.pkg)

precision_validate <- function(Y, m) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)) || nrow(Y) < 4L || ncol(Y) < 2L)
    stop("Y must be a finite numeric matrix with n>=4 and p>=2.")
  if (length(m) != 1L || !is.finite(m) || m != as.integer(m) || m < 1L || m >= min(dim(Y)))
    stop("m must be an integer between 1 and min(n,p)-1.")
  Y
}

precision_glasso_penalty <- "all_entries_including_diagonal"

precision_glasso <- function(S, lambda) {
  S <- (S + t(S))/2
  if (any(!is.finite(S)) || any(diag(S) <= 0) || length(lambda)!=1L || !is.finite(lambda) || lambda <= 0)
    stop("Invalid GLASSO residual or penalty.")
  if (max(abs(S[row(S)!=col(S)])) <= lambda)
    return(diag(1/(diag(S)+lambda),nrow=nrow(S)))
  fit <- glassoFast::glassoFast(S, lambda)
  if (!identical(as.integer(fit$errflag), 0L)) stop("GLASSO failed, errflag=", fit$errflag)
  if (fit$niter>=10000L) stop("GLASSO reached its iteration limit.")
  V <- (fit$wi + t(fit$wi))/2
  if (any(!is.finite(V))) stop("Nonfinite GLASSO precision.")
  chol(V)
  V
}

precision_cv_score_mode <- function() "legacy_log_det"

precision_errors <- function(V, Vu, target, target_u, prefix="") {
  p <- nrow(V); E <- V-target; Eu <- Vu-target_u
  ans <- c(V_F_err=norm(E,"F"), V_u_F_err=norm(Eu,"F"),
           V_normalized_F_err=norm(E,"F")/sqrt(p),
           V_u_normalized_F_err=norm(Eu,"F")/sqrt(p),
           V_max_err=max(abs(E)), V_u_max_err=max(abs(Eu)),
           V_spectral_err=norm(E,"2"), V_u_spectral_err=norm(Eu,"2"))
  names(ans) <- paste0(prefix,names(ans)); ans
}

precision_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path,".tmp_",Sys.getpid())
  saveRDS(object,tmp)
  if (file.exists(path) || !file.rename(tmp,path)) stop("Checkpoint destination exists or rename failed: ",path)
}

# Population sign targets.
poet_population_targets <- function(Sigma0, Sigma_u0, m, rel.tol = 1e-8) {
  p <- nrow(Sigma0)
  stopifnot(is.matrix(Sigma0), is.matrix(Sigma_u0),
            identical(dim(Sigma0), dim(Sigma_u0)), p > 2L,
            length(m) == 1L, is.finite(m), m == as.integer(m), m >= 1L, m < p)
  eig <- eigen((Sigma0 + t(Sigma0)) / 2, symmetric = TRUE)
  lambda <- eig$values
  if (any(!is.finite(lambda)) || min(lambda) <= 0) stop("Sigma0 must be positive definite.")
  # Substitute x=p*t so integrals retain a well-resolved scale at large p.
  laplace <- function(x) vapply(x, function(xx) {
    if (!is.finite(xx)) return(0)
    exp(-0.5 * sum(log1p(2 * xx * lambda / p)))
  }, numeric(1))
  c_d <- integrate(laplace, 0, Inf, rel.tol = rel.tol,
                   subdivisions = 1000L, stop.on.error = TRUE)$value
  theta <- vapply(seq_len(m), function(j) {
    lambda[j] * integrate(function(x) laplace(x) / (1 + 2 * x * lambda[j] / p),
                          0, Inf, rel.tol = rel.tol, subdivisions = 1000L,
                          stop.on.error = TRUE)$value
  }, numeric(1))
  Gamma_m <- eig$vectors[, seq_len(m), drop = FALSE]
  Lambda_m <- diag(lambda[seq_len(m)], nrow = m, ncol = m)
  Lambda_sign <- diag(theta, nrow = m, ncol = m)
  H0u <- c_d * Sigma_u0
  H0 <- Gamma_m %*% Lambda_sign %*% t(Gamma_m) + H0u
  H0 <- (H0 + t(H0)) / 2
  list(c_d = c_d, Gamma_m = Gamma_m, Lambda_m = Lambda_m,
       Lambda_sign = Lambda_sign, theta = theta, H0 = H0, H0u = H0u,
       P0 = chol2inv(chol(H0)), Vu_sign = chol2inv(chol(H0u)),
       Sigma0 = Sigma0, Sigma_u0 = Sigma_u0)
}

# Project a matrix to the PSD cone
eigen_projection <- function(A, epsilon = 1e-5) {
  eig <- eigen(A, symmetric = TRUE)
  values <- pmax(eig$values, epsilon)
  return(eig$vectors %*% diag(values, nrow=length(values)) %*% t(eig$vectors))
}


# Select gamma
select_gamma <- function(A, mu) {
  if (length(mu)!=1L || !is.finite(mu) || mu<=0 || any(!is.finite(A))) stop("Invalid smoothing input.")
  a <- sort(as.vector(abs(A))/mu, decreasing=TRUE)
  if (!length(a)) stop("Empty smoothing input.")
  if (sum(a) <= 1) return(list(gamma=0,u=length(a),threshold=0))
  cs <- cumsum(a); j <- seq_along(a)
  u <- max(which(a > (cs-1)/j))
  gamma <- (cs[u]-1)/u
  list(gamma=gamma,u=u,threshold=mu*gamma)
}
# smoothed l_inf norm approximation
smooth_linf <- function(A, mu) {
  gamma_obj <- select_gamma(A, mu)
  gamma <- gamma_obj$gamma
  U_hat <- sign(A) * pmax(abs(A) / mu - gamma, 0)
  return(sum(U_hat * A) - mu * norm(U_hat, "F")^2 / 2)
}

accelerated_proximal_gradient <- function(S_hat, mu, max_iter = 300, tol = 1e-4) {
  p <- nrow(S_hat)
  S0 <- eigen_projection(S_hat)
  S_prev <- S0
  W_prev <- S0     
  for (t in 1:max_iter) {
    theta_t <- 2 / (1 + t)   
    M_t <- (1 - theta_t) * S_prev + theta_t * W_prev   
    Delta <- S_hat - M_t
    gamma <- select_gamma(Delta, mu)$gamma
    G_t <- -sign(Delta) * pmax(abs(Delta / mu) - gamma, 0)
    eta <- mu
    step_size <- eta / theta_t
    W_temp <- eigen_projection(W_prev - step_size * G_t)
    S_temp <- (1 - theta_t) * S_prev + theta_t * W_temp 
    W_prev <- W_temp
    S_curr <- S_temp    
    gap <- abs(smooth_linf(S_hat - S_curr, mu) - smooth_linf(S_hat - S_prev, mu))
    if (gap <= tol * mu) {
      break
    }    
    S_prev <- S_curr
  }
  return(S_curr)
}



generate_precision_factor_data_joint <- function(B, 
                                                 n = 250, 
                                                 mu = rep(0, nrow(B)),
                                                 mi = 1) {
  p <- nrow(B)
  m <- ncol(B)
  if (length(mu)!=p || any(!is.finite(B)) || any(!is.finite(mu)) || n<4L || m<1L || p<=m) stop("Invalid simulation design.")
  
  # idiosyncratic precision matrix Vu = (0.4^|i-j|)
  indices <- matrix(1:p, nrow = p, ncol = p)
  V_u_raw <- 0.4^abs(indices - t(indices))
  Sigma_u_raw <- solve(V_u_raw)
  
  # covariance of y before normalization
  Sigma_raw <- B %*% t(B) + Sigma_u_raw
  c0 <- p / sum(diag(Sigma_raw))
  
  Sigma0 <- c0 * Sigma_raw
  Sigma_u0 <- c0 * Sigma_u_raw
  
  # covariance of (f, u)
  Sigma_fu <- matrix(0, nrow = m + p, ncol = m + p)
  Sigma_fu[1:m, 1:m] <- c0 * diag(m)
  Sigma_fu[(m + 1):(m + p), (m + 1):(m + p)] <- c0 * Sigma_u_raw
  
  if (mi == 1) {
    # Scenario I: N(0, Sigma_fu)
    joint_sample <- mvtnorm::rmvnorm(
      n, mean = rep(0, m + p), sigma = Sigma_fu
    )
    
  } else if (mi == 2) {
    # Scenario II: t_4
    df <- 4
    joint_sample <- mvtnorm::rmvt(
      n, sigma = Sigma_fu * (df - 2) / df, df = df
    )
    
  } else if (mi == 3) {
    # Scenario III: t_2.2
    df <- 2.2
    joint_sample <- mvtnorm::rmvt(
      n, sigma = Sigma_fu * (df - 2) / df, df = df
    )
    
  } else if (mi == 4) {
    # Scenario IV: 0.8N(0,Sigma_fu)+0.2N(0,10Sigma_fu)
    z <- rbinom(n, size = 1, prob = 0.2)
    joint_sample <- matrix(0, nrow = n, ncol = m + p)
    
    n1 <- sum(z == 0)
    n2 <- sum(z == 1)
    
    if (n1 > 0) {
      joint_sample[z == 0, ] <- mvtnorm::rmvnorm(
        n1, mean = rep(0, m + p), sigma = Sigma_fu
      )
    }
    if (n2 > 0) {
      joint_sample[z == 1, ] <- mvtnorm::rmvnorm(
        n2, mean = rep(0, m + p), sigma = 10 * Sigma_fu
      )
    }
    
    joint_sample <- joint_sample / sqrt(2.8)
    
  } else {
    stop("mi must be 1, 2, 3, or 4.")
  }
  
  f_t <- joint_sample[, 1:m, drop = FALSE]
  u_t <- joint_sample[, (m + 1):(m + p), drop = FALSE]
  Y <- f_t %*% t(B) + u_t
  
  if (length(mu) == p) {
    Y <- Y + matrix(mu, nrow = n, ncol = p, byrow = TRUE)
  }
  
  list(
    data = Y,
    Sigma0 = Sigma0,
    Sigma_u0 = Sigma_u0,
    V0 = solve(Sigma0),
    V_u0 = solve(Sigma_u0),
    V_u_raw = V_u_raw, c0 = c0, scenario = mi,
    scale_convention = "t_and_mixture_covariance_standardized"
  )
}
robust_variance_diag <- function(Y, delta0 = 1e-4) {
  n <- nrow(Y)
  p <- ncol(Y)
  if (any(!is.finite(Y)) || n < 3 || p < 2) stop("Invalid robust variance input.")
  eps <- 1 / max(n, p)^2
  if (n <= 2 * log(1 / eps)) stop("FLW Catoni tuning requires n > 4*log(max(n,p)).")
  h_func <- function(x) sign(x) * log1p(abs(x) + x^2 / 2)
  catoni <- function(z, v) {
    if (!is.finite(v) || v < 0) stop("Invalid Catoni variance bound.")
    if (v == 0 || min(z) == max(z)) return(mean(z))
    alpha <- sqrt(2 * log(1 / eps) /
      (n * (v + 2 * v * log(1 / eps) / (n - 2 * log(1 / eps)))))
    uniroot(function(mu) sum(h_func(alpha * (z - mu))), range(z),
            tol = .Machine$double.eps^0.75 * max(1, max(abs(z))), maxiter = 1000)$root
  }
  v_max <- max(apply(Y, 2, stats::var))
  v2_max <- max(apply(Y^2, 2, stats::var))
  sigma2 <- vapply(seq_len(p), function(j) {
    mu <- catoni(Y[, j], v_max)
    max(catoni(Y[, j]^2, v2_max) - mu^2, delta0)
  }, numeric(1))
  diag(sqrt(sigma2), nrow = p, ncol = p)
}
poet_kendall_tau_a <- function(Y) {
  if (!requireNamespace("pcaPP", quietly = TRUE)) {
    stop("Install pcaPP before running: install.packages('pcaPP')")
  }
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L || any(!is.finite(Y))) {
    stop("Kendall input must be a finite numeric matrix with at least two rows.")
  }
  n <- nrow(Y)
  p <- ncol(Y)
  pairs <- as.double(n) * (n - 1) / 2
  # The old sign-crossproduct computes tau-a, whereas cor.fk computes tau-b.
  # q[j] is the fraction of pairs NOT tied in column j.
  q <- vapply(seq_len(p), function(j) {
    counts <- as.double(rle(sort(Y[, j], method = "radix"))$lengths)
    1 - sum(counts * (counts - 1) / 2) / pairs
  }, numeric(1))
  active <- which(q > 0)
  result <- matrix(0, p, p)
  if (length(active) >= 2L) {
    tau_b <- pcaPP::cor.fk(Y[, active, drop = FALSE])
    result[active, active] <- tau_b * tcrossprod(sqrt(q[active]))
  }
  # The original code explicitly sets ALL diagonal entries to 1,
  # including those belonging to constant columns.
  diag(result) <- 1
  result
}

estimate_Lambda <- function(Y, D_hat, m) {
  T_hat <- poet_kendall_tau_a(Y)
  R_hat <- sin(pi/2*T_hat)
  Sigma1_hat <- D_hat %*% R_hat %*% D_hat
  eig_vals <- eigen(Sigma1_hat,symmetric=TRUE,only.values=TRUE)$values
  list(Lambda_hat=diag(eig_vals[seq_len(m)],nrow=m),Sigma1_hat=Sigma1_hat)
}
poet_huber_location <- function(z, H) {
  z <- as.numeric(z)
  if (!length(z) || any(!is.finite(z)) || length(H) != 1L || !is.finite(H) || H < 0) {
    stop("Huber location requires finite observations and a nonnegative threshold.")
  }
  limits <- range(z)
  if (limits[1] == limits[2]) return(limits[1])
  if (H == 0) return(stats::median(z))
  score <- function(mu) sum(pmin(H, pmax(-H, z - mu)))
  uniroot(score, limits, tol = .Machine$double.eps^0.75 * max(1, max(abs(z))),
          maxiter = 1000L)$root
}

estimate_mu_huber_cv <- function(Y, c_candidates = c(0.5, 1.0, 1.5), nfold = 3) {
  n <- nrow(Y)
  p <- ncol(Y)
  if (n < nfold || nfold < 2 || p < 2 || any(!is.finite(Y))) stop("Invalid Huber CV input.")
  mu_hat <- best_c <- numeric(p)
  for (i in seq_len(p)) {
    y <- Y[, i]
    scale_full <- stats::mad(y)
    # Retain original per-coordinate random folds and candidate grid.
    fold_id <- sample(rep(seq_len(nfold), length.out = n))
    cv_loss <- vapply(c_candidates, function(c_val) {
      H <- c_val * scale_full * sqrt(n / log(p))
      mean(vapply(seq_len(nfold), function(k) {
        mu_train <- poet_huber_location(y[fold_id != k], H)
        mean((y[fold_id == k] - mu_train)^2)
      }, numeric(1)))
    }, numeric(1))
    best_c[i] <- c_candidates[which.min(cv_loss)]
    mu_hat[i] <- poet_huber_location(y, best_c[i] * scale_full * sqrt(n / log(p)))
  }
  list(mu = mu_hat, c_opt = best_c)
}

symmetrize_Y <- function(Y) {
  n_pairs <- floor(nrow(Y) / 2)
  if (n_pairs < 2L) stop("RegTME needs at least four original observations.")
  Y[seq(1L, 2L * n_pairs, 2L), , drop = FALSE] -
    Y[seq(2L, 2L * n_pairs, 2L), , drop = FALSE]
}

regTME_svd <- function(Y, alpha = NULL, tol = 1e-4, maxit = 100, eps = 1e-12) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  if (n < 2L || p < 2L || any(!is.finite(Y))) stop("Invalid RegTME input.")
  # n is the actual number of paired differences, not the original sample size.
  gamma <- p / n
  S_sample <- stats::cov(Y)
  sample_trace <- sum(diag(S_sample))
  if (!is.finite(sample_trace) || sample_trace <= 0) stop("RegTME data have zero variation.")
  S_sample_norm <- p * S_sample / sample_trace
  s_max <- norm(S_sample_norm, type = "2")
  if (is.null(alpha)) alpha <- max(0.1, 1.1 * (gamma - 1 + s_max * (1 + sqrt(gamma))^2))
  if (!is.finite(alpha) || alpha <= max(0, gamma - 1)) stop("RegTME alpha must exceed max(0,p/n-1).")
  use_proj <- n < p
  if (use_proj) {
    decomp <- svd(Y, nu = n, nv = n)
    Z <- sweep(decomp$u, 2, decomp$d, "*")
    proj_mat <- decomp$v
  } else {
    Z <- Y
    proj_mat <- diag(p)
  }
  dim_Z <- ncol(Z)
  shrink <- alpha / (1 + alpha)
  Sigma <- shrink * diag(dim_Z)
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    Sigma_inv <- chol2inv(chol((Sigma + t(Sigma)) / 2))
    quad_forms <- rowSums((Z %*% Sigma_inv) * Z)
    if (any(!is.finite(quad_forms)) || any(quad_forms <= 0)) stop("RegTME has a nonpositive directional denominator.")
    weighted <- Z / sqrt(quad_forms)
    # The ambient dimension p is required even when computations use the row span.
    Sigma_new <- (p / n) * crossprod(weighted) / (1 + alpha) + shrink * diag(dim_Z)
    if (norm(Sigma_new - Sigma, type = "F") < tol) {
      Sigma <- Sigma_new
      converged <- TRUE
      break
    }
    Sigma <- Sigma_new
  }
  if (!converged) warning("RegTME did not reach its tolerance within maxit.")
  A_low <- Sigma - shrink * diag(dim_Z)
  A_full <- if (use_proj) proj_mat %*% A_low %*% t(proj_mat) else A_low
  A_full <- (A_full + t(A_full)) / 2
  den <- sum(diag(A_full))
  if (!is.finite(den) || den <= eps) stop("RegTME unshrunk shape has nonpositive trace.")
  ans <- p * A_full / den
  attr(ans, "regtme_diagnostics") <- list(alpha = alpha, gamma = gamma, iterations = iter,
                                          converged = converged, ambient_dimension = p)
  ans
}

poet_ipsn_shape <- function(Y, m, mu_hat) {
  eig <- eigen(SpatialNP::SSCov(Y), symmetric = TRUE)
  Gamma_pilot <- eig$vectors[, seq_len(m), drop = FALSE]
  centered <- sweep(Y, 2, mu_hat, "-")
  orthogonal <- centered - (centered %*% Gamma_pilot) %*% t(Gamma_pilot)
  radius <- sqrt(rowSums(orthogonal^2))
  if (any(!is.finite(radius)) || any(radius <= 0)) stop("IPSN projected radius is zero or nonfinite.")
  normalized <- sqrt(ncol(Y)) * centered / radius
  raw_shape <- crossprod(normalized) / nrow(Y)
  ncol(Y) * raw_shape / sum(diag(raw_shape))
}

compute_2_part_nonsparse <- function(Y, m, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  Y <- precision_validate(Y, m)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop=FALSE] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_tilde,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c 
    mu <- sqrt(log(p) / n) +1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  
  if (method == "FLW") {
    D_hat <- robust_variance_diag(Y)
    estimate <- estimate_Lambda(Y, D_hat, m)
    Lambda_hat <- estimate$Lambda_hat
    Sigma_hat <- estimate$Sigma1_hat
    HK<-SpatialNP::SSCov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop=FALSE] 
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    R_tilde1 <- eigen_projection(R_tilde,epsilon = 1e-5)
    max_diff <- max(abs(R_tilde - R_tilde1))
    mu <- sqrt(log(p) / n)
    if (max_diff < 2*mu) {
      R_tilde2 <- R_tilde1
    } else {
      R_tilde2 <- accelerated_proximal_gradient(R_tilde, mu, max_iter = 300, tol = 5e-1)
    }
    Sigma_u_hat <- R_tilde2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <- Lambda_hat*c
    mu <- sqrt(log(p) / n) +1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "SPCA") {
    mu_hat <- ICSNP::spatial.median(Y)
    Sigma_hat <- SpatialNP::SCov(Y, location=mu_hat)*p
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop=FALSE] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_tilde,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c
    mu <- sqrt(log(p) / n) + sqrt(log(n) /n)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "TylerSPCA") {
    estimate  <- compute_2_part_nonsparse(Y, m, method = "SPCA")
    Sigma_u_hat0 <- estimate$Sigma_u_hat0
    Gamma_hat    <- estimate$Gamma_hat
    Lambda_hat0  <- estimate$Lambda_hat0
    mu <- sqrt(log(p) / n) + sqrt(log(n) /n)
    lambda_opt0 <- 0.05 * mu
    Omega_u_hat0 <- precision_glasso(Sigma_u_hat0, lambda_opt0)
    Omega_u_hat0 <- (Omega_u_hat0 + t(Omega_u_hat0)) / 2 
    tmp0 <- solve(solve(Lambda_hat0) + t(Gamma_hat) %*% Omega_u_hat0 %*% Gamma_hat)
    Omega_hat0 <- Omega_u_hat0 - Omega_u_hat0 %*% Gamma_hat %*% tmp0 %*% t(Gamma_hat) %*% Omega_u_hat0
    mu_hat <- ICSNP::spatial.median(Y)
    Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")
    Sigma_S_inv <- Omega_hat0
    quad_forms <- rowSums((Y_centered %*% Sigma_S_inv) * Y_centered) 
    if (any(quad_forms<=0) || any(!is.finite(quad_forms))) stop("Invalid or zero Tyler pilot quadratic forms.")
    weights <- 1/quad_forms
    Sigma_T <- (p / n) * crossprod(Y_centered*sqrt(weights))
    Sigma_hat <- (Sigma_T + t(Sigma_T)) / 2  
    eig_Sigma_T <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma_T$vectors[, 1:m, drop = FALSE]
    if(length(eig_Sigma_T$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig_Sigma_T$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig_Sigma_T$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_hat_T <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_hat_T,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "IPSN") {
    mu_hat <- estimate_mu_huber_cv(Y)$mu
    Sigma_0_hat <- poet_ipsn_shape(Y, m, mu_hat)
    eig_Sigma0 <- eigen(Sigma_0_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma0$vectors[, 1:m, drop=FALSE] 
    if(length(eig_Sigma0$values[1:m]) == 1) {
      Lambda_hat0 <- matrix(eig_Sigma0$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eig_Sigma0$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    R_tilde <- Sigma_0_hat - low_rank_part
    Sigma_u_hat0 <- eigen_projection(R_tilde, epsilon = 1e-3)
    Sigma_u_hat0 <- (Sigma_u_hat0 + t(Sigma_u_hat0)) / 2
    mu <- sqrt(log(p) / n) + sqrt(log(p) /p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "RegTME") {
    Y_sym <- symmetrize_Y(Y)
    Sigma_tme <- regTME_svd(Y_sym, alpha = NULL)
    eig <- eigen(Sigma_tme, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    eigvals <- eig$values[1:m]
    eigvals <- pmax(eigvals, 0)
    if (m == 1) {
      Lambda_hat0 <- matrix(eigvals, nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eigvals)
    }
    n_sym <- nrow(Y_sym)
    Sigma_u_hat0 <- Sigma_tme - Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    Sigma_u_hat0 <- eigen_projection(Sigma_u_hat0,epsilon = 1e-3)
    Sigma_u_hat0 <- (Sigma_u_hat0 + t(Sigma_u_hat0)) / 2
    mu <- sqrt(log(p) / n_sym)+1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  return(list(
    cv_grid = cv_grid,
    Lambda_hat0 = Lambda_hat0,
    Gamma_hat = Gamma_hat,
    Sigma_u_hat0 = Sigma_u_hat0
  ))
}
cv_glasso <- function(Y, lambda_grid, K=2L, m=3L,
                      method=c("SAMPLE","SPCA","TylerSPCA","FLW","IPSN","RegTME"), folds=NULL,
                      score_mode=precision_cv_score_mode()) {
  score_mode <- match.arg(score_mode,"legacy_log_det")
  method <- match.arg(method); n <- nrow(Y)
  if (is.null(folds)) folds <- sample(rep(seq_len(K),length.out=n))
  if (length(folds)!=n || anyNA(folds) || !all(seq_len(K)%in%folds)) stop("Invalid CV folds.")
  scores <- matrix(NA_real_,K,length(lambda_grid))
  for (k in seq_len(K)) {
    train <- compute_2_part_nonsparse(Y[folds!=k,,drop=FALSE],m,method)
    valid <- compute_2_part_nonsparse(Y[folds==k,,drop=FALSE],m,method)
    for (j in seq_along(lambda_grid)) {
      V <- precision_glasso(train$Sigma_u_hat0,lambda_grid[j])
      scores[k,j] <- log(det(V))-sum(V*valid$Sigma_u_hat0)
    }
  }
  means <- colMeans(scores)
  selected <- which.max(means)
  if (!length(selected)) stop("No GLASSO tuning score is available.")
  list(lambda_opt=lambda_grid[selected], scores=means, fold_scores=scores,
       selected_index=selected, lambda_grid=lambda_grid, score_mode=score_mode,
       glasso_penalty=precision_glasso_penalty)
}

evaluate_precision_estimators <- function(Y, Sigma_true, Sigma_u_true,
  method=c("SAMPLE","SPCA","TylerSPCA","FLW","IPSN","RegTME"),
  estimate_2_part=NULL, m=3L, folds=NULL, population_targets=NULL,
  score_mode=precision_cv_score_mode()) {
  method <- match.arg(method); p <- ncol(Y)
  if (is.null(estimate_2_part)) estimate_2_part <- compute_2_part_nonsparse(Y,m,method)
  fit <- estimate_2_part
  cv <- cv_glasso(Y,fit$cv_grid,K=2L,m=m,method=method,folds=folds,score_mode=score_mode)
  Vu <- precision_glasso(fit$Sigma_u_hat0,cv$lambda_opt)
  G <- fit$Gamma_hat; L <- fit$Lambda_hat0
  VG <- Vu %*% G
  V <- Vu - VG %*% solve(solve(L)+crossprod(G,VG), t(VG))
  V <- (V+t(V))/2
  truth <- solve(Sigma_true); truth_u <- solve(Sigma_u_true)
  losses <- precision_errors(V,Vu,truth,truth_u)
  if (method=="SPCA") {
    # SS's direct population target is hybrid H0, not Sigma0.
    if (is.null(population_targets)) population_targets <- poet_population_targets(Sigma_true,Sigma_u_true,m)
    names(losses) <- paste0("cross_target_",names(losses))
    losses <- c(losses,precision_errors(V,Vu,population_targets$P0,
                       population_targets$Vu_sign,"matched_target_"))
  }
  answer <- as.list(c(losses,lambda_selected=cv$lambda_opt))
  attr(answer,"cv_diagnostics") <- cv[c("score_mode","glasso_penalty","lambda_grid","fold_scores","scores","selected_index")]
  answer
}

run_precision_simulation_onecase <- function(p=500L,sims=100L,mi=1L,
  n=250L,m=3L,checkpoint_dir=NULL,seed=NULL,verbose=TRUE,loading_sd=NULL,
  score_mode=precision_cv_score_mode()) {
  score_mode <- match.arg(score_mode,"legacy_log_det")
  if (p<=m || n<max(20L,4L*m) || sims<1L || !mi%in%1:4)
    stop("Use p>m, n>=max(20,4m), sims>=1, mi in 1:4; default paper design has m=3.")
  if (is.null(loading_sd)) {
    if (m>3L) stop("For m>3, supply loading_sd explicitly; the paper specifies only 3 factors.")
    loading_sd <- c(1,.75,.5)[seq_len(m)]
  }
  if (length(loading_sd)!=m || any(!is.finite(loading_sd)) || any(loading_sd<=0)) stop("Invalid loading_sd.")
  if (!is.null(seed)) set.seed(seed)
  methods <- c("SAMPLE","SPCA","TylerSPCA","FLW","IPSN","RegTME")
  settings <- list(p=p,n=n,m=m,sims=sims,mi=mi,seed=seed,loading_sd=loading_sd,
    target="trace-p scatter precision; SS hybrid target separately",
    glasso_penalty=precision_glasso_penalty,cv_score_mode=score_mode,
    pilot="FLW/IPSN SpatialNP::SSCov pilot",
    scale_convention="t_and_mixture_covariance_standardized")
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
    precision_atomic_rds(settings,file.path(checkpoint_dir,"settings.rds"))
  }
  simulation_results <- vector("list",sims)
  for (s in seq_len(sims)) {
    B <- sweep(matrix(rnorm(p*m),p,m),2,loading_sd,"*")
    dat <- generate_precision_factor_data_joint(B,n=n,mu=rep(0,p),mi=mi)
    Y <- dat$data
    pop <- poet_population_targets(dat$Sigma0,dat$Sigma_u0,m)
    ans <- setNames(vector("list",length(methods)),methods)
    for (method in methods) {
      started <- proc.time()[["elapsed"]]
      fit <- compute_2_part_nonsparse(Y,m,method)
      ans[[method]] <- evaluate_precision_estimators(Y,dat$Sigma0,dat$Sigma_u0,
        method,fit,m=m,population_targets=pop,score_mode=score_mode)
      ans[[method]]$runtime_sec <- unname(proc.time()[["elapsed"]] - started)
    }
    simulation_results[[s]] <- ans
    if (!is.null(checkpoint_dir)) precision_atomic_rds(
      list(replication=s,loss=ans,settings=settings,random_seed=.Random.seed),
      file.path(checkpoint_dir,sprintf("rep_%04d.rds",s)))
    if (verbose) cat(sprintf("Precision scenario %d: %d/%d complete%s\n",mi,s,sims,
      if (is.null(checkpoint_dir)) "" else " and saved"))
  }
  means <- sds <- setNames(vector("list",length(methods)),methods)
  for (method in methods) {
    mat <- do.call(rbind,lapply(simulation_results,function(x) unlist(x[[method]])))
    means[[method]] <- as.list(setNames(colMeans(mat),paste0("mean_",colnames(mat))))
    sds[[method]] <- as.list(setNames(apply(mat,2,sd),paste0("sd_",colnames(mat))))
  }
  list(simulation_results=simulation_results,mean_results=means,sd_results=sds,settings=settings)
}

if (isTRUE(getOption("poet.autorun", TRUE))) {
  # 运行参数：只跑场景 IV 时，把 mi_vec 改成 4。
  n <- 250L
  p <- 500L
  m <- 3L
  sims <- 100L
  mi_vec <- 1:4                 # 1=正态，2=t4，3=t2.2，4=混合正态
  seed <- 998L

  .poet_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(.poet_file)) {
    .poet_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(.poet_arg) && basename(sub("^--file=", "", .poet_arg[1])) == "precision_matrix_vH.R")
      .poet_file <- sub("^--file=", "", .poet_arg[1])
  }
  .poet_dir <- if (length(.poet_file) == 1L && !is.na(.poet_file) && file.exists(.poet_file))
    dirname(normalizePath(.poet_file, winslash = "/")) else getwd()
  output_dir <- file.path(.poet_dir, "results", "precision")

  if (!length(mi_vec) || anyNA(mi_vec) || any(!mi_vec %in% 1:4) || anyDuplicated(mi_vec))
    stop("mi_vec must contain distinct scenario numbers from 1 to 4.")
  .poet_numbers <- c(n, p, m, sims, seed)
  if (any(!is.finite(.poet_numbers)) || any(.poet_numbers != floor(.poet_numbers)) ||
      n < 20L || p < 4L || m < 1L || m > 3L || m >= min(n, p) || sims < 1L ||
      seed < 0L || seed > .Machine$integer.max) stop("Invalid simulation parameters.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  .poet_index <- 1L
  repeat {
    run_dir <- file.path(output_dir, sprintf("precision_run_%03d", .poet_index))
    if (dir.create(run_dir, showWarnings = FALSE)) break
    if (!file.exists(run_dir)) stop("Cannot create output directory: ", run_dir)
    .poet_index <- .poet_index + 1L
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  settings <- list(n = n, p = p, m = m, sims = sims, scenarios = mi_vec, seed = seed)
  saveRDS(list(settings = settings, session_info = sessionInfo()),
          file.path(run_dir, "run_info.rds"))
  cat("OUTPUT_DIR: ", run_dir, "\n", sep = "")

  results <- tables <- list()
  for (mi in mi_vec) {
    result <- run_precision_simulation_onecase(p = p, n = n, m = m, sims = sims, mi = mi, seed = seed,
                checkpoint_dir = file.path(run_dir, paste0("replicates_mi", mi)))
    results[[paste0("mi", mi)]] <- result
    table <- do.call(rbind, lapply(names(result$mean_results), function(method) {
      means <- unlist(result$mean_results[[method]])
      metrics <- sub("^mean_", "", names(means))
      deviations <- unlist(result$sd_results[[method]])[paste0("sd_", metrics)]
      data.frame(scenario = mi, method = method, metric = metrics,
                 mean = unname(means), sd = unname(deviations))
    }))
    table$value <- sprintf("%.4f(%.4f)", table$mean, table$sd)
    rownames(table) <- NULL
    tables[[paste0("mi", mi)]] <- table
    .poet_path <- file.path(run_dir, paste0("precision_mi", mi))
    saveRDS(list(result = result, table = table), paste0(.poet_path, ".rds"))
    tryCatch(write.csv(table, paste0(.poet_path, ".csv"), row.names = FALSE),
             error = function(e) warning("RDS saved; CSV export failed: ", conditionMessage(e)))
  }
  summary_table <- do.call(rbind, tables)
  rownames(summary_table) <- NULL
  .poet_path <- file.path(run_dir, "precision_all_scenarios")
  saveRDS(list(result = results, table = summary_table), paste0(.poet_path, ".rds"))
  tryCatch(write.csv(summary_table, paste0(.poet_path, ".csv"), row.names = FALSE),
           error = function(e) warning("RDS saved; CSV export failed: ", conditionMessage(e)))
  cat("SAVED: ", run_dir, "\n", sep = "")
}
