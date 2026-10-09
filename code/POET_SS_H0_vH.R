# Source this file to run; options(poet.autorun=FALSE) loads functions only.
for (.poet_package in c("mvtnorm", "SpatialNP", "ICSNP")) {
  if (!requireNamespace(.poet_package, quietly = TRUE)) stop("Missing R package: ", .poet_package)
}
rm(.poet_package)

poet_validate_matrix_input <- function(Y, m) {
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 4L || ncol(Y) < 3L || any(!is.finite(Y))) {
    stop("Y must be a finite numeric matrix with at least four rows and three columns.")
  }
  if (length(m) != 1L || !is.finite(m) || m != as.integer(m) || m < 1L || m >= min(dim(Y))) {
    stop("m must be a positive integer smaller than both n and p.")
  }
  invisible(TRUE)
}

poet_relative_frobenius <- function(estimate, target) {
  p <- nrow(target)
  e <- eigen((target + t(target)) / 2, symmetric = TRUE)
  if (any(e$values <= 0)) stop("Relative Frobenius target must be positive definite.")
  inv_sqrt <- tcrossprod(sweep(e$vectors, 2, 1 / sqrt(e$values), "*"), e$vectors)
  norm(inv_sqrt %*% (estimate - target) %*% inv_sqrt, "F") / sqrt(p)
}

poet_prepare_sim_checkpoints <- function(checkpoint_dir) {
  if (is.null(checkpoint_dir)) return(invisible(NULL))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(checkpoint_dir)) stop("Could not create checkpoint directory.")
  probe <- tempfile("write_check_", tmpdir = checkpoint_dir)
  saveRDS(list(write_check = TRUE), probe)
  unlink(probe)
  invisible(checkpoint_dir)
}

poet_sim_checkpoint <- function(checkpoint_dir, kind, mi, iteration, results, metadata, seed_before) {
  if (is.null(checkpoint_dir)) return(invisible(NULL))
  file <- file.path(checkpoint_dir, sprintf("%s_mi%d_rep%04d.rds", kind, mi, iteration))
  if (file.exists(file)) stop("Checkpoint already exists; use a new output directory: ", file)
  saveRDS(list(iteration = iteration, results = results, metadata = metadata,
               seed_before = seed_before, seed_after = get(".Random.seed", envir = .GlobalEnv),
               saved_at = Sys.time()), file)
  invisible(file)
}

poet_summarize_simulations <- function(simulation_results, methods) {
  means <- deviations <- setNames(vector("list", length(methods)), methods)
  for (method in methods) {
    metric_names <- names(simulation_results[[1L]][[method]])
    vals <- vapply(simulation_results, function(one) unlist(one[[method]], use.names = FALSE),
                   numeric(length(metric_names)))
    dim(vals) <- c(length(metric_names), length(simulation_results))
    means[[method]] <- as.list(setNames(rowMeans(vals), paste0("mean_", metric_names)))
    deviations[[method]] <- as.list(setNames(apply(vals, 1, stats::sd), paste0("sd_", metric_names)))
  }
  list(mean_results = means, sd_results = deviations)
}

poet_population_targets <- function(Sigma0, Sigma_u0, m, rel.tol = 1e-8) {
  p <- nrow(Sigma0)
  stopifnot(is.matrix(Sigma0), is.matrix(Sigma_u0),
            identical(dim(Sigma0), dim(Sigma_u0)), p > 2L,
            length(m) == 1L, is.finite(m), m == as.integer(m), m >= 1L, m < p)
  eig <- eigen((Sigma0 + t(Sigma0)) / 2, symmetric = TRUE)
  lambda <- eig$values
  if (any(!is.finite(lambda)) || min(lambda) <= 0) stop("Sigma0 must be positive definite.")
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
       Sigma0 = Sigma0, Sigma_u0 = Sigma_u0)
}

eigen_projection <- function(A, epsilon = 1e-5) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || any(!is.finite(A)) ||
      length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("PSD projection requires a finite square matrix and positive epsilon.")
  }
  eig <- eigen((A + t(A)) / 2, symmetric = TRUE)
  ans <- tcrossprod(sweep(eig$vectors, 2, sqrt(pmax(eig$values, epsilon)), "*"))
  (ans + t(ans)) / 2
}

adaptive_threshold <- function(S, tau) {
  diag_S <- diag(S)
  diag_pos <- pmax(diag_S, 1e-8)

  tau_mat <- outer(diag_pos, diag_pos, function(x, y) tau * sqrt(x * y))
  S_offdiag <- S
  diag(S_offdiag) <- 0

  S_thresh <- sign(S_offdiag) * pmax(abs(S_offdiag) - tau_mat, 0)
  diag(S_thresh) <- diag_S
  return(S_thresh)
}

generate_factor_data_elliptical <- function(B,
                                            n = 250,
                                            dist = c("normal", "t4", "t2.2", "mixnorm")) {
  dist <- match.arg(dist)

  d <- nrow(B)
  m <- ncol(B)

  idx <- matrix(1:d, d, d)
  Sigma_u <- 0.9 ^ abs(idx - t(idx))

  Sigma_y_raw <- B %*% t(B) + Sigma_u

  c0 <- d / sum(diag(Sigma_y_raw))

  Sigma_fu <- matrix(0, nrow = m + d, ncol = m + d)
  Sigma_fu[1:m, 1:m] <- c0 * diag(m)
  Sigma_fu[(m + 1):(m + d), (m + 1):(m + d)] <- c0 * Sigma_u

  Sigma0 <- c0 * Sigma_y_raw
  Sigma_u0 <- c0 * Sigma_u

  if (dist == "normal") {
    joint_sample <- mvtnorm::rmvnorm(n, mean = rep(0, m + d), sigma = Sigma_fu)

  } else if (dist == "t4") {
    nu <- 4
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_fu * (nu - 2) / nu, df = nu)

  } else if (dist == "t2.2") {
    nu <- 2.2
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_fu * (nu - 2) / nu, df = nu)

  } else if (dist == "mixnorm") {
    z <- rbinom(n, size = 1, prob = 0.2)
    joint_sample <- matrix(0, n, m + d)
    n1 <- sum(z == 0)
    n2 <- sum(z == 1)
    if (n1 > 0) {
      joint_sample[z == 0, ] <- mvtnorm::rmvnorm(n1, mean = rep(0, m + d), sigma = Sigma_fu)
    }
    if (n2 > 0) {
      joint_sample[z == 1, ] <- mvtnorm::rmvnorm(n2, mean = rep(0, m + d), sigma = 10 * Sigma_fu)
    }
    joint_sample <- joint_sample / sqrt(2.8)
  }

  f_mat <- joint_sample[, 1:m, drop = FALSE]
  u_mat <- joint_sample[, (m + 1):(m + d), drop = FALSE]

  Y <- f_mat %*% t(B) + u_mat

  return(list(
    data = Y,
    Sigma0 = Sigma0,
    Sigma_u0 = Sigma_u0
  ))
}

compute_2_part_nonsparse <- function(Y, m, method = "SPCA") {
  method <- match.arg(method)
  poet_validate_matrix_input(Y, m)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SPCA") {
    mu_hat <- ICSNP::spatial.median(Y)
    S_hat <- SpatialNP::SCov(Y, location = mu_hat) * p
    eig <- eigen(S_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, seq_len(m), drop = FALSE]
    Lambda_hat <- diag(eig$values[seq_len(m)], nrow = m, ncol = m)
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- S_hat - low_rank_part
    mu <- sqrt(log(p) / n)
    tau <- 0.8*(mu +sqrt(log(n)/n))
    R_hat <- adaptive_threshold(R_tilde, tau = tau)
    R_hat <- eigen_projection(R_hat,epsilon = 1e-3)
    Sigma_hat <- low_rank_part+R_hat
    Sigma_u_hat <- R_hat
    Sigma_hat  <- (Sigma_hat  + t(Sigma_hat )) / 2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    # SS estimates H0 rather than trace-normalized Sigma0.
    Sigma_hat0 <- Sigma_hat
    Sigma_u_hat0 <- Sigma_u_hat
    Lambda_hat0 <- Lambda_hat
  }
  return(list(
    Lambda_hat0 = Lambda_hat0,
    Gamma_hat = Gamma_hat,
    Sigma_u_hat0 = Sigma_u_hat0,
    Sigma_hat0 =Sigma_hat0
  ))
}

evaluate_poet_nonsparse <- function(Y, Sigma_true, Sigma_u_true,
                                     method = "SPCA",
                                     estimate_2_part = NULL, m = 3, population_targets = NULL) {
  method <- match.arg(method)
  p <- ncol(Y)
  eig <- eigen(Sigma_true, symmetric = TRUE)
  Gamma_true <- eig$vectors[, seq_len(m), drop = FALSE]
  Lambda_true <- diag(eig$values[seq_len(m)], nrow = m, ncol = m)
  if (is.null(estimate_2_part)) estimate_2_part <- compute_2_part_nonsparse(Y, m, method)
  Gamma <- estimate_2_part$Gamma_hat
  for (k in seq_len(m)) if (sum(Gamma[, k] * Gamma_true[, k]) < 0) Gamma[, k] <- -Gamma[, k]
  Sigma_hat <- estimate_2_part$Sigma_hat0
  Sigma_u_hat <- estimate_2_part$Sigma_u_hat0
  Lambda_hat <- estimate_2_part$Lambda_hat0
  gamma_error <- sqrt(p) * max(abs(Gamma - Gamma_true))
  if (method == "SPCA") {
    if (is.null(population_targets)) population_targets <- poet_population_targets(Sigma_true, Sigma_u_true, m)
    target <- population_targets
    return(list(
      Gamma_err = gamma_error,
      Lambda_sign_err = max(abs(Lambda_hat %*% solve(target$Lambda_sign) - diag(m))),
      H0_max_err = max(abs(Sigma_hat - target$H0)),
      H0_rF_err = poet_relative_frobenius(Sigma_hat, target$H0),
      H0u_spectral_err = norm(Sigma_u_hat - target$H0u, "2"),
      H0_inverse_spectral_err = norm(solve(Sigma_hat) - solve(target$H0), "2"),
      cross_target_Lambda0_discrepancy = max(abs(Lambda_hat %*% solve(Lambda_true) - diag(m))),
      cross_target_Sigma0_max_discrepancy = max(abs(Sigma_hat - Sigma_true)),
      cross_target_Sigma0_rF_discrepancy = poet_relative_frobenius(Sigma_hat, Sigma_true),
      cross_target_Sigma_u0_spectral_discrepancy = norm(Sigma_u_hat - Sigma_u_true, "2"),
      population_H0_vs_Sigma0_max = max(abs(target$H0 - Sigma_true)),
      population_H0_vs_Sigma0_rF = poet_relative_frobenius(target$H0, Sigma_true)))
  }
}

run_poet_nonsparse_simulation <- function(p = 500, sims = 100, n = 250, m = 3,
                                          mi = 1, checkpoint_dir = NULL, seed = NULL,
                                          methods = "SPCA") {
  if (!identical(methods, "SPCA")) stop("This script runs POET-SS only.")
  if (!(mi %in% 1:4) || length(mi) != 1L || sims < 1 || sims != as.integer(sims)) stop("Invalid scenario or repetitions.")
  if (m < 1 || m > 3 || m >= min(n,p)) stop("Simulation loading variances are specified for m=1,2,3 only.")
  if (!is.null(seed)) set.seed(seed)
  poet_prepare_sim_checkpoints(checkpoint_dir)
  dist <- c("normal", "t4", "t2.2", "mixnorm")[mi]
  metadata <- list(kind = "ss_hybrid", n = n, p = p, m = m, mi = mi, dist = dist,
                   sims = sims, seed = seed, methods = methods,
                   SS_target = "H0", student_t = "covariance standardized",
                   mixture = "(0.8*N(0,Sigma0)+0.2*N(0,10*Sigma0))/sqrt(2.8)",
                   thresholds = "existing finite-sample adaptive thresholds; not literal theorem b_n",
                   SS_pilot = "p * SpatialNP::SCov(Y, location = ICSNP::spatial.median(Y))",
                   session = utils::sessionInfo())
  simulation_results <- vector("list", sims)
  for (s in seq_len(sims)) {
    seed_before <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    sm <- c(1, 0.75^2, 0.5^2)
    B <- matrix(0, p, m)
    for (k in seq_len(m)) B[, k] <- rnorm(p, sd = sqrt(sm[k]))
    dat <- generate_factor_data_elliptical(B, n = n, dist = dist)
    population <- if ("SPCA" %in% methods) poet_population_targets(dat$Sigma0, dat$Sigma_u0, m) else NULL
    one <- setNames(vector("list", length(methods)), methods)
    for (method in methods) {
      elapsed <- unname(system.time({ estimate <- compute_2_part_nonsparse(dat$data, m, method) })[["elapsed"]])
      one[[method]] <- evaluate_poet_nonsparse(dat$data, dat$Sigma0, dat$Sigma_u0,
                      method, estimate, m = m, population_targets = population)
      one[[method]]$runtime_sec <- elapsed
    }
    simulation_results[[s]] <- one
    poet_sim_checkpoint(checkpoint_dir, "ss_hybrid", mi, s, one, metadata, seed_before)
    if (isTRUE(getOption("poet.progress", TRUE))) cat(sprintf("SS-H0 scenario %d: %d/%d completed\n", mi, s, sims))
  }
  summary <- poet_summarize_simulations(simulation_results, methods)
  c(list(simulation_results = simulation_results, metadata = metadata), summary)
}

poet_ss_csv <- function(table, file) {
  for (column in intersect(c("mean", "sd"), names(table))) table[[column]] <- sprintf("%.4f", table[[column]])
  utils::write.csv(table, file, row.names = FALSE)
}

poet_ss_latex <- function(table, settings, file) {
  metrics <- c("Lambda_sign_err", "Gamma_err", "H0_rF_err", "H0u_spectral_err", "H0_inverse_spectral_err")
  rows <- vapply(settings$scenarios, function(mi) {
    one <- table[table$scenario == mi, ]
    values <- one$value[match(metrics, one$metric)]
    if (anyNA(values)) stop("Missing SS-H0 table metrics.")
    paste0(c("I", "II", "III", "IV")[mi], " & ", paste(values, collapse = " & "), " \\\\")
  }, character(1))
  caption <- sprintf(paste0("Performance of POET-SS for the population hybrid target $\\H_0$. ",
                    "Entries are means (standard deviations) over %d replications, with $n=%d$, $d=%d$, and $m=%d$."),
                    settings$sims, settings$n, settings$p, settings$m)
  lines <- c("% Requires xcolor and the paper's matrix macros.", "\\begin{table}[htbp]", "\\centering", "\\color{red}",
             paste0("\\caption{\\textcolor{red}{", caption, "}}"), "\\label{tab:ss_hybrid}",
             "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{lccccc}", "\\hline",
             "Scenario & $\\|\\hat\\L_m^S(\\L_m^S)^{-1}-\\I_m\\|_{\\max}$ & $\\sqrt d\\|\\hat\\G_m-\\G_m\\|_{\\max}$ & $\\|\\hat\\H_0^\\tau-\\H_0\\|_{\\H_0}$ & $\\|\\hat\\H_{0u}^\\tau-c_d\\bms_u\\|_2$ & $\\|(\\hat\\H_0^\\tau)^{-1}-\\P_0\\|_2$ \\\\",
             "\\hline", rows, "\\hline", "\\end{tabular}", "\\par\\medskip",
             "{\\footnotesize Scenarios I--IV are normal, $t_4$, $t_{2.2}$, and contaminated normal, respectively.}",
             "\\end{table}")
  writeLines(lines, file, useBytes = TRUE)
}

if (isTRUE(getOption("poet.autorun", TRUE))) {
  # Parameters: set mi_vec <- 4 for scenario IV only.
  n <- 250L
  p <- 500L
  m <- 3L
  sims <- 100L
  mi_vec <- 1:4                 # 1=normal, 2=t4, 3=t2.2, 4=contaminated normal
  seed <- 988L

  .poet_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(.poet_file)) {
    .poet_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(.poet_arg) && basename(sub("^--file=", "", .poet_arg[1])) == "POET_SS_H0_vH.R")
      .poet_file <- sub("^--file=", "", .poet_arg[1])
  }
  .poet_dir <- if (length(.poet_file) == 1L && !is.na(.poet_file) && file.exists(.poet_file))
    dirname(normalizePath(.poet_file, winslash = "/")) else getwd()
  output_dir <- file.path(.poet_dir, "results", "ss_hybrid")

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
    run_dir <- file.path(output_dir, sprintf("ss_hybrid_run_%03d", .poet_index))
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
    result <- run_poet_nonsparse_simulation(p = p, n = n, m = m, sims = sims, mi = mi, seed = seed,
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
    .poet_path <- file.path(run_dir, paste0("ss_hybrid_mi", mi))
    saveRDS(list(result = result, table = table), paste0(.poet_path, ".rds"))
    tryCatch(poet_ss_csv(table, paste0(.poet_path, ".csv")),
             error = function(e) warning("RDS saved; CSV export failed: ", conditionMessage(e)))
  }
  summary_table <- do.call(rbind, tables)
  rownames(summary_table) <- NULL
  .poet_path <- file.path(run_dir, "ss_hybrid_all_scenarios")
  saveRDS(list(result = results, table = summary_table), paste0(.poet_path, ".rds"))
  tryCatch(poet_ss_csv(summary_table, paste0(.poet_path, ".csv")),
           error = function(e) warning("RDS saved; CSV export failed: ", conditionMessage(e)))
  appendix_metrics <- c("Lambda_sign_err", "Gamma_err", "H0_rF_err", "H0u_spectral_err", "H0_inverse_spectral_err")
  appendix_table <- summary_table[summary_table$metric %in% appendix_metrics, ]
  appendix_table <- appendix_table[order(match(appendix_table$scenario, mi_vec),
                                          match(appendix_table$metric, appendix_metrics)), ]
  rownames(appendix_table) <- NULL
  poet_ss_csv(appendix_table, file.path(run_dir, "POET_SS_H0_table.csv"))
  poet_ss_latex(appendix_table, settings, file.path(run_dir, "POET_SS_H0_table.tex"))
  cat("SAVED: ", run_dir, "\n", sep = "")
}
