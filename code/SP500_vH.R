# SP500 rolling portfolios. TylerSPCA = POET-TME; SPCA = POET-SS.
# Set options(poet.autorun = FALSE) before source() to load functions only.
poet_source_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                             error = function(e) "")
if (length(poet_source_file) != 1L || !nzchar(poet_source_file)) {
  entry <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  entry <- if (length(entry)) sub("^--file=", "", entry[[1L]]) else ""
  poet_source_file <- if (nzchar(entry) && identical(basename(entry), "SP500_vH.R") && file.exists(entry))
    normalizePath(entry, winslash = "/", mustWork = TRUE) else ""
}

poet_require_packages <- function() {
  packages <- c("SpatialNP", "ICSNP", "pcaPP", "readr", "lubridate", "dplyr")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install these R packages before running: ", paste(absent, collapse = ", "))
}

poet_md5_object <- function(object) {
  path <- tempfile("poet_fingerprint_", fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(object, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

poet_atomic_rds <- function(object, destination) {
  temp <- tempfile("checkpoint_", tmpdir = dirname(destination), fileext = ".tmp")
  on.exit(if (file.exists(temp)) unlink(temp))
  saveRDS(object, temp, compress = FALSE)
  if (file.exists(destination)) stop("Refusing to replace an existing checkpoint: ", destination)
  if (!file.rename(temp, destination)) stop("Cannot move checkpoint into place: ", destination)
  invisible(destination)
}

poet_new_run_directory <- function(output_dir, prefix) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  index <- 1L
  repeat {
    run_dir <- file.path(output_dir, sprintf("%s_run_%03d", prefix, index))
    if (!file.exists(run_dir)) {
      if (dir.create(run_dir, showWarnings = FALSE))
        return(normalizePath(run_dir, winslash = "/", mustWork = TRUE))
      if (!file.exists(run_dir)) stop("Cannot create output directory: ", run_dir)
    }
    index <- index + 1L
  }
}

poet_prepare_output <- function(output_dir, config, resume_dir = NULL) {
  if (!is.null(resume_dir)) {
    run_dir <- normalizePath(resume_dir, winslash = "/", mustWork = TRUE)
    info <- readRDS(file.path(run_dir, "run_info.rds"))
    if (!identical(info$config[names(config)], config) ||
        !identical(info$config$gr_result_md5, config$gr_result_md5)) {
      stop("Resume settings or input data differ. Start a new run.")
    }
  } else {
    run_dir <- poet_new_run_directory(output_dir, "sp500")
    dir.create(file.path(run_dir, "months"), showWarnings = FALSE)
    info <- list(config = config, started_at = Sys.time(), run_dir = run_dir, session_info = sessionInfo())
    poet_atomic_rds(info, file.path(run_dir, "run_info.rds"))
  }
  probe <- tempfile("write_probe_", tmpdir = file.path(run_dir, "months"))
  saveRDS(TRUE, probe)
  unlink(probe)
  cat("SP500_OUTPUT_DIR: ", run_dir, "\n", sep = "")
  flush.console()
  run_dir
}

poet_checkpoint_month <- function(month, monthly_returns, monthly_mhat, run_dir, config_hash) {
  assign("SP500_progress", list(run_dir = run_dir, completed_month = month,
    monthly_returns = monthly_returns, monthly_mhat = monthly_mhat), envir = .GlobalEnv)
  keep <- startsWith(names(monthly_returns), paste0(month, "_"))
  payload <- list(month = month, monthly_returns = monthly_returns[keep],
    monthly_mhat = monthly_mhat[[month]], config_hash = config_hash, completed_at = Sys.time())
  destination <- file.path(run_dir, "months", paste0("month_", month, ".rds"))
  tryCatch(poet_atomic_rds(payload, destination), error = function(e) {
    stop("Monthly save failed. Results remain in SP500_progress; do not close R. ", conditionMessage(e))
  })
}

poet_finish_output <- function(final_result, run_dir) {
  info <- list(run_dir = run_dir, completed_at = Sys.time(), complete = TRUE)
  assign("SP500_latest_result", final_result, envir = .GlobalEnv)
  assign("SP500_latest_result_info", info, envir = .GlobalEnv)
  destination <- file.path(run_dir, "SP500_results.rds")
  if (file.exists(destination)) {
    saved <- tryCatch(readRDS(destination), error = function(e) {
      stop("Existing final RDS cannot be read. Complete result remains in SP500_latest_result; the existing path will not be overwritten. ", conditionMessage(e))
    })
    if (!is.list(saved) || !is.list(saved$result) ||
        !identical(saved$result$config[names(final_result$config)], final_result$config) ||
        !identical(saved$result$all_returns, final_result$all_returns)) {
      stop("Existing final RDS does not match this completed run. Complete result remains in SP500_latest_result; the existing file will not be overwritten.")
    }
  } else {
    tryCatch(poet_atomic_rds(list(result = final_result, info = info), destination), error = function(e) {
      stop("Complete result is retained in SP500_latest_result, but RDS save failed. Do not close R. ", conditionMessage(e))
    })
  }
  cat("FULL_RESULT_SAVED: ", destination, "\n", sep = "")
  tables <- c("yearly_sd", "yearly_pval", "overall_sd", "overall_pval", "monthly_mhat")
  for (name in tables) {
    tryCatch(readr::write_csv(final_result[[name]], file.path(run_dir, paste0(name, ".csv"))), error = function(e) {
      warning("CSV export failed for ", name, "; complete RDS is already saved. ", conditionMessage(e), call. = FALSE)
    })
  }
  final_result
}

poet_read_returns <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  x <- readr::read_csv(path, col_types = readr::cols(date = readr::col_date(), .default = readr::col_double()),
                       show_col_types = FALSE, progress = FALSE, name_repair = "check_unique")
  if (nrow(readr::problems(x))) stop("CSV contains values that could not be parsed: ", path)
  if (!"date" %in% names(x) || anyNA(x$date) || anyDuplicated(x$date)) stop("Missing, invalid or duplicated dates in ", path)
  stocks <- setdiff(names(x), "date")
  if (!length(stocks) || any(vapply(x[stocks], function(y) any(is.infinite(y)), logical(1)))) stop("Invalid stock columns in ", path)
  x <- as.data.frame(x[order(x$date), ])
  x$DATE <- x$date
  x$year <- lubridate::year(x$date)
  x$month <- lubridate::month(x$date)
  x
}

poet_month_seed <- function(seed, month, method) {
  h <- as.double(seed) %% 2147483646
  for (value in utf8ToInt(paste(month, method, sep = ":"))) h <- (h * 131 + value) %% 2147483646
  as.integer(h + 1)
}

poet_with_seed <- function(seed, expression) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  force(expression)
}

poet_spd_shape <- function(A, stage) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || any(!is.finite(A))) stop("Nonfinite matrix at ", stage)
  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE)
  floor <- 1e-8 * max(1, max(eig$values))
  original_min <- min(eig$values)
  corrected <- original_min < floor
  if (corrected) A <- tcrossprod(sweep(eig$vectors, 2, sqrt(pmax(eig$values, floor)), "*"))
  tr <- sum(diag(A))
  if (!is.finite(tr) || tr <= 0) stop("Nonpositive trace at ", stage)
  A <- A * nrow(A) / tr
  A
}

poet_shape_inverse <- function(A, stage) {
  A <- poet_spd_shape(A, stage)
  list(Sigma_hat = A, Sigma_hat_inv = chol2inv(chol(A)))
}

poet_ipsn_hac_test <- function(ri, rn) {
  ri <- as.numeric(ri); rn <- as.numeric(rn)
  T_len <- length(ri)
  if (length(rn) != T_len || T_len <= 4L || any(!is.finite(c(ri, rn)))) stop("HAC needs more than four aligned finite return pairs.")
  mu_i <- mean(ri); mu_n <- mean(rn)
  gamma_i <- mean(ri^2); gamma_n <- mean(rn^2)
  denom_i <- mean((ri - mu_i)^2); denom_n <- mean((rn - mu_n)^2)
  if (min(denom_i, denom_n) <= 0) stop("HAC variance comparison is undefined for constant returns.")
  Delta_hat <- log(denom_i) - log(denom_n)
  y_t <- cbind(ri - mu_i, rn - mu_n, ri^2 - gamma_i, rn^2 - gamma_n)
  fit <- apply(y_t, 2, function(y) {
    x <- y[-length(y)]; z <- y[-1L]
    ss <- sum(x^2)
    rho <- if (ss > 0) sum(x * z) / ss else 0
    rho <- max(-0.999, min(0.999, rho))
    c(rho = rho, innovation_variance = mean((z - rho * x)^2))
  })
  rho_vec <- fit[1L, ]; sigma_sq_vec <- fit[2L, ]
  denominator <- sum(sigma_sq_vec^2 / (1 - rho_vec)^4)
  numerator <- sum(4 * rho_vec^2 * sigma_sq_vec^2 / (1 - rho_vec)^8)
  alpha2_hat <- if (denominator > 0) numerator / denominator else 0
  bw <- 1.3221 * (alpha2_hat * T_len)^(1/5)
  kern_QS <- function(x) {
    if (abs(x) < 1e-8) return(1)
    z <- 6 * pi * x / 5
    25 / (12 * pi^2 * x^2) * (sin(z) / z - cos(z))
  }
  Psi_hat <- crossprod(y_t) / T_len
  if (is.finite(bw) && bw > 0) for (j in seq_len(T_len - 1L)) {
    weight <- kern_QS(j / bw)
    if (abs(weight) < 1e-10) next
    cov_j <- crossprod(y_t[(j + 1L):T_len, , drop = FALSE], y_t[seq_len(T_len - j), , drop = FALSE]) / T_len
    Psi_hat <- Psi_hat + weight * (cov_j + t(cov_j))
  }
  Psi_hat <- T_len * Psi_hat / (T_len - 4)
  gradient <- c(-2 * mu_i / denom_i, 2 * mu_n / denom_n, 1 / denom_i, -1 / denom_n)
  var_delta <- as.numeric(crossprod(gradient, Psi_hat %*% gradient)) / T_len
  if (!is.finite(var_delta) || var_delta < -1e-10) stop("Invalid HAC variance estimate.")
  sd_delta <- sqrt(max(0, var_delta))
  statistic <- if (sd_delta > 1e-14) Delta_hat / sd_delta else if (abs(Delta_hat) <= 1e-12) 0 else sign(Delta_hat) * Inf
  list(Delta_hat = Delta_hat, test_statistic = statistic, p_value = pnorm(statistic),
       sd_Delta_hat = sd_delta, bw = bw, rho = rho_vec, sigma2 = sigma_sq_vec,
       alpha2_hat = alpha2_hat)
}

tune_factor_number <- function(X, kmax, method = c("SAMPLE", "FLW", "SPCA", "TylerSPCA", "IPSN", "RegTME"),
                                factor_rule = c("GR", "ER", "fixed"), fixed_m = NULL) {
  method <- match.arg(method); factor_rule <- match.arg(factor_rule)
  upper_rank <- min(nrow(X), ncol(X)) - 1L
  if (factor_rule == "fixed") {
    if (length(fixed_m) != 1L || !is.finite(fixed_m) || fixed_m != as.integer(fixed_m) || fixed_m < 1L || fixed_m > upper_rank) stop("fixed_m must lie between 1 and min(n,p)-1.")
    selected <- as.integer(fixed_m)
    ratios <- numeric()
  } else {
    kmax <- min(as.integer(kmax), upper_rank - 1L)
    if (!is.finite(kmax) || kmax < 1L) stop("Too few observations/assets for factor-ratio selection.")
    pilot <- Sigma_hat_estimate(X, method = method)
    eigenvalues <- eigen(pilot, symmetric = TRUE, only.values = TRUE)$values[seq_len(upper_rank)]
    numerical_floor <- max(eigenvalues) * .Machine$double.eps * max(dim(X))
    if (!is.finite(numerical_floor) || numerical_floor <= 0) stop("Factor pilot has no positive spectrum.")
    eigenvalues <- pmax(eigenvalues, numerical_floor)
    index <- seq_len(kmax)
    if (factor_rule == "ER") ratios <- eigenvalues[index] / eigenvalues[index + 1L]
    else {
      tails <- rev(cumsum(rev(eigenvalues)))
      ratios <- log1p(eigenvalues[index] / tails[index]) / log1p(eigenvalues[index + 1L] / tails[index + 1L])
    }
    if (any(!is.finite(ratios))) stop("Invalid factor selection ratios.")
    selected <- which.max(ratios)
  }
  poet_log("factor count", paste0("mhat=", selected))
  list(best_m_mktcr = selected, method = method, ratio_seq = ratios,
       factor_rule = factor_rule)
}

compute_weights_nonparse <- function(Y, kmax = 10L,
                                      method = c("SAMPLE", "SPCA", "FLW", "TylerSPCA", "IPSN", "RegTME", "EW"),
                                      factor_rule = c("GR", "ER", "fixed"), fixed_m = NULL) {
  method <- match.arg(method); factor_rule <- match.arg(factor_rule)
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)) || nrow(Y) < 5L || ncol(Y) < 2L) stop("Portfolio input must be a finite numeric matrix, n>=5, p>=2.")
  p <- ncol(Y)
  if (method == "EW") {
    weights <- matrix(1 / p, p, 1L)
    selected <- 0L
  } else {
    tuning <- poet_step("factor selection", tune_factor_number(Y, kmax, method, factor_rule, fixed_m))
    selected <- tuning$best_m_mktcr
    estimate <- poet_step("matrix estimation", compute_2_part_nonsparse(Y, selected, method))
    ones <- rep(1, p)
    numerator <- as.vector(estimate$Sigma_hat_inv %*% ones)
    denominator <- sum(numerator)
    if (!is.finite(denominator) || denominator <= 0) stop("Nonpositive MVP normalization denominator.")
    weights <- matrix(numerator / denominator, p, 1L)
  }
  rownames(weights) <- colnames(Y); colnames(weights) <- "weight"
  if (any(!is.finite(weights)) || abs(sum(weights) - 1) > 1e-8) stop("Invalid portfolio weights.")
  attr(weights, "factor_number") <- selected
  weights
}

poet_preflight_months <- function(train_data_raw, test_data_raw, start_date, end_date, window_months) {
  first <- as.Date(start_date); last <- as.Date(end_date)
  if (is.na(first) || is.na(last) || first > last || format(first, "%d") != "01") stop("start_date must be a month start and must not exceed end_date.")
  last_month <- as.Date(format(last, "%Y-%m-01"))
  if (last != lubridate::`%m+%`(last_month, lubridate::period(month = 1)) - 1L) stop("end_date must be a calendar month end.")
  starts <- seq(first, last_month, by = "month")
  stocks <- setdiff(names(train_data_raw), c("date", "DATE", "year", "month"))
  if (!all(stocks %in% names(test_data_raw))) stop("Test data are missing training asset columns.")
  train_matrix <- as.matrix(train_data_raw[stocks]); test_matrix <- as.matrix(test_data_raw[stocks])
  plans <- lapply(starts, function(day) {
    day <- as.Date(day, origin = "1970-01-01")
    train_start <- lubridate::`%m-%`(day, lubridate::period(month = window_months))
    next_month <- lubridate::`%m+%`(day, lubridate::period(month = 1))
    tr <- which(train_data_raw$DATE >= train_start & train_data_raw$DATE < day)
    te <- which(test_data_raw$DATE >= day & test_data_raw$DATE < next_month)
    month_label <- format(day, "%Y-%m")
    if (length(tr) < 5L || min(train_data_raw$DATE[tr]) > train_start + 7L ||
        max(train_data_raw$DATE[tr]) < day - 7L) stop("Incomplete training window for ", month_label)
    known_test_dates <- train_data_raw$DATE[train_data_raw$DATE >= day & train_data_raw$DATE < next_month]
    test_dates <- test_data_raw$DATE[te]
    missing_dates <- unique(known_test_dates[!known_test_dates %in% test_dates])
    if (length(missing_dates)) {
      shown <- paste(format(head(missing_dates, 10L), "%Y-%m-%d"), collapse = ", ")
      suffix <- if (length(missing_dates) > 10L) paste0(" (and ", length(missing_dates) - 10L, " more)") else ""
      stop("Incomplete test window for ", month_label,
           ": dates present in the training CSV are missing from the test CSV: ", shown, suffix)
    }
    if (!length(te)) stop("Incomplete test window for ", month_label, ": no test dates.")
    if (min(test_dates) > day + 7L || max(test_dates) < next_month - 7L)
      stop("Incomplete test window for ", month_label, ": observed dates ",
           format(min(test_dates), "%Y-%m-%d"), " to ", format(max(test_dates), "%Y-%m-%d"),
           " do not cover both month boundaries within the seven-day allowance.")
    active <- which(colSums(!is.finite(train_matrix[tr, , drop = FALSE])) == 0L)
    if (length(active) < 2L) stop("Fewer than two eligible assets for ", format(day, "%Y-%m"))
    x <- test_matrix[te, active, drop = FALSE]
    list(month = format(day, "%Y-%m"), train_rows = tr, test_rows = te, stocks = stocks[active],
         train_start = train_start, train_end = day - 1L,
         missing_cells = sum(is.na(x)))
  })
  plans
}

poet_summarise_returns <- function(all_returns, annualization_days = NULL) {
  all_returns <- all_returns[order(all_returns$date, all_returns$method), ]
  all_returns$year <- lubridate::year(all_returns$date)
  yearly <- dplyr::summarise(dplyr::group_by(all_returns, year, method),
                             annualized_sd = stats::sd(portfolio_return) * sqrt(if (is.null(annualization_days)) dplyr::n() else annualization_days),
                             .groups = "drop")
  days_per_year <- table(lubridate::year(unique(all_returns$date)))
  annual_days <- if (is.null(annualization_days)) round(mean(as.numeric(days_per_year))) else annualization_days
  overall <- dplyr::summarise(dplyr::group_by(all_returns, method),
                              annualized_sd = stats::sd(portfolio_return) * sqrt(annual_days), .groups = "drop")
  compare <- function(frame) {
    others <- setdiff(unique(frame$method), "TylerSPCA")
    if (!"TylerSPCA" %in% frame$method || !length(others)) return(data.frame(method = character(), pval = numeric()))
    base <- frame[frame$method == "TylerSPCA", c("date", "portfolio_return")]
    dplyr::bind_rows(lapply(others, function(method) {
      other <- frame[frame$method == method, c("date", "portfolio_return")]
      pair <- merge(base, other, by = "date", suffixes = c("_tme", "_other"), sort = TRUE)
      if (nrow(pair) != nrow(base) || nrow(pair) != nrow(other)) stop("Methods do not have matching return dates.")
      pval <- if (nrow(pair) > 4L) poet_ipsn_hac_test(pair$portfolio_return_tme, pair$portfolio_return_other)$p_value else NA_real_
      data.frame(method = method, pval = pval)
    }))
  }
  yearly_pval <- dplyr::bind_rows(lapply(sort(unique(all_returns$year)), function(y) {
    tab <- compare(all_returns[all_returns$year == y, ])
    tab$year <- rep(y, nrow(tab)); tab[, c("year", "method", "pval")]
  }))
  list(yearly_sd = yearly, yearly_pval = yearly_pval, overall_sd = overall,
       overall_pval = compare(all_returns), all_returns = all_returns,
       overall_annualization_days = annual_days)
}

compute_yearly_risk <- function(train_data_raw, test_data_raw, methods, config, output_dir,
                                resume_dir = NULL, stop_after_months = NULL,
                                gr_result = NULL, pilot_cache_dir = NULL) {
  old_cache_dir <- poet_state$pilot_cache_dir
  on.exit(poet_state$pilot_cache_dir <- old_cache_dir)
  poet_state$pilot_cache_dir <- pilot_cache_dir
  if (!is.null(pilot_cache_dir)) {
    dir.create(pilot_cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(pilot_cache_dir)) stop("Cannot create SSCov cache directory.")
  }
  plans <- poet_preflight_months(train_data_raw, test_data_raw, config$start_date, config$end_date, config$window_months)
  if (config$missing_policy == "error" && any(vapply(plans, function(x) x$missing_cells > 0L, logical(1)))) {
    stop("Preflight found missing out-of-sample returns.")
  }
  run_dir <- poet_prepare_output(output_dir, config, resume_dir)
  config_hash <- poet_md5_object(readRDS(file.path(run_dir, "run_info.rds"))$config)
  monthly_returns <- monthly_mhat <- list()
  completed <- character()
  for (plan in plans) {
    path <- file.path(run_dir, "months", paste0("month_", plan$month, ".rds"))
    if (!file.exists(path)) next
    payload <- readRDS(path)
    if (!identical(payload$config_hash, config_hash) || !identical(payload$month, plan$month) ||
        !setequal(names(payload$monthly_returns), paste(plan$month, methods, sep = "_"))) stop("Checkpoint identity does not match this run: ", path)
    monthly_returns <- c(monthly_returns, payload$monthly_returns)
    monthly_mhat[[plan$month]] <- payload$monthly_mhat
    completed <- c(completed, plan$month)
  }
  if (length(completed) && !identical(completed, vapply(plans[seq_along(completed)], `[[`, character(1), "month"))) stop("Checkpoint months are not a contiguous completed prefix.")
  new_completed <- 0L
  for (i in seq_along(plans)) {
    plan <- plans[[i]]; current_ym <- plan$month
    if (current_ym %in% completed) next
    poet_begin_month(current_ym)
    month_started <- proc.time()[["elapsed"]]
    Y <- as.matrix(train_data_raw[plan$train_rows, plan$stocks, drop = FALSE])
    X_test <- as.matrix(test_data_raw[plan$test_rows, plan$stocks, drop = FALSE])
    X_test[is.na(X_test)] <- 0  # Preserve the original zero-excess-return convention.
    dates <- test_data_raw$DATE[plan$test_rows]
    count_rows <- list()
    reused_this_month <- 0L
    for (method in methods) {
      poet_state$method <- method
      method_seed <- poet_month_seed(config$seed, current_ym, method)
      key <- paste(current_ym, method, sep = "_")
      from_gr <- NULL
      if (!is.null(gr_result)) {
        candidate <- gr_result$counts[gr_result$counts$month == current_ym & gr_result$counts$method == method, ]
        if (nrow(candidate) != 1L || is.na(candidate$selected_m) ||
            candidate$seed != method_seed || candidate$n_train != nrow(Y) || candidate$n_assets != ncol(Y))
          stop("GR monthly data or seed mismatch: ", key)
        if (method == "EW" || candidate$selected_m == config$fixed_m) from_gr <- candidate
      }
      if (!is.null(from_gr)) {
        copied <- gr_result$returns[[key]]
        if (is.null(copied) || !identical(as.Date(copied$date), as.Date(dates))) stop("GR return dates differ: ", key)
        values <- copied$portfolio_return
        selected <- from_gr$selected_m
        reused_this_month <- reused_this_month + 1L
      } else {
        weights <- poet_with_seed(method_seed,
          poet_step("weights TOTAL", compute_weights_nonparse(Y, config$kmax, method, config$factor_rule, config$fixed_m)))
        w <- as.numeric(weights[, 1L])
        values <- as.vector(X_test %*% w)
        selected <- attr(weights, "factor_number")
      }
      if (any(!is.finite(values))) stop("Nonfinite portfolio returns for ", current_ym, " ", method)
      monthly_returns[[key]] <- data.frame(date = dates, method = method, portfolio_return = values)
      count_rows[[method]] <- data.frame(month = current_ym, method = method, selected_m = selected,
                                           factor_rule = if (method == "EW") "none" else config$factor_rule,
                                           seed = method_seed,
                                           n_train = nrow(Y), n_assets = ncol(Y))
    }
    monthly_mhat[[current_ym]] <- dplyr::bind_rows(count_rows)
    poet_checkpoint_month(current_ym, monthly_returns, monthly_mhat, run_dir, config_hash)
    reuse_note <- if (is.null(gr_result)) "" else sprintf("; %d/%d methods reused from GR", reused_this_month, length(methods))
    cat(sprintf("[%d/%d] %s completed and saved (%.1f s%s)\n", i, length(plans), current_ym, proc.time()[["elapsed"]] - month_started, reuse_note))
    flush.console()
    new_completed <- new_completed + 1L
    if (!is.null(stop_after_months) && new_completed >= stop_after_months && i < length(plans)) {
      return(invisible(list(complete = FALSE, run_dir = run_dir, completed_month = current_ym)))
    }
  }
  final_result <- poet_summarise_returns(dplyr::bind_rows(monthly_returns), config$annualization_days)
  final_result$monthly_mhat <- dplyr::bind_rows(monthly_mhat)
  final_result$config <- config
  final_result$complete <- TRUE
  final_result$run_dir <- run_dir
  poet_finish_output(final_result, run_dir)
}

poet_reuse_file <- function(reuse_dir, label) {
  if (is.null(reuse_dir)) return(NULL)
  paths <- list.files(file.path(reuse_dir, label), pattern = "^SP500_results\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(paths) != 1L) stop("Expected one completed baseline result for ", label, " in ", reuse_dir)
  paths[[1L]]
}

poet_read_reuse_result <- function(path, config) {
  saved <- readRDS(path)
  old <- saved$result
  if (!isTRUE(saved$info$complete) || !isTRUE(old$complete)) stop("Reuse requires completed results: ", path)
  fields <- c("train_file", "test_file", "input_md5", "start_date", "end_date", "window_months",
              "factor_rule", "fixed_m", "kmax", "seed", "missing_policy", "annualization_days")
  if (!all(fields %in% names(old$config)) || !identical(old$config[fields], config[fields]))
    stop("Saved data or experiment settings differ from this rerun: ", path)
  if (!all(c("FLW", "IPSN") %in% config$methods) || !all(config$methods %in% old$config$methods))
    stop("A partial method update must rerun both FLW and IPSN and use methods present in the saved results.")
  list(result = old, path = normalizePath(path, winslash = "/", mustWork = TRUE))
}

poet_read_gr_result <- function(path, config) {
  saved <- readRDS(path)
  gr <- saved$result
  if (!isTRUE(saved$info$complete) || !isTRUE(gr$complete) ||
      !identical(gr$config$factor_rule, "GR") || config$factor_rule != "fixed")
    stop("Fixed-factor reuse requires completed GR results.")
  fields <- c("train_file", "test_file", "input_md5", "start_date", "end_date", "window_months",
              "kmax", "seed", "missing_policy", "annualization_days")
  if (any(config$methods %in% c("FLW", "IPSN"))) fields <- c(fields, "flw_ipsn_pilot")
  if (!all(fields %in% names(gr$config)) || !identical(gr$config[fields], config[fields]) ||
      !all(config$methods %in% gr$config$methods)) stop("GR settings, data or pilot differ from this run.")
  counts <- as.data.frame(gr$monthly_mhat[gr$monthly_mhat$method %in% config$methods, ])
  returns <- as.data.frame(gr$all_returns[gr$all_returns$method %in% config$methods, c("date", "method", "portfolio_return")])
  if (anyNA(counts[c("month", "method", "selected_m", "seed", "n_train", "n_assets")]) ||
      anyNA(returns$date) || any(!is.finite(returns$portfolio_return)) ||
      anyDuplicated(counts[c("month", "method")]) || anyDuplicated(returns[c("date", "method")]))
    stop("Invalid GR observations for fixed-factor reuse.")
  returns <- returns[order(returns$date, returns$method), ]
  list(counts = counts, returns = split(returns, paste(format(returns$date, "%Y-%m"), returns$method, sep = "_")))
}

poet_merge_reused_results <- function(fresh, saved) {
  old <- saved$result
  updated <- fresh$config$methods
  validate <- function(run) {
    returns <- run$all_returns; counts <- run$monthly_mhat
    if (!setequal(unique(returns$method), run$config$methods) ||
        !setequal(unique(counts$method), run$config$methods) ||
        anyNA(returns$date) || anyNA(counts$month) || any(!is.finite(returns$portfolio_return)) ||
        anyDuplicated(returns[c("date", "method")]) || anyDuplicated(counts[c("month", "method")]))
      stop("Invalid or duplicated observations in results to combine.")
    dates <- sort(unique(returns$date)); months <- sort(unique(counts$month))
    if (!identical(months, sort(unique(format(dates, "%Y-%m"))))) stop("Return dates and factor-selection months differ.")
    for (method in run$config$methods) {
      if (!identical(sort(returns$date[returns$method == method]), dates) ||
          !identical(sort(counts$month[counts$method == method]), months))
        stop("Methods have different date or month coverage.")
    }
    list(dates = dates, months = months)
  }
  if (!identical(validate(old), validate(fresh))) stop("Old and new results have different date coverage.")
  cols <- c("month", "method", "factor_rule", "seed", "n_train", "n_assets")
  comparable <- function(x) {
    x <- as.data.frame(x[x$method %in% updated, cols])
    x <- x[order(x$month, x$method), ]; rownames(x) <- NULL; x
  }
  if (!identical(comparable(old$monthly_mhat), comparable(fresh$monthly_mhat)))
    stop("Old and new monthly sample sizes, seeds or factor rules differ.")
  returns <- dplyr::bind_rows(old$all_returns[!old$all_returns$method %in% updated, ], fresh$all_returns)
  combined <- poet_summarise_returns(returns, fresh$config$annualization_days)
  combined$monthly_mhat <- dplyr::bind_rows(old$monthly_mhat[!old$monthly_mhat$method %in% updated, ], fresh$monthly_mhat)
  combined$monthly_mhat <- combined$monthly_mhat[order(combined$monthly_mhat$month, combined$monthly_mhat$method), ]
  combined$config <- fresh$config
  combined$config$methods <- old$config$methods
  combined$complete <- TRUE
  combined$run_dir <- file.path(fresh$run_dir, "combined")
  combined$reused_result_file <- saved$path
  combined$rerun_methods <- updated
  combined$rerun_dir <- fresh$run_dir
  dir.create(combined$run_dir, showWarnings = FALSE)
  poet_finish_output(combined, combined$run_dir)
}

run_sp500 <- function(train_file, test_file, output_dir,
                       methods = c("EW", "SAMPLE", "FLW", "TylerSPCA", "SPCA", "IPSN", "RegTME"),
                       start_date = "2005-01-01", end_date = "2023-12-31", window_months = 120L,
                       factor_rule = c("GR", "ER", "fixed"), fixed_m = NULL, kmax = 10L,
                       seed = 20261003L, missing_policy = c("cash_zero", "error"),
                       annualization_days = NULL, resume_dir = NULL, stop_after_months = NULL,
                       reuse_result_file = NULL, gr_result_file = NULL, pilot_cache_dir = NULL) {
  poet_require_packages()
  factor_rule <- match.arg(factor_rule); missing_policy <- match.arg(missing_policy)
  allowed <- c("EW", "SAMPLE", "FLW", "TylerSPCA", "SPCA", "IPSN", "RegTME")
  if (!length(methods) || anyDuplicated(methods) || !all(methods %in% allowed)) stop("Invalid or duplicated methods.")
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > 2147483646) stop("Invalid seed.")
  if (length(window_months) != 1L || !is.finite(window_months) || window_months < 1 || window_months != as.integer(window_months)) stop("window_months must be a positive integer.")
  if (length(kmax) != 1L || !is.finite(kmax) || kmax < 1 || kmax != as.integer(kmax)) stop("kmax must be a positive integer.")
  if (!is.null(annualization_days) && (length(annualization_days) != 1L || !is.finite(annualization_days) || annualization_days <= 0)) stop("annualization_days must be positive.")
  if (!is.null(stop_after_months) && (length(stop_after_months) != 1L || !is.finite(stop_after_months) || stop_after_months < 1)) stop("stop_after_months must be positive.")
  train_file <- normalizePath(train_file, winslash = "/", mustWork = TRUE)
  test_file <- normalizePath(test_file, winslash = "/", mustWork = TRUE)
  if (!is.null(resume_dir) && is.null(gr_result_file)) {
    gr_result_file <- readRDS(file.path(resume_dir, "run_info.rds"))$config$gr_result_file
  }
  if (factor_rule != "fixed") fixed_m <- NULL
  if (factor_rule == "fixed" && (length(fixed_m) != 1L || !is.finite(fixed_m) || fixed_m < 1 || fixed_m != as.integer(fixed_m))) stop("fixed_m must be a positive integer for fixed factors.")
  config <- list(train_file = train_file, test_file = test_file,
    input_md5 = unname(tools::md5sum(c(train_file, test_file))),
    methods = methods, start_date = as.character(as.Date(start_date)), end_date = as.character(as.Date(end_date)),
    window_months = as.integer(window_months), factor_rule = factor_rule,
    fixed_m = if (is.null(fixed_m)) NULL else as.integer(fixed_m), kmax = as.integer(kmax), seed = as.integer(seed),
    missing_policy = missing_policy, annualization_days = annualization_days,
    flw_ipsn_pilot = "SSCov")
  reused <- if (!is.null(reuse_result_file)) poet_read_reuse_result(reuse_result_file, config) else NULL
  gr_result <- if (!is.null(gr_result_file)) poet_read_gr_result(gr_result_file, config) else NULL
  if (!is.null(gr_result_file)) {
    config$gr_result_file <- normalizePath(gr_result_file, winslash = "/", mustWork = TRUE)
    config$gr_result_md5 <- unname(tools::md5sum(gr_result_file))
  }
  result <- compute_yearly_risk(poet_read_returns(train_file), poet_read_returns(test_file), methods,
                                config, output_dir, resume_dir, stop_after_months, gr_result, pilot_cache_dir)
  if (isTRUE(result$complete)) {
    if (!is.null(reused)) result <- poet_merge_reused_results(result, reused)
    print(result$overall_sd)
    if (isTRUE(getOption("poet.verbose", FALSE))) print(result$overall_pval)
  }
  invisible(result)
}

run_sp500_sensitivity <- function(train_file, test_file, output_dir,
                                   methods = c("EW", "SAMPLE", "FLW", "TylerSPCA", "SPCA", "IPSN", "RegTME"), factor_grid = c(1L, 2L, 3L),
                                   include_gr = TRUE, reuse_dir = NULL, reuse_gr = TRUE,
                                   gr_resume_dir = NULL, pilot_cache_dir = file.path(output_dir, "pilot_cache"), ...) {
  if (!length(factor_grid) || any(!is.finite(factor_grid)) || any(factor_grid < 1) || any(factor_grid != as.integer(factor_grid)) || anyDuplicated(factor_grid)) stop("factor_grid must contain distinct positive integers.")
  extra <- list(...)
  if (any(c("factor_rule", "fixed_m", "resume_dir", "stop_after_months", "reuse_result_file", "gr_result_file") %in% names(extra))) stop("Sensitivity owns factor_rule/fixed_m and requires completed configurations.")
  if (length(reuse_gr) != 1L || is.na(reuse_gr) || !is.logical(reuse_gr)) stop("reuse_gr must be TRUE or FALSE.")
  if (!include_gr && !is.null(gr_resume_dir)) stop("gr_resume_dir requires include_gr=TRUE.")
  labels <- c(if (include_gr) "GR", paste0("m", factor_grid))
  expected <- list()
  if (!is.null(reuse_dir)) {
    defaults <- formals(run_sp500)
    for (field in c("start_date", "end_date", "window_months", "kmax", "seed", "missing_policy", "annualization_days")) {
      value <- if (field %in% names(extra)) extra[[field]] else eval(defaults[[field]], environment(run_sp500))
      if (field %in% c("window_months", "kmax", "seed")) value <- as.integer(value)
      if (field %in% c("start_date", "end_date")) value <- as.character(as.Date(value))
      if (field == "missing_policy") value <- match.arg(value, c("cash_zero", "error"))
      expected[field] <- list(value)
    }
    expected$train_file <- normalizePath(train_file, winslash = "/", mustWork = TRUE)
    expected$test_file <- normalizePath(test_file, winslash = "/", mustWork = TRUE)
    expected$input_md5 <- unname(tools::md5sum(c(expected$train_file, expected$test_file)))
    expected$methods <- methods
  }
  reuse_files <- if (!is.null(reuse_dir)) setNames(lapply(labels, function(label) {
    path <- poet_reuse_file(reuse_dir, label); saved <- readRDS(path)
    if (!isTRUE(saved$info$complete) || !isTRUE(saved$result$complete)) stop("Incomplete baseline: ", path)
    expected_m <- if (label == "GR") NULL else as.integer(sub("m", "", label))
    if (!identical(saved$result$config$factor_rule, if (label == "GR") "GR" else "fixed") ||
        !identical(saved$result$config$fixed_m, expected_m)) stop("Baseline factor setting mismatch: ", path)
    expected$factor_rule <- if (label == "GR") "GR" else "fixed"
    expected["fixed_m"] <- list(expected_m)
    poet_read_reuse_result(path, expected)
    path
  }), labels) else NULL
  collection <- poet_new_run_directory(output_dir, "sensitivity")
  results <- setNames(vector("list", length(labels)), labels)
  gr_file <- NULL
  for (label in labels) {
    args <- c(list(train_file = train_file, test_file = test_file, output_dir = file.path(collection, label),
                    methods = methods, factor_rule = if (label == "GR") "GR" else "fixed",
                    fixed_m = if (label == "GR") NULL else as.integer(sub("m", "", label)),
                    reuse_result_file = if (is.null(reuse_files)) NULL else reuse_files[[label]],
                    gr_result_file = if (label != "GR" && reuse_gr) gr_file else NULL,
                    resume_dir = if (label == "GR") gr_resume_dir else NULL,
                    pilot_cache_dir = pilot_cache_dir), extra)
    run <- do.call(run_sp500, args)
    if (label == "GR") gr_file <- file.path(run$run_dir, "SP500_results.rds")
    tab <- as.data.frame(run$overall_sd); tab$factor_setting <- label; tab$run_dir <- run$run_dir
    results[[label]] <- tab
  }
  long <- dplyr::bind_rows(results)
  compact <- reshape(long[, c("factor_setting", "method", "annualized_sd")],
                       idvar = "factor_setting", timevar = "method", direction = "wide")
  names(compact) <- sub("^annualized_sd\\.", "", names(compact))
  summary_result <- list(compact = compact, runs = long, output_dir = collection)
  assign("SP500_latest_sensitivity", summary_result, envir = .GlobalEnv)
  poet_atomic_rds(summary_result, file.path(collection, "sensitivity_results.rds"))
  tables <- list(sensitivity_overall_long = long, sensitivity_overall_compact = compact)
  for (name in names(tables)) tryCatch(
    readr::write_csv(tables[[name]], file.path(collection, paste0(name, ".csv"))),
    error = function(e) warning("Sensitivity CSV export failed; summary RDS and individual runs are saved. ", conditionMessage(e), call. = FALSE))
  print(compact)
  invisible(summary_result)
}

poet_state <- new.env(parent = emptyenv())
poet_state$month <- ""
poet_state$method <- ""
poet_state$Y <- NULL
poet_state$cache <- new.env(parent = emptyenv())

poet_log <- function(stage, status, elapsed = NA_real_) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  suffix <- if (is.na(elapsed)) "" else sprintf(" | %.3f s", elapsed)
  if (isTRUE(getOption("poet.verbose", FALSE)) || identical(status, "ERROR")) {
    cat(sprintf("[%s] month=%s method=%s | %s | %s%s\n",
                stamp, poet_state$month, poet_state$method, stage, status, suffix))
    flush.console()
  }
  path <- Sys.getenv("POET_TIMING_FILE", unset = "")
  if (nzchar(path) && status != "START") {
    row <- data.frame(timestamp = stamp, month = poet_state$month,
                      method = poet_state$method, stage = stage,
                      status = status, elapsed_seconds = elapsed)
    write.table(row, path, sep = ",", row.names = FALSE,
                col.names = !file.exists(path), append = file.exists(path),
                qmethod = "double")
  }
}

poet_step <- function(stage, expression) {
  poet_log(stage, "START")
  started <- proc.time()[["elapsed"]]
  completed <- FALSE
  on.exit(poet_log(stage, if (completed) "DONE" else "ERROR",
                    proc.time()[["elapsed"]] - started))
  value <- force(expression)
  completed <- TRUE
  value
}

poet_begin_month <- function(month) {
  poet_state$month <- month
  poet_state$method <- ""
  poet_state$Y <- NULL
  poet_state$cache <- new.env(parent = emptyenv())
}

poet_cached_matrix <- function(key, X, calculate) {
  if (!identical(X, poet_state$Y)) {
    poet_state$Y <- X
    poet_state$cache <- new.env(parent = emptyenv())
  }
  if (exists(key, envir = poet_state$cache, inherits = FALSE)) {
    poet_log(key, "CACHE", 0)
    return(get(key, envir = poet_state$cache, inherits = FALSE))
  }
  value <- poet_step(key, calculate())
  assign(key, value, envir = poet_state$cache)
  value
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
  q <- vapply(seq_len(p), function(j) {
    counts <- as.double(rle(sort(Y[, j], method = "radix"))$lengths)
    1 - sum(counts * (counts - 1) / 2) / pairs
  }, numeric(1))
  active <- which(q > 0)
  result <- matrix(0, p, p)
  if (length(active) >= 2L) {
    tau_b <- pcaPP::cor.fk(Y[, active, drop = FALSE])
    result[active, active] <- tau_b * tcrossprod(sqrt(q[active]))  # Convert tau-b to tau-a.
  }
  diag(result) <- 1
  result
}

poet_sscov <- function(X) {
  poet_cached_matrix("SSCov", X, function() {
    cache_dir <- poet_state$pilot_cache_dir
    if (is.null(cache_dir)) return(SpatialNP::SSCov(X))
    if (length(cache_dir) != 1L || is.na(cache_dir) || !dir.exists(cache_dir)) {
      stop("SSCov cache directory does not exist.")
    }
    key <- poet_md5_object(list(X = X, SpatialNP = as.character(utils::packageVersion("SpatialNP"))))
    path <- file.path(cache_dir, paste0("SSCov_", key, ".rds"))
    valid_matrix <- function(value) {
      is.matrix(value) && is.numeric(value) &&
        identical(dim(value), c(ncol(X), ncol(X))) && all(is.finite(value))
    }
    if (file.exists(path)) {
      saved <- readRDS(path)
      if (!is.list(saved) || !identical(saved$key, key) || !valid_matrix(saved$matrix)) {
        stop("Invalid SSCov cache: ", path)
      }
      poet_log("SSCov disk", "CACHE", 0)
      return(saved$matrix)
    }
    value <- SpatialNP::SSCov(X)
    if (!valid_matrix(value)) stop("SSCov returned an invalid matrix.")
    poet_atomic_rds(list(key = key, matrix = value), path)
    value
  })
}

poet_median <- function(X) {
  poet_cached_matrix("spatial median", X, function() ICSNP::spatial.median(X))
}

poet_scov <- function(X) {
  poet_cached_matrix("SCov", X, function() SpatialNP::SCov(X, location = poet_median(X)))
}

poet_regtme <- function(Y) {
  poet_cached_matrix("RegTME", Y, function() regTME_svd(symmetrize_Y(Y), alpha = NULL))
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

select_gamma <- function(A, mu) {
  if (length(mu) != 1L || !is.finite(mu) || mu <= 0 || any(!is.finite(A))) {
    stop("The smoothing parameter must be positive and A finite.")
  }
  a <- sort(as.vector(abs(A)) / mu, decreasing = TRUE)
  if (!length(a)) stop("A must not be empty.")
  if (sum(a) <= 1) return(list(gamma = 0, u = length(a), threshold = 0))
  candidates <- (cumsum(a) - 1) / seq_along(a)
  u <- max(which(a > candidates))
  gamma <- max(0, candidates[u])
  list(gamma = gamma, u = u, threshold = mu * gamma)
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
  ans
}

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
estimate_Lambda <- function(Y, D_hat, m) {
  T_hat <- poet_step("FLW/Kendall tau-a", poet_kendall_tau_a(Y))
  R_hat <- sin(pi / 2 * T_hat)
  Sigma1_hat <- D_hat %*% R_hat %*% D_hat
  eig_vals <- poet_step("FLW/pilot eigenvalues",
                        eigen(Sigma1_hat, symmetric = TRUE, only.values = TRUE)$values)
  list(Lambda_hat = diag(eig_vals[seq_len(m)], nrow = m, ncol = m),
       Sigma1_hat = Sigma1_hat)
}

Sigma_hat_estimate <- function(Y, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
  }
  if (method == "FLW") {
    Sigma_hat <- poet_sscov(Y)
  }
  if (method == "SPCA") {
    Sigma_hat <- poet_scov(Y)*p
  }
  if (method == "TylerSPCA") {
    Sigma_hat <- poet_scov(Y)*p
  }
  if (method == "IPSN") {
    Sigma_hat <- poet_sscov(Y)
  }
  if (method == "RegTME") {
    Sigma_hat <- poet_regtme(Y)
  }
  return(Sigma_hat)
}

compute_2_part_nonsparse <- function(Y, m, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat <- Sigma_hat*c
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    tau <- 0.1*(sqrt(log(p) / n) +1/sqrt(p))
    Sigma_u_hat <- adaptive_threshold(Sigma_u_hat, tau = tau)
    Sigma_hat <- Sigma_u_hat + low_rank_part
  }
  
  if (method == "FLW") {
    D_hat <- poet_step("FLW/robust variances", robust_variance_diag(Y))
    estimate <- estimate_Lambda(Y, D_hat, m)
    Lambda_hat <- estimate$Lambda_hat
    Sigma_hat <- estimate$Sigma1_hat
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat <- Sigma_hat*c
    Lambda_hat <- Lambda_hat*c
    HK<-poet_sscov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE] 
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    Sigma_u_hat1 <- eigen_projection(Sigma_u_hat,epsilon = 1e-5)
    max_diff <- max(abs(Sigma_u_hat - Sigma_u_hat1))
    mu <- sqrt(log(p) / n)
    if (max_diff < 2*mu) {
      Sigma_u_hat <- Sigma_u_hat1
    } else {
      Sigma_u_hat <- poet_step("FLW/residual iterative correction",
                                accelerated_proximal_gradient(Sigma_u_hat, mu, max_iter = 300, tol = 5e-3))
    }
    tau <- 0.1*(sqrt(log(p) / n) +1/sqrt(p))
    Sigma_u_hat <- adaptive_threshold(Sigma_u_hat, tau = tau)
    Sigma_hat <- low_rank_part + Sigma_u_hat
    
  }
  if (method == "SPCA") {
    Sigma_hat <- poet_scov(Y)*p
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    mu <- sqrt(log(p) / n)
    tau <- 0.1*(mu +sqrt(log(n)/n))
    Sigma_u_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
  }
  if (method == "TylerSPCA") {
    Sigma_hat <- poet_scov(Y)*p
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    mu <- sqrt(log(p) / n)
    tau <- 0.1*(mu +sqrt(log(n)/n))
    Sigma_u_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
    pilot <- poet_shape_inverse(Sigma_hat, "TME/SS pilot")
    Sigma_hat_inv <- pilot$Sigma_hat_inv
    mu_hat <- poet_median(Y)
    Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")
    quad_forms <- rowSums((Y_centered %*% Sigma_hat_inv) * Y_centered) 
    if (any(!is.finite(quad_forms)) || any(quad_forms <= 0)) stop("TME has a nonpositive directional denominator.")
    Sigma_T <- (p / n) * crossprod(Y_centered / sqrt(quad_forms))
    c<-p/sum(diag(Sigma_T))
    Sigma_hat <- Sigma_T*c
    eig_Sigma_T <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma_T$vectors[, 1:m, drop = FALSE]
    if(length(eig_Sigma_T$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig_Sigma_T$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig_Sigma_T$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    Sigma_u_hat <- adaptive_threshold(Sigma_u_hat, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
  }
  if (method == "IPSN") {
    mu_hat <- estimate_mu_huber_cv(Y)$mu
    HK<-poet_sscov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    Y_centered <- sweep(Y, 2, mu_hat, "-")
    projected <- Y_centered - (Y_centered %*% Gamma_hat) %*% t(Gamma_hat)
    norms <- sqrt(rowSums(projected^2))
    if (any(!is.finite(norms)) || any(norms <= .Machine$double.eps * pmax(1, sqrt(rowSums(Y_centered^2))))) {
      stop("IPSN projected radius is zero or numerically degenerate.")
    }
    X_hat <- sqrt(p) * Y_centered / norms
    Sigma_0_raw <- crossprod(X_hat) / n
    eta_hat <- p / sum(diag(Sigma_0_raw))
    Sigma_0_hat <- Sigma_0_raw * eta_hat
    eig_Sigma0 <- eigen(Sigma_0_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma0$vectors[, 1:m, drop = FALSE] 
    if(length(eig_Sigma0$values[1:m]) == 1) {
      Lambda_hat0 <- matrix(eig_Sigma0$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eig_Sigma0$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    Sigma_u_hat0 <- Sigma_0_hat - low_rank_part
    tau <- 0.1*(sqrt(log(p) / n) + sqrt(log(p) /p))
    Sigma_u_hat0 <- adaptive_threshold(Sigma_u_hat0, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat0
    
  }
  if (method == "RegTME") {
    Y_sym <- symmetrize_Y(Y)
    Sigma_tme <- poet_regtme(Y)
    eig <- eigen(Sigma_tme, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    eigvals <- eig$values[1:m]
    if (m == 1) {
      Lambda_hat0 <- matrix(eigvals, nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eigvals)
    }
    n_sym <- nrow(Y_sym)
    Sigma_u_hat0 <- Sigma_tme - Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    tau <- 0.1*(sqrt(log(p) / n_sym)+sqrt(log(p) /p))
    Sigma_u_hat0 <- adaptive_threshold(Sigma_u_hat0, tau = tau)
    Sigma_hat <- Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)+Sigma_u_hat0
  }
  poet_shape_inverse(Sigma_hat, paste0(method, "/final"))
}
# 运行参数：一般只需修改下面的文件路径和日期。
if (isTRUE(getOption("poet.autorun", TRUE))) {
  train_file <- "D:/R/somework/POET/data/sp500_excess_return_wide.csv"
  test_file <- "D:/R/somework/POET/data/sp500_anymember_excess_return_wide.csv"
  project_dir <- if (nzchar(poet_source_file)) dirname(poet_source_file) else getwd()
  output_dir <- file.path(project_dir, "results", "sp500_sscov")
  start_date <- "2005-01-01"
  end_date <- "2023-12-31"
  methods <- c("EW", "SAMPLE", "FLW", "TylerSPCA", "SPCA", "IPSN", "RegTME")
  window_months <- 120L
  factor_rule <- "GR"  # 可改为 "ER" 或 "fixed"。
  fixed_m <- NULL      # factor_rule = "fixed" 时填写 1、2 或 3。
  seed <- 20261003L
  resume_dir <- NULL   # 单次续跑时填写已有 sp500_run_001 等文件夹。

  run_sensitivity <- TRUE  # TRUE：运行因子数敏感性实验。
  sensitivity_methods <- methods
  factor_grid <- c(1L, 2L, 3L)
  reuse_dir <- NULL    # 只重算部分方法时，填写含 GR、m1、m2、m3 的已完成结果目录，并相应修改 methods。
  reuse_gr <- TRUE
  gr_resume_dir <- NULL  # 敏感性实验续跑 GR 时填写已有 GR 运行目录；默认新建。
  pilot_cache_dir <- file.path(output_dir, "pilot_cache")

  if (run_sensitivity) {
    result <- run_sp500_sensitivity(
      train_file = train_file, test_file = test_file, output_dir = output_dir,
      methods = sensitivity_methods, factor_grid = factor_grid, reuse_dir = reuse_dir,
      reuse_gr = reuse_gr, gr_resume_dir = gr_resume_dir, pilot_cache_dir = pilot_cache_dir,
      start_date = start_date, end_date = end_date, window_months = window_months,
      seed = seed)
  } else {
    reuse_result_file <- poet_reuse_file(reuse_dir,
      if (factor_rule == "fixed") paste0("m", fixed_m) else factor_rule)
    result <- run_sp500(
      train_file = train_file, test_file = test_file, output_dir = output_dir,
      methods = methods, start_date = start_date, end_date = end_date,
      window_months = window_months, factor_rule = factor_rule, fixed_m = fixed_m,
      seed = seed, resume_dir = resume_dir, reuse_result_file = reuse_result_file,
      pilot_cache_dir = pilot_cache_dir)
  }
}
