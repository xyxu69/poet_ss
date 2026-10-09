# options(poet.autorun=FALSE) 可只载入函数。
.poet_factor_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.poet_factor_source_file)) {
  .poet_factor_cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(.poet_factor_cli) == 1L) {
    .poet_factor_cli <- sub("^--file=", "", .poet_factor_cli)
    if (tolower(basename(.poet_factor_cli)) == "factor_number_simulation.r")
      .poet_factor_source_file <- .poet_factor_cli
  }
}
.poet_factor_dir <- if (!is.null(.poet_factor_source_file))
  dirname(normalizePath(.poet_factor_source_file, winslash = "/", mustWork = TRUE)) else getwd()
for (.poet_package in c("mvtnorm", "SpatialNP", "ICSNP")) {
  if (!requireNamespace(.poet_package, quietly = TRUE)) stop("Missing R package: ", .poet_package)
}
rm(.poet_package)

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

poet_select_sign_factors <- function(Y, kmax = 10L) {
  n <- nrow(Y)
  p <- ncol(Y)
  limit <- min(n, p) - 1L
  kmax <- min(as.integer(kmax), limit - 1L)
  if (kmax < 1L) stop("Factor selection requires min(n,p) >= 3")
  center <- ICSNP::spatial.median(Y)
  S <- p * SpatialNP::SCov(Y, location = center)
  values <- pmax(eigen(S, symmetric = TRUE, only.values = TRUE)$values[seq_len(limit)], 0)
  floor <- max(values[1], 1) * .Machine$double.eps
  values <- pmax(values, floor)
  tail_sum <- rev(cumsum(rev(values)))
  indices <- seq_len(kmax)
  er <- values[indices] / values[indices + 1L]
  gr <- log1p(values[indices] / tail_sum[indices]) /
        log1p(values[indices + 1L] / tail_sum[indices + 1L])
  c(ER = which.max(er), GR = which.max(gr))
}

# source_file 仅兼容旧调用并记录来源；计算使用本文件内的生成器。
run_factor_number_simulation <- function(source_file = NULL,
                                         output_dir = file.path(.poet_factor_dir, "results", "factor"),
                                         dimensions = c(200L, 300L, 400L, 500L),
                                         n = 250L, m = 3L, sims = 1000L,
                                         scenarios = 1:4, kmax = 10L, seed = 988L) {
  stopifnot(all(dimensions > m), m >= 1L, m <= 3L, sims >= 1L,
            n > m + 1L, all(scenarios %in% 1:4), kmax >= m)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  .poet_index <- 1L
  repeat {
    run_dir <- file.path(output_dir, sprintf("factor_run_%03d", .poet_index))
    if (dir.create(run_dir, showWarnings = FALSE)) break
    if (!file.exists(run_dir)) stop("Cannot create output directory: ", run_dir)
    .poet_index <- .poet_index + 1L
  }
  settings <- list(n = n, m = m, dimensions = dimensions, sims = sims,
                   scenarios = scenarios, kmax = kmax, seed = seed)
  saveRDS(list(settings = settings,
               source_file_reference = source_file, generator = "embedded",
               session_info = sessionInfo()), file.path(run_dir, "run_info.rds"))
  set.seed(seed)
  distributions <- c("normal", "t4", "t2.2", "mixnorm")
  all_rows <- list()
  for (scenario in scenarios) for (p in dimensions) {
    checkpoint_dir <- file.path(run_dir, "replicates", paste0("mi", scenario, "_p", p))
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    selections <- matrix(NA_integer_, nrow = sims, ncol = 2L, dimnames = list(NULL, c("ER", "GR")))
    for (replicate in seq_len(sims)) {
      variances <- c(1, 0.75^2, 0.5^2)[seq_len(m)]
      B <- sapply(variances, function(v) rnorm(p, sd = sqrt(v)))
      B <- matrix(B, nrow = p, ncol = m)
      generated <- generate_factor_data_elliptical(B, n = n, dist = distributions[scenario])
      selections[replicate, ] <- poet_select_sign_factors(generated$data, kmax)
      saveRDS(list(selection = selections[replicate, ], scenario = scenario, p = p,
                   completed = replicate, rng_state = .Random.seed),
              file.path(checkpoint_dir, sprintf("rep_%04d.rds", replicate)))
      if (replicate %% 50L == 0L || replicate == sims)
        cat(sprintf("factor mi=%d p=%d %d/%d saved\n", scenario, p, replicate, sims))
    }
    for (selector in colnames(selections)) {
      counts <- tabulate(selections[, selector], nbins = kmax)
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        scenario = scenario, p = p, n = n, true_m = m, selector = selector,
        selected_m = seq_len(kmax), count = counts, frequency = counts / sims)
    }
  }
  table <- do.call(rbind, all_rows)
  saveRDS(list(table = table, settings = settings, session_info = sessionInfo()),
          file.path(run_dir, "factor_number_results.rds"))
  tryCatch(write.csv(table, file.path(run_dir, "factor_number_frequencies.csv"), row.names = FALSE),
           error = function(e) warning("RDS saved; CSV failed: ", conditionMessage(e)))
  cat("SAVED: ", normalizePath(run_dir, winslash = "/"), "\n", sep = "")
  invisible(list(table = table, run_dir = run_dir))
}

if (isTRUE(getOption("poet.autorun", TRUE))) {
  # 参数区：按需修改，然后直接 Source 或全选运行本文件。
  n <- 250L
  m <- 3L
  dimensions <- c(200L, 300L, 400L, 500L)
  sims <- 1000L
  scenarios <- 1:4
  kmax <- 10L
  seed <- 988L
  output_dir <- file.path(.poet_factor_dir, "results", "factor")

  result <- run_factor_number_simulation(output_dir = output_dir,
    dimensions = dimensions, n = n, m = m, sims = sims,
    scenarios = scenarios, kmax = kmax, seed = seed)
}
