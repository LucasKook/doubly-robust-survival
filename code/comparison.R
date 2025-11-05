### Comparison of survGCM against DRLRT
### LK 2025

set.seed(1234)
save <- TRUE

### Dependencies

library("comets")
library("survival")
library("tidyverse")
library("ggpubr")

### Re-implementation of DLRLT
lrm <- comets:::lrm
drlrt <- function(Y, X, Z, B = 499, modelX = "lrm", ...) {
  N <- NROW(Y)
  m0 <- survival::coxph(Y ~ Z)
  compute_lr <- function(Y, X, Z) {
    m1 <- survival::coxph(Y ~ X + Z)
    2 * (stats::logLik(m1) - stats::logLik(m0))
  }
  LR <- compute_lr(Y, X, Z)
  mX <- do.call(modelX, c(list(y = X, x = Z), list(...)))
  rX <- stats::residuals(mX, response = X, data = Z)
  smpl <- sapply(seq_len(B), \(b) {
    newx <- stats::predict(mX, data = Z) + sample(rX, size = N, replace = TRUE)
    compute_lr(Y, newx, Z)
  })
  pv <- (1 + sum(LR < smpl)) / (1 + B)
  df <- c("df" = NCOL(X))

  structure(
    list(
      statistic = c("LR" = LR), p.value = pv, parameter = df,
      alternative = "X and Y are dependent given Z",
      data.name = paste0(deparse(match.call()), collapse = "\n"),
      method = "Doubly Robust Cox Partial Likelihood Ratio Test"
    ),
    class = c("drlrt", "htest")
  )
}

### DGP
dgp <- function(n = 1000, censoring_rate = 0.1,
                x_type = c("continuous", "binary", "MX2"),
                beta = 0) {
  x_type <- match.arg(x_type)

  # Generate continuous confounder Z
  Z <- rnorm(n)

  # Generate X influenced by Z (with some noise)
  X <- switch(x_type,
    "continuous" = 0.5 * Z + rnorm(n, sd = 0.5),
    "MX2" = rgamma(n, 2, exp(0.5 * Z)),
    "binary" = {
      prob_X <- plogis(0.8 * Z) # probability depends on Z
      rbinom(n, size = 1, prob = prob_X)
    }
  )

  # Generate survival times affected by Z but not X
  # Weibull distribution for baseline hazard
  lambda <- 0.1 # scale parameter
  rho <- 2 # shape parameter

  # Hazard depends on Z but not X
  hazard_ratio_Z <- exp(0.8 * Z + beta * X)
  true_times <- (-log(runif(n)) / (lambda * hazard_ratio_Z))^(1 / rho)

  # Generate censoring times affected by Z and X
  censoring_times <- rexp(n, rate = censoring_rate * exp(0.5 * Z + 0.5 * X))

  # Determine observed times and event indicators
  observed_times <- pmin(true_times, censoring_times)
  event <- as.numeric(true_times <= censoring_times)

  # Create Surv object
  Y <- survival::Surv(time = observed_times, event = event)

  # Return data.frame
  data.frame(Z = Z, X = X, Y = Y)
}

### Run
nsim <- 100
pb <- txtProgressBar(0, nsim, style = 3, width = 60)
res <- lapply(c(0, 0.5), \(tb) {
  lapply(seq_len(nsim), \(iter) {
    setTxtProgressBar(pb, iter)
    dat <- dgp(n = 300, x_type = "continuous", beta = tb)
    tic <- Sys.time()
    DRLRT <- drlrt(dat$Y, dat$X, dat$Z, B = 199)$p.value
    toc <- Sys.time()
    time_DRLRT <- toc - tic
    tic <- Sys.time()
    survGCM <- comets(Y ~ X | Z,
      data = dat, reg_YonZ = "cox",
      reg_XonZ = "lrm"
    )$p.value
    toc <- Sys.time()
    time_survGCM <- toc - tic
    tic <- Sys.time()
    naive <- comets(Y[, 2] ~ X | Z,
      data = dat, reg_XonZ = "lrm"
    )$p.value
    toc <- Sys.time()
    time_naive <- toc - tic
    c(
      beta = tb, "DRLRT" = DRLRT, "survGCM" = survGCM, "naive" = naive,
      "time_DRLRT" = time_DRLRT, "time_survGCM" = time_survGCM, "time_naive" = time_naive
    )
  }) |> dplyr::bind_rows()
}) |> dplyr::bind_rows()


### Plot
p1 <- res |>
  pivot_longer(
    cols = DRLRT:naive, names_to = "method",
    values_to = "p.value"
  ) |>
  ggplot(aes(x = p.value, color = method)) +
  facet_wrap(~beta, labeller = label_bquote(beta == .(beta))) +
  stat_ecdf() +
  labs(x = "p-value", y = "ECDF", tag = "A") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  theme(text = element_text(size = 13.5)) +
  scale_color_brewer(palette = "Dark2", labels = c("DRLRT" = "Yang et al.", "naive" = "binary GCM", "survGCM" = "TRAM-GCM"))

p2 <- res |>
  pivot_longer(
    cols = time_DRLRT:time_naive, names_to = "method",
    values_to = "runtime"
  ) |>
  ggplot(aes(x = runtime, color = method)) +
  stat_ecdf(pad = FALSE) +
  coord_flip() +
  labs(y = "relative rank", x = "runtime in seconds", tag = "B") +
  theme_bw() +
  scale_x_log10() +
  theme(text = element_text(size = 13.5)) +
  scale_color_brewer(palette = "Dark2", labels = c("time_DRLRT" = "Yang et al.", "time_naive" = "binary GCM", "time_survGCM" = "TRAM-GCM"))

ggarrange(p1, p2, common.legend = TRUE, widths = c(3, 2))

### Save
if (save) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  if (!dir.exists("./figures")) {
    dir.create("./figures")
  }
  saveRDS(res, "./results/results.rds")
  ggsave("./figures/comparison.pdf", height = 3.5, width = 11)
}
