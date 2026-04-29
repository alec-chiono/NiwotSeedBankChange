# Self-contained posterior retrodictive check function for cmdstanr fits
# Adapted from Michael Betancourt's plot_hist_quantiles
# (https://github.com/betanalpha/mcmc_visualization_tools)
# License: BSD-3-Clause

#' Plot posterior retrodictive check histograms
#'
#' Compares observed data against posterior predictive distributions from a
#' cmdstanr fit, faceted by year × habitat. For each panel, nested quantile
#' ribbons of histogram bin counts across posterior draws are shown in red,
#' with the observed histogram overlaid in black.
#'
#' @param fit A `CmdStanMCMC` object from cmdstanr.
#' @param data_df A data.frame with the original observations (one row per
#'   species × plot × year × habitat).
#' @param predicted_var Name of the generated quantities variable in the Stan
#'   model. Must be indexed as `[N_years, N_habitats, N_species]`.
#' @param year_labels Character vector of year labels (length N_years).
#' @param habitat_labels Character vector of habitat labels (length N_habitats).
#' @param count_col Column name for the count variable in `data_df`.
#' @param year_col Column name for the year variable in `data_df`.
#' @param habitat_col Column name for the habitat variable in `data_df`.
#' @param species_col Column name for the species variable in `data_df`.
#' @param n_bins Number of histogram bins.
#' @param max_count Optional upper cap for predicted and observed counts. Values
#'   exceeding this are set to `max_count`. Useful for highly right-skewed data.
#' @param probs Quantile probabilities for the nested ribbons. Must be
#'   symmetric (e.g., 0.1–0.9) with odd length so the middle element is the
#'   median.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of producing a
#'   base R plot.
plot_prc_hist <- function(
  fit,
  data_df,
  predicted_var = "predicted_seed_counts",
  year_labels = NULL,

  habitat_labels = NULL,
  count_col = "count",
  year_col = "year",
  habitat_col = "habitat",
  species_col = "USDA_code",
  n_bins = 25,

  max_count = NULL,
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
) {
  # ── Betancourt color palette ──
  c_light           <- "#DCBCBC"
  c_light_highlight <- "#C79999"
  c_mid             <- "#B97C7C"
  c_mid_highlight   <- "#A25050"
  c_dark            <- "#8F2727"

  # ── Dimension labels ──
  years    <- sort(unique(data_df[[year_col]]))
  habitats <- sort(unique(data_df[[habitat_col]]))
  N_years    <- length(years)
  N_habitats <- length(habitats)

  if (is.null(year_labels))    year_labels    <- as.character(years)
  if (is.null(habitat_labels)) habitat_labels <- habitats

  # ── Extract posterior draws as matrix: [draw × (y*h*s)] ──
  draws_mat <- fit$draws(predicted_var, format = "matrix")
  N_draws   <- nrow(draws_mat)
  N_species <- ncol(draws_mat) / (N_years * N_habitats)

  # ── Observed counts summed across plots per species × year × habitat ──
  obs_summary <- data_df |>
    dplyr::group_by(
      .data[[year_col]], .data[[habitat_col]], .data[[species_col]]
    ) |>
    dplyr::summarize(total = sum(.data[[count_col]]), .groups = "drop")

  # ── Helper: configure bins ──
  configure_bins <- function(values, baseline, n_bins) {
    lo <- min(c(values, baseline))
    hi <- max(c(values, baseline))
    delta <- (hi - lo) / n_bins
    if (delta == 0) delta <- 1
    N <- (hi - lo) / delta
    excess <- N - floor(N)
    if (excess > 1e-15) {
      lo <- lo - 0.5 * delta * excess
      hi <- hi + 0.5 * delta * excess
    }
    list(lo = lo, hi = hi, delta = delta, breaks = seq(lo, hi, delta))
  }

  # ── Helper: bin plotting coords (step-function style) ──
  configure_bin_plotting <- function(breaks) {
    B    <- length(breaks) - 1
    idxs <- rep(seq_len(B), each = 2)
    xs   <- vapply(seq_along(idxs), function(b) {
      if (b %% 2 == 1) breaks[idxs[b]] else breaks[idxs[b] + 1]
    }, numeric(1))
    list(idxs = idxs, xs = xs)
  }

  # ── Set up plot layout ──
  old_par <- par(
    mfrow = c(N_years, N_habitats),
    mar = c(4, 4, 3, 1)
  )
  on.exit(par(old_par), add = TRUE)

  # ── Loop over year × habitat ──
  for (y in seq_len(N_years)) {
    for (h in seq_len(N_habitats)) {
      # Column names for this year × habitat slice
      col_names <- paste0(
        predicted_var, "[",
        y, ",", h, ",", seq_len(N_species), "]"
      )
      # pred_mat_cell: [N_draws × N_species]
      pred_mat_cell <- draws_mat[, col_names, drop = FALSE]

      # Observed values for this cell (one per species)
      obs_cell <- obs_summary |>
        dplyr::filter(
          .data[[year_col]]    == years[y],
          .data[[habitat_col]] == habitats[h]
        ) |>
        dplyr::pull(total)

      # Apply truncation if requested
      pred_vals <- as.vector(pred_mat_cell)
      if (!is.null(max_count)) {
        pred_vals <- pmin(pred_vals, max_count)
        pred_mat_cell[] <- pmin(as.vector(pred_mat_cell), max_count)
        obs_cell <- pmin(obs_cell, max_count)
      }

      # Configure bins
      bc <- configure_bins(pred_vals, obs_cell, n_bins)
      breaks <- bc$breaks
      B      <- length(breaks) - 1

      pc <- configure_bin_plotting(breaks)
      plot_idxs <- pc$idxs
      plot_xs   <- pc$xs

      # ── Bin counts per posterior draw ──
      bin_counts_mat <- matrix(NA_real_, nrow = N_draws, ncol = B)
      for (b in seq_len(B)) {
        bin_counts_mat[, b] <- rowSums(
          pred_mat_cell >= breaks[b] & pred_mat_cell < breaks[b + 1]
        )
      }
      # Last bin is inclusive on the right
      bin_counts_mat[, B] <- bin_counts_mat[, B] + rowSums(
        pred_mat_cell == breaks[B + 1]
      )

      # ── Quantiles of bin counts across draws ──
      quantiles <- apply(bin_counts_mat, 2, quantile, probs = probs)
      # quantiles: [length(probs) × B]

      # Expand to step-function coordinates
      plot_quantiles <- quantiles[, plot_idxs, drop = FALSE]

      # ── Observed bin counts ──
      obs_counts <- vapply(seq_len(B), function(b) {
        sum(obs_cell >= breaks[b] & obs_cell < breaks[b + 1])
      }, numeric(1))
      obs_counts[B] <- obs_counts[B] +
        sum(obs_cell == breaks[B + 1])
      plot_obs <- obs_counts[plot_idxs]

      # ── y-axis range ──
      ylim <- c(0, max(c(quantiles[length(probs), ], obs_counts)))

      # ── Plot ──
      title <- paste0(year_labels[y], " \u00D7 ", habitat_labels[h])
      xlab_text <- if (!is.null(max_count)) {
        paste0("Seed count (capped at ", max_count, ")")
      } else {
        "Seed count"
      }

      plot(
        1, type = "n", main = title,
        xlim = c(bc$lo, bc$hi), xlab = xlab_text,
        ylim = ylim, ylab = "Number of species"
      )

      # Nested quantile ribbons (10-90%, 20-80%, 30-70%, 40-60%)
      polygon(
        c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[1, ], rev(plot_quantiles[9, ])),
        col = c_light, border = NA
      )
      polygon(
        c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[2, ], rev(plot_quantiles[8, ])),
        col = c_light_highlight, border = NA
      )
      polygon(
        c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[3, ], rev(plot_quantiles[7, ])),
        col = c_mid, border = NA
      )
      polygon(
        c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[4, ], rev(plot_quantiles[6, ])),
        col = c_mid_highlight, border = NA
      )
      # Median line segments
      for (b in seq_len(B)) {
        idx1 <- 2 * b - 1
        idx2 <- 2 * b
        lines(
          plot_xs[idx1:idx2],
          plot_quantiles[5, idx1:idx2],
          col = c_dark, lwd = 2
        )
      }

      # Observed histogram overlaid
      lines(plot_xs, plot_obs, col = "white", lty = 1, lwd = 4)
      lines(plot_xs, plot_obs, col = "black", lty = 1, lwd = 2)
    }
  }

  invisible(NULL)
}
