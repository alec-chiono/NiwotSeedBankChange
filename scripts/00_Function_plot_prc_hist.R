# Self-contained posterior retrodictive check functions for cmdstanr fits
# Adapted from Michael Betancourt's plot_hist_quantiles
# (https://github.com/betanalpha/mcmc_visualization_tools)
# License: BSD-3-Clause
#
# All functions return a ggplot2/patchwork object; save with ggsave().
# Requires: ggplot2, patchwork, dplyr

library(ggplot2)
library(patchwork)
library(dplyr)

# ── Shared palette and helpers ────────────────────────────────────────────────

.cols <- c(
  light = "#DCBCBC",
  light_highlight = "#C79999",
  mid = "#B97C7C",
  mid_highlight = "#A25050",
  dark = "#8F2727"
)

.betancourt_theme <- function() {
  theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      plot.title = element_text(size = 11)
    )
}

# Build a tidy data frame of quantile ribbon rectangles + median segments
# for one set of posterior draws and one observed vector.
#
# Returns list(ribbons, medians, observed) — data frames ready for ggplot layers.
.prc_bin_data <- function(
  pred_mat,
  obs_vals,
  n_bins,
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
) {
  # Breaks
  all_vals <- c(as.vector(pred_mat), obs_vals)
  lo <- min(all_vals, na.rm = TRUE)
  hi <- max(all_vals, na.rm = TRUE)
  # Guard: if all values identical (e.g. everything capped or all zero),
  # expand the range by ±0.5 so bins have non-zero width
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
    ctr <- if (is.finite(lo)) lo else 0
    lo <- ctr - 0.5
    hi <- ctr + 0.5
  }
  delta <- (hi - lo) / n_bins
  if (!is.finite(delta) || delta <= 0) {
    delta <- 1
  }
  excess <- ((hi - lo) / delta) - floor((hi - lo) / delta)
  if (excess > 1e-15) {
    lo <- lo - 0.5 * delta * excess
    hi <- hi + 0.5 * delta * excess
  }
  breaks <- seq(lo, hi, by = delta)
  B <- max(length(breaks) - 1L, 1L)

  # Per-draw bin counts matrix [N_draws × B]
  bcm <- matrix(NA_real_, nrow(pred_mat), B)
  for (b in seq_len(B)) {
    bcm[, b] <- rowSums(pred_mat >= breaks[b] & pred_mat < breaks[b + 1L])
  }
  bcm[, B] <- bcm[, B] + rowSums(pred_mat == breaks[B + 1L])

  # Quantiles [length(probs) × B]
  quants <- apply(bcm, 2L, quantile, probs = probs)

  # Ribbon data frame (one row per bin per band)
  bands <- list(
    list(ymin_p = 1L, ymax_p = 9L, fill = .cols[["light"]]),
    list(ymin_p = 2L, ymax_p = 8L, fill = .cols[["light_highlight"]]),
    list(ymin_p = 3L, ymax_p = 7L, fill = .cols[["mid"]]),
    list(ymin_p = 4L, ymax_p = 6L, fill = .cols[["mid_highlight"]])
  )
  ribbons <- do.call(
    rbind,
    lapply(bands, function(bnd) {
      data.frame(
        xmin = breaks[seq_len(B)],
        xmax = breaks[seq_len(B) + 1L],
        ymin = quants[bnd$ymin_p, ],
        ymax = quants[bnd$ymax_p, ],
        fill = bnd$fill,
        stringsAsFactors = FALSE
      )
    })
  )
  ribbons$fill <- factor(
    ribbons$fill,
    levels = unname(.cols[c(
      "light",
      "light_highlight",
      "mid",
      "mid_highlight"
    )])
  )

  # Median segments (one row per bin)
  medians <- data.frame(
    x = breaks[seq_len(B)],
    xend = breaks[seq_len(B) + 1L],
    y = quants[5L, ],
    yend = quants[5L, ]
  )

  # Observed bin counts
  obs_bc <- vapply(
    seq_len(B),
    function(b) sum(obs_vals >= breaks[b] & obs_vals < breaks[b + 1L]),
    numeric(1L)
  )
  obs_bc[B] <- obs_bc[B] + sum(obs_vals == breaks[B + 1L])
  observed <- data.frame(
    x = breaks[seq_len(B)],
    xend = breaks[seq_len(B) + 1L],
    y = obs_bc,
    yend = obs_bc
  )

  list(
    ribbons = ribbons,
    medians = medians,
    observed = observed,
    breaks = breaks,
    lo = lo,
    hi = hi,
    ylim_top = max(c(quants[9L, ], obs_bc), na.rm = TRUE)
  )
}

# Build a single ggplot panel from .prc_bin_data output
.prc_panel <- function(bd, title, xlab, ylab = "Number of species") {
  # Build a staircase path for the observed histogram so both horizontal
  # and vertical connectors between bins are drawn.
  # x alternates: left edge of bin 1, then for each subsequent break the
  # right edge appears twice (once closing the current bin, once opening
  # the next), finished by the right edge of the last bin.
  # y follows: 0 baseline cap, count for each bin repeated twice, 0 cap.
  B <- nrow(bd$observed)
  obs_step <- data.frame(
    x = rep(bd$breaks[seq_len(B + 1L)], each = 2L),
    y = c(0, rep(bd$observed$y, each = 2L), 0)
  )

  ggplot() +
    geom_rect(
      data = bd$ribbons,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      colour = NA
    ) +
    scale_fill_identity() +
    geom_segment(
      data = bd$medians,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour = .cols[["dark"]],
      linewidth = 0.9
    ) +
    geom_path(
      data = obs_step,
      aes(x = x, y = y),
      colour = "black",
      linewidth = 0.75,
      linejoin = "mitre"
    ) +
    coord_cartesian(xlim = c(bd$lo, bd$hi), ylim = c(0, bd$ylim_top * 1.05)) +
    labs(title = title, x = xlab, y = ylab) +
    .betancourt_theme()
}

# ── Extract b_veg posterior summaries ────────────────────────────────────────

.extract_bveg_summaries <- function(fit, USDA_lookup) {
  bveg_draws <- fit$draws("b_veg", format = "matrix")
  var_names <- colnames(bveg_draws)

  do.call(
    rbind,
    lapply(seq_len(nrow(USDA_lookup)), function(ki) {
      code <- USDA_lookup$USDA_code[ki]
      id <- as.integer(USDA_lookup$USDA_code_id[ki])
      do.call(
        rbind,
        lapply(
          list(
            list(h_idx = 1L, hab = "mesic"),
            list(h_idx = 2L, hab = "xeric")
          ),
          function(hh) {
            vname <- sprintf("b_veg[%d,%d]", id, hh$h_idx)
            if (!vname %in% var_names) {
              return(NULL)
            }
            d <- bveg_draws[, vname]
            data.frame(
              USDA_code = code,
              habitat = hh$hab,
              median = median(d),
              lo90 = quantile(d, 0.05),
              hi90 = quantile(d, 0.95),
              lo50 = quantile(d, 0.25),
              hi50 = quantile(d, 0.75),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    })
  )
}

.ols_slopes <- function(vegf) {
  do.call(
    rbind,
    by(vegf, list(vegf$USDA_code, vegf$habitat), function(df) {
      if (nrow(df) < 3 || var(df$year) == 0) {
        return(NULL)
      }
      b <- coef(lm(count_scaled ~ year, data = df))[["year"]]
      data.frame(
        USDA_code = df$USDA_code[1],
        habitat = df$habitat[1],
        ols_slope = b,
        stringsAsFactors = FALSE
      )
    })
  )
}


# ── plot_prc_hist2() ──────────────────────────────────────────────────────────
#' @return A patchwork object (year × habitat grid). Save with ggsave().

plot_prc_hist2 <- function(
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
  years <- sort(unique(data_df[[year_col]]))
  habitats <- sort(unique(data_df[[habitat_col]]))
  N_years <- length(years)
  N_habs <- length(habitats)

  if (is.null(year_labels)) {
    year_labels <- as.character(years)
  }
  if (is.null(habitat_labels)) {
    habitat_labels <- as.character(habitats)
  }

  draws_mat <- fit$draws("pred_count", format = "matrix")

  xlab_text <- if (!is.null(max_count)) {
    paste0("Seed count (capped at ", max_count, ")")
  } else {
    "Seed count"
  }

  panels <- vector("list", N_years * N_habs)
  k <- 1L
  for (y in seq_len(N_years)) {
    for (h in seq_len(N_habs)) {
      row_idx <- which(
        data_df[[year_col]] == years[y] &
          data_df[[habitat_col]] == habitats[h]
      )
      title_str <- paste0(year_labels[y], " \u00D7 ", habitat_labels[h])

      if (length(row_idx) == 0L) {
        panels[[k]] <- ggplot() + labs(title = title_str) + .betancourt_theme()
      } else {
        pm <- draws_mat[, row_idx, drop = FALSE]
        obs <- data_df[[count_col]][row_idx]
        if (!is.null(max_count)) {
          pm[] <- pmin(as.vector(pm), max_count)
          obs <- pmin(obs, max_count)
        }
        bd <- .prc_bin_data(pm, obs, n_bins, probs)
        panels[[k]] <- .prc_panel(
          bd,
          title = title_str,
          xlab = xlab_text,
          ylab = "Number of species"
        )
      }
      k <- k + 1L
    }
  }

  wrap_plots(panels, nrow = N_years, ncol = N_habs)
}


# ── plot_bveg_vs_ols() ────────────────────────────────────────────────────────
#' @return A patchwork of two ggplot panels (mesic | xeric). Save with ggsave().

plot_bveg_vs_ols <- function(fit, vegf, USDA_lookup) {
  bveg <- .extract_bveg_summaries(fit, USDA_lookup)
  ols <- .ols_slopes(vegf)
  merged <- merge(bveg, ols, by = c("USDA_code", "habitat"))

  make_panel <- function(hab) {
    d <- merged[merged$habitat == hab, ]
    lims <- range(c(d$ols_slope, d$lo90, d$hi90), na.rm = TRUE)
    title <- paste0("Veg slopes \u00D7 ", tools::toTitleCase(hab))

    ggplot(d, aes(x = ols_slope, y = median)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      geom_segment(
        aes(xend = ols_slope, y = lo90, yend = hi90),
        colour = .cols[["light_highlight"]],
        linewidth = 0.8
      ) +
      geom_segment(
        aes(xend = ols_slope, y = lo50, yend = hi50),
        colour = .cols[["mid"]],
        linewidth = 1.5
      ) +
      geom_point(colour = .cols[["dark"]], size = 1.8) +
      coord_cartesian(xlim = lims, ylim = lims) +
      labs(
        title = title,
        x = "OLS slope (observed)",
        y = "Posterior median b_veg"
      ) +
      .betancourt_theme()
  }

  make_panel("mesic") + make_panel("xeric")
}


# ── plot_seed_prc() ───────────────────────────────────────────────────────────
#' @return A single ggplot panel. Save with ggsave().

plot_seed_prc <- function(
  fit,
  seedf,
  seed_divisor = 2023 - 1989,
  n_bins = 25,
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
) {
  seed_draws <- fit$draws("seed_change_rep", format = "matrix")
  obs_change <- seedf$change / seed_divisor

  bd <- .prc_bin_data(seed_draws, obs_change, n_bins, probs)
  .prc_panel(
    bd,
    title = paste0("Seed change (n=", nrow(seedf), ")"),
    xlab = "Seed change (per year)",
    ylab = "Count"
  )
}
