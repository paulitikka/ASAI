# Replot Figures 1-4 from the bundled source-data CSVs.
# Command line, from the package folder: Rscript plot_figures.R
# RStudio: source('plot_figures.R'); plot_asai_figures('C:/path/to/package')
# The Python pipeline calculates the full-precision PCA coordinates and source data.

plot_asai_figures <- function(root = '.') {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop("Install ggplot2 first: install.packages('ggplot2')")
  }
  root <- normalizePath(root, mustWork = TRUE)
  src <- file.path(root, 'source_data')
  dest <- file.path(root, 'outputs', 'figures_R')
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  read_source <- function(name) read.csv(file.path(src, name), check.names = FALSE)
  colour <- '#477B8A'
  base_theme <- ggplot2::theme_classic(base_size = 10, base_family = 'sans') +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = 'bold'))
  save_plot <- function(p, name, width, height) {
    ggplot2::ggsave(file.path(dest, paste0(name, '.pdf')), p,
                   width = width, height = height, device = grDevices::pdf,
                   useDingbats = FALSE)
    ggplot2::ggsave(file.path(dest, paste0(name, '.png')), p,
                   width = width, height = height, dpi = 600, bg = 'white')
    ggplot2::ggsave(file.path(dest, paste0(name, '.svg')), p,
                   width = width, height = height, device = grDevices::svg)
  }

  d1 <- read_source('Figure1_histogram_bins.csv')
  d1$aa <- factor(d1$aa, levels = sort(unique(d1$aa)))
  p1 <- ggplot2::ggplot(d1, ggplot2::aes(bin_mid, density)) +
    ggplot2::geom_col(width = 0.025, fill = colour) +
    ggplot2::facet_wrap(~aa, ncol = 5, scales = 'free_y') +
    ggplot2::scale_x_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1, 1.5),
                               expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(x = 'Residue relative solvent accessibility (RSA)', y = 'Probability density') +
    base_theme
  save_plot(p1, 'Figure1_residue_RSA', 7.1, 5.5)

  d2 <- read_source('Figure2_PC2_dispersion.csv')
  line <- read_source('Figure2_fit_line.csv')
  r <- stats::cor(d2$sd_chain_mean_RSA, d2$PC2, method = 'pearson')
  p <- stats::cor.test(d2$sd_chain_mean_RSA, d2$PC2, method = 'pearson')$p.value
  p2 <- ggplot2::ggplot(d2, ggplot2::aes(sd_chain_mean_RSA, PC2)) +
    ggplot2::geom_line(data = line, ggplot2::aes(sd_chain_mean_RSA, fitted_PC2),
                       inherit.aes = FALSE, colour = 'grey50', linewidth = 0.4) +
    ggplot2::geom_point(colour = colour, size = 2.3) +
    ggplot2::geom_text(ggplot2::aes(x = label_x, y = label_y, label = aa), size = 2.6) +
    ggplot2::labs(x = 'SD of chain-specific mean RSA', y = 'PC2 coordinate',
                  subtitle = sprintf('Pearson r = %.3f; P = %.3g; n = 20 types', r, p)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.17)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.16)) + base_theme
  save_plot(p2, 'Figure2_PC2_dispersion', 4.2, 3.7)

  d3 <- read_source('Figure3_divergence.csv')
  d3 <- d3[order(d3$interval_order), ]
  d3$interval_label <- factor(d3$interval_label, levels = d3$interval_label)
  p3 <- ggplot2::ggplot(d3, ggplot2::aes(interval_label, CV_across_types, group = 1)) +
    ggplot2::geom_line(colour = colour, linewidth = 0.5) +
    ggplot2::geom_point(colour = colour, size = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 3.12)) +
    ggplot2::labs(x = 'Chain-specific mean RSA interval (%)',
                  y = 'Coefficient of variation across types') + base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_plot(p3, 'Figure3_divergence', 7.1, 3.3)

  d4 <- read_source('Figure4_kingdom_coordinates.csv')
  panels <- do.call(rbind, lapply(1:2, function(k) {
    data.frame(aa = d4$aa, component = paste0(if (k == 1) 'a  ' else 'b  ', 'PC', k),
               x = d4[[paste0('PC', k, '_prok')]], y = d4[[paste0('PC', k, '_euk')]],
               label_x = d4[[paste0('label_x_PC', k)]],
               label_y = d4[[paste0('label_y_PC', k)]])
  }))
  # Common x/y limits preserve a 1:1 comparison in both panels.
  limits <- range(c(panels$x, panels$y, panels$label_x, panels$label_y))
  padding <- diff(limits) * 0.12
  limits <- limits + c(-padding, padding)
  p4 <- ggplot2::ggplot(panels, ggplot2::aes(x, y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 'dashed',
                         colour = 'grey60', linewidth = 0.35) +
    ggplot2::geom_point(colour = colour, size = 2.1) +
    ggplot2::geom_text(ggplot2::aes(x = label_x, y = label_y, label = aa), size = 2.6) +
    ggplot2::facet_wrap(~component, nrow = 1) +
    ggplot2::coord_fixed(xlim = limits, ylim = limits) +
    ggplot2::labs(x = 'Prokaryotic component coordinate', y = 'Eukaryotic component coordinate') +
    base_theme
  save_plot(p4, 'Figure4_kingdom_coordinates', 7.1, 3.6)
  invisible(list(Figure1 = p1, Figure2 = p2, Figure3 = p3, Figure4 = p4))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  plot_asai_figures(if (length(args)) args[[1L]] else '.')
}
