# dev/plot_holdout_report.R
#
# Histogram of the per-configuration means in a saved holdout report, contenders against seeds.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/plot_holdout_report.R                                  # the default report
#   Rscript dev/plot_holdout_report.R state/holdoutReportX.txt out.png
#
# Input is section 2 of `Rscript dev/holdout_test.R --report`, saved as text: the printed tibble
# of label / group / n / mean / sd. Lines that are not a tibble row are ignored, so the whole
# report can be saved and passed in unedited.

suppressMessages(library(tidyverse))
here::i_am("dev/plot_holdout_report.R")

args <- commandArgs(trailingOnly = TRUE)
in_path  <- if (length(args) >= 1) args[1] else here::here("state", "holdoutReport2vSVREC.txt")
out_path <- if (length(args) >= 2) args[2] else sub("[.][^.]*$", ".png", in_path)

# A tibble row prints as: <row number> <label> <group> <n> <mean> <sd>. The header, the type
# row (<chr> <chr> ...) and any surrounding report text do not match.
rx <- "^\\s*\\d+\\s+(\\S+)\\s+(contender|seed)\\s+(\\d+)\\s+(-?[\\d.]+)\\s+(-?[\\d.]+)\\s*$"
rows <- readLines(in_path) |> stringr::str_match(rx)
rows <- rows[!is.na(rows[, 1]), , drop = FALSE]
if (!nrow(rows)) stop("no configuration rows found in ", in_path)

cfg <- tibble::tibble(label = rows[, 2],
                      group = factor(rows[, 3], levels = c("contender", "seed")),
                      n     = as.integer(rows[, 4]),
                      mean  = as.numeric(rows[, 5]),
                      sd    = as.numeric(rows[, 6]))

# Stacked, not overlaid: with 20 configurations against 5, a stack shows both the shape of the
# contender distribution and where the seeds sit within it. One bin per 0.025 of accuracy.
p <- ggplot2::ggplot(cfg, ggplot2::aes(x = mean, fill = group)) +
  ggplot2::geom_histogram(binwidth = 0.025, boundary = 0, colour = "white", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = c(contender = "#2C6FB5", seed = "#D98E32"), name = NULL) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_width(1), expand = c(0, 0, 0.05, 0)) +
  ggplot2::labs(x = "Mean accuracy", y = "Configurations") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "top")

ggplot2::ggsave(out_path, p, width = 6, height = 4, dpi = 300)
cat(nrow(cfg), " configurations (", sum(cfg$group == "contender"), " contender, ",
    sum(cfg$group == "seed"), " seed) -> ", out_path, "\n", sep = "")
