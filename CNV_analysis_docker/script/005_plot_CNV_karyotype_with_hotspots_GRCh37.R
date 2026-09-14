#!/usr/bin/env Rscript

# ==============================================================================
# Genome-wide CNV plot with recurrent hotspot regions overlaid (GRCh37)
# ==============================================================================
#
# Purpose
#   Same genome-wide plot as plot_CNV_karyotype_from_filtered_GRCh37.R (each
#   sample on its own row, each retained CNV drawn across its full genomic
#   interval), PLUS the 20 recurrent hPSC hotspot regions used by
#   003_karyotype_cnv_filtering_GRCh37_conf75_validated_chrX2filtered.R drawn
#   as shaded background bands, so it's visually obvious which retained CNVs
#   fall inside a known recurrent region and which don't.
#
# Expected inputs
#   1) The filtered calls CSV, same as the other plotting script:
#        cnvPartition_3.2.0_filtered_data_GRCh37_corrected_conf75_chrX2filtered.csv
#      Required columns: sample_id, chr_num, base_start_pos, base_end_pos, value
#   2) The hotspot reference CSV written by script 003 (one row per region):
#        cnvPartition_3.2.0_hotspot_reference_GRCh37.csv
#      Required columns: hotspot_rank, hotspot_chr, hotspot_start, hotspot_end
#
# Input format
#   Both inputs are the European CSV format produced by write.csv2():
#     field separator = ";"
#     decimal separator = ","
#
# Usage
#   Rscript 005_plot_CNV_karyotype_with_hotspots_GRCh37.R \
#     [filtered_calls.csv] [hotspot_reference.csv] [output.pdf] [sample_selection]
#
#   sample_selection (optional, 4th argument): restrict the plot to a subset
#   of samples instead of plotting all of them (the hotspot bands are always
#   drawn in full regardless of this selection -- they are reference regions,
#   not tied to any one sample). Accepts EITHER a comma-separated list of
#   sample_id and/or sample_name values (any mix of the two -- each token is
#   matched against both), OR the path to a plain-text file listing one value
#   per line (same convention as --samples in 002_build_CNV_sample_metadata.R).
#   Leave empty/omit to plot every sample, exactly as before.
#     Rscript 005_....R filtered.csv hotspots.csv out.pdf "W1.6,K1.7"
#
# If no arguments are supplied, the defaults below are used.
#
# Output
#   Editable vector PDF.
#
# Genome assembly
#   GRCh37 / hg19
# ==============================================================================

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(viridis)

# ------------------------------------------------------------------------------
# Input / output
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[1]
} else {
  "cnvPartition_3.2.0_filtered_data_GRCh37_corrected_conf75_chrX2filtered.csv"
}

hotspot_file <- if (length(args) >= 2) {
  args[2]
} else {
  "cnvPartition_3.2.0_hotspot_reference_GRCh37.csv"
}

output_pdf <- if (length(args) >= 3) {
  args[3]
} else {
  "CNV_karyotype_plot_with_hotspots_GRCh37.pdf"
}

if (!file.exists(input_file)) stop("Filtered calls file not found: ", input_file, call. = FALSE)
if (!file.exists(hotspot_file)) stop("Hotspot reference file not found: ", hotspot_file, call. = FALSE)

# ------------------------------------------------------------------------------
# Read filtered calls and hotspot reference table
#
# read.csv2() correctly handles the European CSV format used by the filtering
# pipeline: semicolon-separated fields and comma decimal separators.
# ------------------------------------------------------------------------------

cnv <- read.csv2(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hotspots <- read.csv2(
  hotspot_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------------------------

required_columns <- c(
  "sample_id",
  "chr_num",
  "base_start_pos",
  "base_end_pos",
  "value"
)

missing_columns <- setdiff(required_columns, colnames(cnv))
if (length(missing_columns) > 0) {
  stop(
    "Filtered calls file is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

required_hotspot_columns <- c("hotspot_rank", "hotspot_chr", "hotspot_start", "hotspot_end")
missing_hotspot_columns <- setdiff(required_hotspot_columns, colnames(hotspots))
if (length(missing_hotspot_columns) > 0) {
  stop(
    "Hotspot reference file is missing required column(s): ",
    paste(missing_hotspot_columns, collapse = ", ")
  )
}

if (nrow(cnv) == 0) {
  stop("The filtered input contains no retained CNV calls.")
}
if (nrow(hotspots) == 0) {
  stop("The hotspot reference file contains no regions.")
}

# `sample_id` is deliberately used by NAME rather than by numeric column index,
# even though it is column B in the locked-filter output. This makes the script
# safer if additional audit columns are later inserted or reordered.
cnv$sample_id <- trimws(as.character(cnv$sample_id))

if (any(is.na(cnv$sample_id) | cnv$sample_id == "")) {
  stop("At least one retained call has a missing or blank sample_id.")
}

# ------------------------------------------------------------------------------
# Y-axis label: sample_name when available, sample_id otherwise
#
# Grouping/ordering of CNV calls onto rows is still done by `sample_id`
# (guaranteed unique and non-blank, checked above); `sample_name` is only used
# for what gets PRINTED on the axis, so a missing/blank name for one sample
# never breaks the plot for the others. If two different sample_id values
# happen to share the same sample_name, the sample_id is appended to keep
# them visually distinguishable.
# ------------------------------------------------------------------------------

if ("sample_name" %in% names(cnv)) {
  cnv$sample_name <- trimws(as.character(cnv$sample_name))
} else {
  cnv$sample_name <- NA_character_
}

sample_key_to_label <- cnv %>%
  distinct(sample_id, sample_name) %>%
  mutate(
    display_label = ifelse(is.na(sample_name) | sample_name == "", sample_id, sample_name)
  )
dup_labels <- sample_key_to_label$display_label[duplicated(sample_key_to_label$display_label)]
if (length(dup_labels)) {
  sample_key_to_label <- sample_key_to_label %>%
    mutate(
      display_label = ifelse(
        display_label %in% dup_labels,
        paste0(display_label, " (", sample_id, ")"),
        display_label
      )
    )
}
cnv <- cnv %>%
  left_join(sample_key_to_label %>% select(sample_id, display_label), by = "sample_id")

# ------------------------------------------------------------------------------
# Optional sample selection: restrict the plot to a subset of samples
#
# sample_selection (4th CLI argument) lets you plot only some of the retained
# samples instead of all of them. It accepts EITHER a comma-separated list of
# sample_id and/or sample_name values (any mix of the two, each token matched
# against both), OR the path to a plain-text file listing one value per line
# (same convention as --samples in 002_build_CNV_sample_metadata.R). Leave it
# empty/omit it to plot every sample, exactly as before this feature existed.
# The hotspot reference bands are NOT filtered by this selection -- they are
# fixed genomic regions, independent of which samples are shown.
# ------------------------------------------------------------------------------

sample_selection_arg <- if (length(args) >= 4 && nzchar(trimws(args[4]))) trimws(args[4]) else NA_character_

if (!is.na(sample_selection_arg)) {
  if (file.exists(sample_selection_arg)) {
    requested_samples <- trimws(readLines(sample_selection_arg, warn = FALSE))
  } else {
    requested_samples <- trimws(strsplit(sample_selection_arg, ",", fixed = TRUE)[[1]])
  }
  requested_samples <- unique(requested_samples[nzchar(requested_samples)])
  if (!length(requested_samples)) {
    stop("A sample selection was supplied but contained no usable sample_id/sample_name.", call. = FALSE)
  }

  # Match each requested token against sample_id first, then against the
  # display label (sample_name, already disambiguated above) -- either form
  # is accepted so you can copy-paste straight from a previous plot's axis.
  matched_ids <- unique(sample_key_to_label$sample_id[
    sample_key_to_label$sample_id %in% requested_samples |
    sample_key_to_label$display_label %in% requested_samples
  ])
  matched_tokens <- unique(requested_samples[
    requested_samples %in% sample_key_to_label$sample_id |
    requested_samples %in% sample_key_to_label$display_label
  ])
  unmatched_tokens <- setdiff(requested_samples, matched_tokens)

  if (length(unmatched_tokens)) {
    stop(
      "Sample selection contains value(s) that match neither a sample_id nor a sample_name ",
      "in this filtered CSV: ", paste(unmatched_tokens, collapse = ", "), "\n",
      "Available sample_id values: ", paste(sample_key_to_label$sample_id, collapse = ", "), "\n",
      "Available sample_name values: ", paste(sample_key_to_label$display_label, collapse = ", "),
      call. = FALSE
    )
  }

  cnv <- cnv %>% filter(sample_id %in% matched_ids)
  sample_key_to_label <- sample_key_to_label %>% filter(sample_id %in% matched_ids)

  if (!nrow(cnv)) {
    stop(
      "After applying the sample selection, no retained CNV calls remain to plot for: ",
      paste(matched_ids, collapse = ", "),
      call. = FALSE
    )
  }

  message("Sample selection applied: plotting ", length(matched_ids), " of the samples in the filtered CSV.")
}

# ------------------------------------------------------------------------------
# GRCh37 chromosome lengths
# ------------------------------------------------------------------------------

chr_info <- data.frame(
  chr_num = c(as.character(1:22), "X", "Y"),
  chr_length = c(
    249250621, 243199373, 198022430, 191154276,
    180915260, 171115067, 159138663, 146364022,
    141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753,
    81195210, 78077248, 59128983, 63025520,
    48129895, 51304566, 155270560, 59373566
  ),
  stringsAsFactors = FALSE
)

chr_info <- chr_info %>%
  mutate(
    chr_offset = lag(cumsum(chr_length), default = 0),
    chr_end_global = chr_offset + chr_length,
    chr_mid_global = chr_offset + chr_length / 2
  )

# ------------------------------------------------------------------------------
# Normalize chromosome names (CNV calls)
# ------------------------------------------------------------------------------

cnv <- cnv %>%
  mutate(
    chr_num = toupper(
      sub(
        "^chr",
        "",
        trimws(as.character(chr_num)),
        ignore.case = TRUE
      )
    )
  )

# GenomeStudio may occasionally report a pseudoautosomal `XY` chromosome label.
# Such a call cannot be positioned unambiguously on this X/Y axis and therefore
# must be reviewed rather than silently converted.
if (any(cnv$chr_num == "XY", na.rm = TRUE)) {
  stop(
    "The filtered input contains chromosome 'XY' calls. ",
    "Review these calls before plotting; they are not automatically assigned to X or Y."
  )
}

unknown_chr <- setdiff(unique(cnv$chr_num), chr_info$chr_num)
if (length(unknown_chr) > 0) {
  stop(
    "Unrecognized chromosome value(s) in filtered calls: ",
    paste(unknown_chr, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Normalize chromosome names (hotspot regions) and validate coordinates
# ------------------------------------------------------------------------------

hotspots <- hotspots %>%
  mutate(
    hotspot_chr = toupper(
      sub("^chr", "", trimws(as.character(hotspot_chr)), ignore.case = TRUE)
    ),
    hotspot_start = as.numeric(hotspot_start),
    hotspot_end = as.numeric(hotspot_end)
  )

unknown_hotspot_chr <- setdiff(unique(hotspots$hotspot_chr), chr_info$chr_num)
if (length(unknown_hotspot_chr) > 0) {
  stop(
    "Unrecognized chromosome value(s) in hotspot reference: ",
    paste(unknown_hotspot_chr, collapse = ", ")
  )
}

if (any(is.na(hotspots$hotspot_start)) || any(is.na(hotspots$hotspot_end))) {
  stop("Some hotspot coordinates could not be interpreted as numeric.")
}
if (any(hotspots$hotspot_end < hotspots$hotspot_start)) {
  stop("At least one hotspot region has end coordinate < start coordinate.")
}

# ------------------------------------------------------------------------------
# Validate coordinates / values (CNV calls)
# ------------------------------------------------------------------------------

cnv$base_start_pos <- as.numeric(cnv$base_start_pos)
cnv$base_end_pos   <- as.numeric(cnv$base_end_pos)
cnv$value          <- as.numeric(cnv$value)

if (any(is.na(cnv$base_start_pos)) || any(is.na(cnv$base_end_pos))) {
  stop("Some CNV coordinates could not be interpreted as numeric.")
}
if (any(cnv$base_start_pos < 1)) {
  stop("At least one CNV has a start coordinate < 1.")
}
if (any(cnv$base_end_pos < cnv$base_start_pos)) {
  stop("At least one CNV has end coordinate < start coordinate.")
}
if (any(is.na(cnv$value))) {
  stop("At least one CNV value could not be interpreted as numeric.")
}

# ------------------------------------------------------------------------------
# Convert chromosome-local coordinates to one continuous genome-wide axis
# ------------------------------------------------------------------------------

cnv <- cnv %>%
  left_join(
    chr_info %>% select(chr_num, chr_length, chr_offset),
    by = "chr_num"
  ) %>%
  mutate(
    genome_start = chr_offset + base_start_pos,
    genome_end   = chr_offset + base_end_pos
  )

if (any(cnv$base_end_pos > cnv$chr_length)) {
  stop("At least one CNV extends beyond the GRCh37 chromosome length.")
}

hotspots <- hotspots %>%
  left_join(
    chr_info %>% select(chr_num, chr_length, chr_offset),
    by = c("hotspot_chr" = "chr_num")
  ) %>%
  mutate(
    genome_start = chr_offset + hotspot_start,
    genome_end   = chr_offset + hotspot_end,
    genome_mid   = (genome_start + genome_end) / 2,
    hotspot_label = paste0("#", hotspot_rank)
  )

if (any(hotspots$hotspot_end > hotspots$chr_length)) {
  stop("At least one hotspot region extends beyond the GRCh37 chromosome length.")
}

# ------------------------------------------------------------------------------
# Minimum visual width
#
# On a linear genome-wide axis (~3.1 Gb across the whole plot width), a CNV of
# 100-350 kb -- exactly the size range this pipeline is built to detect around
# recurrent hotspots -- is only 0.003-0.01% of the axis and renders as a
# sub-pixel sliver, effectively invisible. Several of the 20 hotspot regions
# themselves are similarly small (down to ~100 kb). This does NOT change any
# underlying data (the CSV/audit/hotspot-reference files still record the true
# coordinates/size); it only guarantees every retained call and every hotspot
# band is drawn at least min_visual_width_bp wide on this specific plot, so
# nothing is lost to visual scale. Features already wider than this are drawn
# at their true size, unchanged.
total_genome_length <- sum(chr_info$chr_length)
min_visual_width_bp <- total_genome_length * 0.0035  # ~0.35% of genome width

cnv <- cnv %>%
  mutate(
    genome_mid = (genome_start + genome_end) / 2,
    display_half_width = pmax(genome_end - genome_start, min_visual_width_bp) / 2,
    display_start = pmax(chr_offset, genome_mid - display_half_width),
    display_end   = pmin(chr_offset + chr_length, genome_mid + display_half_width)
  )

hotspots <- hotspots %>%
  mutate(
    display_half_width = pmax(genome_end - genome_start, min_visual_width_bp) / 2,
    display_start = pmax(chr_offset, genome_mid - display_half_width),
    display_end   = pmin(chr_offset + chr_length, genome_mid + display_half_width)
  )

# ------------------------------------------------------------------------------
# Sample order
#
# Samples are ordered by first appearance in the filtered file. Thus the user can
# control the biological display order simply by ordering sample blocks in the
# input table before running this plotting script.
# ------------------------------------------------------------------------------

sample_order <- unique(cnv$sample_id)
n_samples <- length(sample_order)
# Axis labels follow the same order as sample_order, but show display_label
# (sample_name, falling back to sample_id) instead of the raw sample_id.
label_order <- sample_key_to_label$display_label[match(sample_order, sample_key_to_label$sample_id)]

cnv <- cnv %>%
  mutate(
    sample_y = n_samples - match(sample_id, sample_order) + 1
  )

# ------------------------------------------------------------------------------
# CNV-state labels
#
# `value = 2` is labeled cnLOH, consistent with the cnvPartition retained-call
# interpretation used in the locked filtering workflow.
# ------------------------------------------------------------------------------

cnv <- cnv %>%
  mutate(
    CNV_class = case_when(
      value == 0 ~ "CN 0",
      value == 1 ~ "CN 1",
      value == 2 ~ "cnLOH",
      value == 3 ~ "CN 3",
      value == 4 ~ "CN 4",
      TRUE       ~ paste0("CN ", value)
    )
  )

known_levels <- c("CN 0", "CN 1", "cnLOH", "CN 3", "CN 4")
observed_levels <- unique(cnv$CNV_class)
legend_levels <- c(
  known_levels[known_levels %in% observed_levels],
  setdiff(observed_levels, known_levels)
)

cnv$CNV_class <- factor(cnv$CNV_class, levels = legend_levels)

# ------------------------------------------------------------------------------
# Plot
#
# Hotspot regions are drawn FIRST as full-height shaded background bands (so
# CNV calls are always drawn on top of them, never hidden), with a small
# rank-number label above each band. Each CNV call is then a rectangle
# spanning its complete start-to-end interval on its sample's row, exactly as
# in the non-hotspot plotting script.
# ------------------------------------------------------------------------------

band_half_height <- 0.22
hotspot_band_color <- "#F4A259"   # warm orange, distinct from the viridis CNV fill scale

p <- ggplot() +

  # Recurrent hotspot regions -- full-height shaded background bands, widened
  # to at least min_visual_width_bp where needed so small regions (down to
  # ~100 kb) stay visible; true coordinates are unchanged in the CSV.
  geom_rect(
    data = hotspots,
    aes(xmin = display_start, xmax = display_end, ymin = 0.5, ymax = n_samples + 0.5),
    fill = hotspot_band_color,
    alpha = 0.18,
    color = NA
  ) +

  # Hotspot rank labels along the top edge
  geom_text(
    data = hotspots,
    aes(x = genome_mid, y = n_samples + 0.5, label = hotspot_label),
    color = hotspot_band_color,
    size = 2.3,
    vjust = -0.4,
    angle = 90,
    hjust = 0
  ) +

  # Light sample-row guides
  geom_hline(
    yintercept = seq_len(n_samples),
    linewidth = 0.3,
    color = "grey90"
  ) +

  # Chromosome boundaries
  geom_vline(
    data = chr_info[-nrow(chr_info), ],
    aes(xintercept = chr_end_global),
    linewidth = 0.25,
    color = "grey75"
  ) +

  # Retained CNV intervals -- widened to at least min_visual_width_bp where
  # needed so small calls stay visible; true coordinates are unchanged in the
  # CSV/audit outputs.
  geom_rect(
    data = cnv,
    aes(
      xmin = display_start,
      xmax = display_end,
      ymin = sample_y - band_half_height,
      ymax = sample_y + band_half_height,
      fill = CNV_class
    ),
    color = NA
  ) +

  # Chromosome labels
  scale_x_continuous(
    breaks = chr_info$chr_mid_global,
    labels = chr_info$chr_num,
    expand = expansion(mult = c(0.002, 0.03))
  ) +

  # sample_name labels on the Y axis (falls back to sample_id when a name is
  # missing; see display_label above)
  scale_y_continuous(
    breaks = seq_len(n_samples),
    labels = rev(label_order),
    limits = c(0.5, n_samples + 1.3),
    expand = c(0, 0)
  ) +

  # Cividis is perceptually uniform and color-blind safe
  scale_fill_viridis_d(
    option = "C",
    begin = 0.08,
    end = 0.92,
    direction = 1,
    name = "CNV state"
  ) +

  labs(
    x = "Chromosome",
    y = NULL,
    caption = paste0(
      "Shaded bands = ", nrow(hotspots),
      " recurrent hPSC hotspot regions (GRCh37; numbers = hotspot_rank from cnvPartition_3.2.0_hotspot_reference_GRCh37.csv). ",
      "Features narrower than ~", round(min_visual_width_bp / 1e6, 1),
      " Mb are widened on this plot to stay visible; true coordinates/size are unchanged in the underlying CSV files."
    )
  ) +

  theme_classic(base_size = 10) +

  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.4),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.caption = element_text(size = 7, color = "grey40", hjust = 0),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# ------------------------------------------------------------------------------
# Export editable vector PDF
# ------------------------------------------------------------------------------

plot_height <- max(3.5, 1.5 + 0.42 * n_samples)

# Make sure the destination directory exists (ggsave/cairo_pdf silently
# reroute to the current working directory otherwise).
output_dir_for_pdf <- dirname(output_pdf)
if (!dir.exists(output_dir_for_pdf)) {
  dir.create(output_dir_for_pdf, recursive = TRUE, showWarnings = FALSE)
}

# NOTE: we deliberately do NOT use ggsave(..., device = cairo_pdf, ...) here.
# In some R/ggplot2 version combinations, passing `device` as a bare function
# causes ggsave() to fail to forward `filename` correctly to the underlying
# device -- the destination *directory* ends up right, but the file is saved
# under cairo_pdf()'s own default name ("Rplots.pdf") instead of the intended
# name. Opening the device explicitly and printing the plot to it sidesteps
# that argument-forwarding entirely and always writes to the exact path given.
#
# Cairo produces a slightly nicer vector PDF (better font/transparency
# handling), but requires the system Cairo/X11 libraries -- on macOS that
# means XQuartz must be installed, which is often not the case. Rather than
# require that, we try cairo_pdf() first and, if it is not usable on this
# machine, fall back automatically to base R's pdf() device, which is also
# fully vector-based (editable in Illustrator/Affinity/Inkscape) and needs no
# external system library at all.
open_pdf_device <- function(filename, width, height) {
  # Probe Cairo availability on a throwaway temp file, NOT on the real output
  # file: drawing a test rectangle directly on `filename` would otherwise
  # leave a spurious blank first page in front of the real plot once we
  # print() it (ggplot2's print method starts a fresh page, so the probe's
  # page and the real plot's page would both end up in the same PDF).
  probe_file <- tempfile(fileext = ".pdf")
  cairo_ok <- tryCatch({
    grDevices::cairo_pdf(filename = probe_file, width = 2, height = 2)
    # cairo_pdf() can "succeed" (no error) while silently failing to load its
    # native library, in which case the device it opened is unusable. Force a
    # trivial draw so a real Cairo failure surfaces here, inside tryCatch.
    grid::grid.rect(gp = grid::gpar(col = NA, fill = NA))
    grDevices::dev.off()
    TRUE
  }, error = function(e) {
    if (!is.null(grDevices::dev.list())) grDevices::dev.off()
    FALSE
  }, warning = function(w) {
    if (!is.null(grDevices::dev.list())) grDevices::dev.off()
    FALSE
  })
  unlink(probe_file)

  if (cairo_ok) {
    grDevices::cairo_pdf(filename = filename, width = width, height = height, onefile = TRUE)
  } else {
    message("cairo_pdf() is not usable on this R installation (Cairo/X11 not available) -- using base R's pdf() device instead.")
    grDevices::pdf(file = filename, width = width, height = height, onefile = TRUE)
  }
}

open_pdf_device(output_pdf, width = 14.5, height = plot_height)
print(p)
invisible(grDevices::dev.off())

if (!file.exists(output_pdf)) {
  stop("PDF was not created at the expected path: ", output_pdf, call. = FALSE)
}

message("Plotting complete.")
message("Filtered calls input: ", input_file)
message("Hotspot reference input: ", hotspot_file)
message("Samples: ", n_samples)
message("Retained CNV calls plotted: ", nrow(cnv))
message("Hotspot regions plotted: ", nrow(hotspots))
message("Vector PDF: ", output_pdf)
