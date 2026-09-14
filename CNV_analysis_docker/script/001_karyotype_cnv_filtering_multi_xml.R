#!/usr/bin/env Rscript

# ==============================================================================
# Multi-XML version of 001_karyotpe_cnv_filtering.R
# ==============================================================================
#
# Purpose
#   Same extraction/filtering logic as 001_karyotpe_cnv_filtering.R, but reads
#   bookmarks from SEVERAL GenomeStudio XML exports at once and writes ONE
#   combined pair of output CSVs instead of one pair per XML file.
#
#   Use this instead of 001_karyotpe_cnv_filtering.R when an analysis spans
#   more than one XML export (e.g. several chips/batches run together). For a
#   single XML, keep using 001_karyotpe_cnv_filtering.R / test_CNV.Rmd exactly
#   as before -- this script is purely additive, nothing about the single-XML
#   pipeline changes.
#
#   Every bookmark row keeps a `source_xml` column recording which XML file it
#   came from, so the combined CSVs remain fully traceable back to the
#   original per-chip exports even after merging.
#
#   Filtering (confidence/size/hotspot-overlap thresholds) is applied ONCE to
#   the combined table rather than once per XML and then merged. This gives
#   an identical result to filtering each XML separately and row-binding the
#   filtered outputs, because every filtering rule here is per-row (it never
#   compares one row to another), but it is simpler and avoids ever having
#   mismatched threshold settings between files of the same run.
#
# Usage
#   Rscript 001_karyotype_cnv_filtering_multi_xml.R --xml="file1.xml,file2.xml,..." \
#     [--comment_threshold=40] [--hotspot_size_threshold=100000] \
#     [--non_hotspot_size_threshold=350000] [--output_dir=.] [--combined_name=combined]
#
#   --xml (required) accepts EITHER a comma-separated inline list of XML
#   paths, OR the path to a plain-text file listing one XML path per line
#   (same convention as --samples in 002_build_CNV_sample_metadata.R).
#
#   --combined_name (optional, default "combined") sets the shared prefix of
#   the two output files: <combined_name>_extracted_data.csv and
#   <combined_name>_filtered_data.csv. Give it a run-specific name (e.g.
#   "--combined_name=batch1_batch2") if you keep several multi-XML runs in
#   the same output_dir and want to tell them apart.
#
#   Examples:
#     Rscript 001_karyotype_cnv_filtering_multi_xml.R \
#       --xml="input/EB1.xml,input/EB2.xml,input/EB3.xml" \
#       --comment_threshold=75 --hotspot_size_threshold=100000 \
#       --non_hotspot_size_threshold=350000 --output_dir=output/run1 \
#       --combined_name=combined
#
#     Rscript 001_karyotype_cnv_filtering_multi_xml.R --xml=xml_list.txt --output_dir=output/run1
#       (xml_list.txt: one XML path per line)
#
# ==============================================================================

library(dplyr)
library(xml2)
library(tools) # For file_path_sans_ext (kept for parity with 001; not used directly here)

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

# Pull one "--flag=value" out of an args vector; returns the value (or NA if
# absent) and the args vector with that flag removed. Identical helper to the
# one used in 002_build_CNV_sample_metadata.R, for a consistent CLI style.
extract_flag <- function(args, flag_prefix) {
  pattern <- paste0("^", flag_prefix)
  hit <- grep(pattern, args, value = TRUE)
  value <- if (length(hit)) sub(pattern, "", hit[1]) else NA_character_
  list(value = value, remaining = args[!grepl(pattern, args)])
}

step_xml  <- extract_flag(args, "--xml=")
step_ct   <- extract_flag(step_xml$remaining, "--comment_threshold=")
step_hs   <- extract_flag(step_ct$remaining, "--hotspot_size_threshold=")
step_nhs  <- extract_flag(step_hs$remaining, "--non_hotspot_size_threshold=")
step_out  <- extract_flag(step_nhs$remaining, "--output_dir=")
step_name <- extract_flag(step_out$remaining, "--combined_name=")

xml_arg <- step_xml$value
if (is.na(xml_arg) || !nzchar(trimws(xml_arg))) {
  stop(
    "No XML files supplied. Usage:\n",
    "  Rscript 001_karyotype_cnv_filtering_multi_xml.R --xml=\"file1.xml,file2.xml,...\" ",
    "[--comment_threshold=40] [--hotspot_size_threshold=100000] ",
    "[--non_hotspot_size_threshold=350000] [--output_dir=.] [--combined_name=combined]\n",
    "  --xml accepts a comma-separated inline list, or a path to a text file with one XML path per line.",
    call. = FALSE
  )
}

# --xml: comma-separated inline list, OR a text file with one path per line
# (same pattern as --samples in 002_build_CNV_sample_metadata.R).
if (file.exists(xml_arg)) {
  xml_filenames <- trimws(readLines(xml_arg, warn = FALSE))
} else {
  xml_filenames <- trimws(strsplit(xml_arg, ",", fixed = TRUE)[[1]])
}
xml_filenames <- unique(xml_filenames[nzchar(xml_filenames)])
if (!length(xml_filenames)) {
  stop("--xml was given but contained no usable XML path.", call. = FALSE)
}
missing_xml <- xml_filenames[!file.exists(xml_filenames)]
if (length(missing_xml)) {
  stop("XML file(s) not found: ", paste(missing_xml, collapse = ", "), call. = FALSE)
}

comment_threshold          <- if (!is.na(step_ct$value))   as.numeric(step_ct$value)   else 40
hotspot_size_threshold     <- if (!is.na(step_hs$value))   as.numeric(step_hs$value)   else 100000
non_hotspot_size_threshold <- if (!is.na(step_nhs$value))  as.numeric(step_nhs$value)  else 350000
output_dir    <- if (!is.na(step_out$value)  && nzchar(trimws(step_out$value)))  step_out$value  else getwd()
combined_name <- if (!is.na(step_name$value) && nzchar(trimws(step_name$value))) step_name$value else "combined"

# Make sure the output directory exists (creates it, and any missing parent
# directories, if it doesn't -- harmless if it already exists).
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("XML files to process (", length(xml_filenames), "):\n", sep = "")
cat(paste0("  - ", xml_filenames), sep = "\n")

# Define hotspot locations directly in the script -- identical table to
# 001_karyotpe_cnv_filtering.R (the source of truth); if that table is ever
# revised, mirror the change here too.
hotspot_locations <- data.frame(
  hot_chr = c(20, 12, 17, 1, 5, 18, 17, 7, 9, 11, 13, 1, 8, 14, 6, 15, 3, 22),
  hot_start = c(31260580, 11784484, 53204134, 172930860, 105164299, 58532768, 7307685, 132915240, 40965786, 2778770, 87047745, 15873505, 92287772, 18375541, 129978855, 66907662, 19121, 23925838),
  hot_end = c(32166810, 25403186, 54216532, 185830860, 118068501, 63932768, 8140855, 134100999, 87165786, 10678770, 101047745, 16748447, 126218590, 106573063, 138678855, 67007662, 26358322, 49228316)
)

# ------------------------------------------------------------------------------
# Step 1: parse every XML's bookmarks into one combined raw data frame, each
# row tagged with source_xml for traceability back to its original file.
# ------------------------------------------------------------------------------

read_one_xml <- function(xml_filename) {
  xml_data <- read_xml(xml_filename)
  bookmarks <- xml_find_all(xml_data, ".//bookmark")

  sample_ids <- c()
  chr_nums <- c()
  base_start_positions <- c()
  base_end_positions <- c()
  values <- c()
  comments <- c()

  for (bookmark in bookmarks) {
    sample_id <- xml_text(xml_find_first(bookmark, "sample_id"))
    chr_num <- xml_text(xml_find_first(bookmark, "chr_num"))
    base_start_pos <- xml_text(xml_find_first(bookmark, "base_start_pos"))
    base_end_pos <- xml_text(xml_find_first(bookmark, "base_end_pos"))
    value <- xml_text(xml_find_first(bookmark, "value"))
    comment_text <- xml_text(xml_find_first(bookmark, "comment"))
    comment_number <- as.numeric(sub(".*?(\\d+\\.\\d+).*", "\\1", comment_text))

    sample_ids <- c(sample_ids, sample_id)
    chr_nums <- c(chr_nums, chr_num)
    base_start_positions <- c(base_start_positions, base_start_pos)
    base_end_positions <- c(base_end_positions, base_end_pos)
    values <- c(values, value)
    comments <- c(comments, comment_number)
  }

  data.frame(
    sample_id = sample_ids,
    chr_num = chr_nums,
    base_start_pos = base_start_positions,
    base_end_pos = base_end_positions,
    value = values,
    comment = comments,
    source_xml = basename(xml_filename),
    stringsAsFactors = FALSE
  )
}

per_xml_counts <- integer(0)
extracted_list <- lapply(xml_filenames, function(f) {
  d <- read_one_xml(f)
  per_xml_counts[basename(f)] <<- nrow(d)
  d
})
extracted_data <- do.call(rbind, extracted_list)

if (!nrow(extracted_data)) {
  stop("No <bookmark> elements were found in any of the supplied XML files.", call. = FALSE)
}

cat("Bookmarks read per XML file:\n")
for (nm in names(per_xml_counts)) cat("  - ", nm, ": ", per_xml_counts[[nm]], "\n", sep = "")
cat("Total bookmarks combined: ", nrow(extracted_data), "\n", sep = "")

# Rename column before saving the extracted data (same as 001)
extracted_data <- extracted_data %>%
  rename(confidence = comment)

# Captured BEFORE the filtering steps below mutate/subset `extracted_data`,
# so the final summary can report the true combined raw bookmark count.
total_bookmarks_combined <- nrow(extracted_data)

# Duplicate-key sanity check ACROSS files: the same
# (sample_id, chr_num, base_start_pos, base_end_pos, value) bookmark appearing
# in more than one supplied XML is very likely the same physical bookmark
# exported twice (e.g. an XML listed twice by mistake, or overlapping export
# ranges) rather than two genuinely independent calls. Both/all copies are
# kept (nothing is dropped silently) but flagged so it can be checked by hand.
dup_key <- paste(
  extracted_data$sample_id, extracted_data$chr_num,
  extracted_data$base_start_pos, extracted_data$base_end_pos, extracted_data$value,
  sep = "|"
)
cross_file_dups <- duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE)
if (any(cross_file_dups)) {
  dup_report <- extracted_data[cross_file_dups, c("sample_id", "chr_num", "source_xml")]
  warning(
    sum(cross_file_dups), " bookmark row(s) are exact duplicates across the supplied XML files ",
    "(identical sample_id/chr_num/coordinates/value) -- all copies are kept in the combined output; ",
    "check whether the same XML was listed twice or export ranges overlapped:\n",
    paste(utils::capture.output(print(dup_report, row.names = FALSE)), collapse = "\n"),
    call. = FALSE
  )
}

# Save the combined extracted data, same write.csv2() convention (European
# CSV: semicolon fields, comma decimals) as 001_karyotpe_cnv_filtering.R, so
# it remains a drop-in --extracted input for 002_build_CNV_sample_metadata.R
# and a drop-in `extracted_input`/`historical_input` for 003.
extracted_out_path <- file.path(output_dir, paste0(combined_name, "_extracted_data.csv"))
write.csv2(extracted_data, extracted_out_path, row.names = FALSE, quote = TRUE, na = "")
print(paste("Combined extracted data has been successfully saved to", extracted_out_path))
str(extracted_data)

# ------------------------------------------------------------------------------
# Steps 2-7: same filtering logic as 001_karyotpe_cnv_filtering.R, applied
# once to the combined table (see header note on why this is equivalent to
# filtering per-XML then merging).
# ------------------------------------------------------------------------------

extracted_data <- extracted_data %>%
  mutate(base_end_pos = as.numeric(base_end_pos),
         base_start_pos = as.numeric(base_start_pos)) %>%
  filter(!is.na(base_end_pos) & !is.na(base_start_pos)) %>%
  mutate(CNV_size = base_end_pos - base_start_pos)

extracted_data$hotspot_region <- FALSE

extracted_data <- extracted_data %>%
  filter(confidence > comment_threshold)

extracted_data <- extracted_data %>%
  filter(!(value == 2 & chr_num == "X"))

for (i in seq_len(nrow(extracted_data))) {
  if (extracted_data$value[i] != 2) {
    chr_matches <- hotspot_locations %>%
      filter(hot_chr == extracted_data$chr_num[i])

    if (nrow(chr_matches) > 0) {
      for (j in seq_len(nrow(chr_matches))) {
        if (!is.na(extracted_data$base_start_pos[i]) && !is.na(extracted_data$base_end_pos[i]) &&
            !is.na(chr_matches$hot_start[j]) && !is.na(chr_matches$hot_end[j])) {
          if ((extracted_data$base_start_pos[i] > chr_matches$hot_start[j] & extracted_data$base_start_pos[i] < chr_matches$hot_end[j]) |
              (extracted_data$base_end_pos[i] > chr_matches$hot_start[j] & extracted_data$base_end_pos[i] < chr_matches$hot_end[j])) {
            extracted_data$hotspot_region[i] <- TRUE
            break
          }
        }
      }
    }
  }
}

filtered_data <- extracted_data %>%
  filter((hotspot_region == TRUE & CNV_size > hotspot_size_threshold) |
           (hotspot_region == FALSE & CNV_size > non_hotspot_size_threshold))

hotspot_true_count <- filtered_data %>%
  filter(hotspot_region == TRUE) %>%
  nrow()
print(paste("Number of rows with hotspot_region = TRUE:", hotspot_true_count))

filtered_out_path <- file.path(output_dir, paste0(combined_name, "_filtered_data.csv"))
write.csv2(filtered_data, filtered_out_path, row.names = FALSE, quote = TRUE, na = "")
print(paste("Combined filtered data has been successfully saved to", filtered_out_path))

cat("\nSummary\n")
cat("  XML files combined: ", length(xml_filenames), "\n", sep = "")
cat("  Total bookmarks (all XMLs): ", total_bookmarks_combined, "\n", sep = "")
cat("  Rows retained after confidence/size/hotspot filtering: ", nrow(filtered_data), "\n", sep = "")
cat("  Combined extracted CSV: ", extracted_out_path, "\n", sep = "")
cat("  Combined filtered CSV: ", filtered_out_path, "\n", sep = "")
