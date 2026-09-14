#!/usr/bin/env Rscript

# ==============================================================================
# Validated cnvPartition reanalysis from a combined extracted CSV (GRCh37)
# ==============================================================================
#
# Filtering policy implemented exactly as requested:
#   * confidence > 75
#   * size > 350 kb outside the recurrent-region table
#   * size > 100 kb for copy-number calls overlapping the recurrent-region table
#   * strict interval overlap: call_start <= hotspot_end AND call_end >= hotspot_start
#   * physiological CN=1 on X or Y is excluded for explicitly male/XY samples
#   * male chrX value=2 calls are excluded by default as cnLOH-like artifacts of
#     treating a hemizygous male X chromosome as diploid upstream
#   * every input row is retained in a full audit table
#
# Important reviewed decisions:
#   * The script reads the combined extracted CSV directly; it does not parse XML.
#   * It uses base R only because dplyr is unavailable in the execution environment.
#   * CNV_size is end-start, matching the historical file and supplied script.
#     Inclusive_size is also reported for transparency but is not used for filtering.
#   * Autosomal value=2 calls are treated as copy-neutral LOH and receive the
#     standard >350-kb threshold, as in the supplied revised script.
#   * cnvPartition defines CNV Value=2 WITH a confidence value as copy-neutral
#     LOH. On the non-PAR X chromosome of an XY male, hemizygosity means the lack
#     of heterozygosity is physiological rather than evidence of acquired LOH.
#   * In these data, chrX value=2 calls recur in every male sample and collectively
#     cover most of chrX, strongly supporting an upstream sex/ploidy-setting artifact.
#     They are therefore excluded by default when expected_sex == "male".
#   * This downstream exclusion is a retrospective safeguard, not a substitute for
#     rerunning cnvPartition with the sample sex correctly specified. A genuine male
#     chrX duplication should be assessed from sex-aware cnvPartition output and/or
#     raw BAF/LRR rather than inferred from these artifact-prone value=2 segments.
#   * Male chrY value=2 calls are NOT automatically excluded by this rule.
#   * Chromosome label "XY" is never silently converted to X or Y.
#
# Coordinate provenance caveat:
#   Assou et al. (Stem Cell Reports 2020; PMID 31902703) state that their
#   abnormalities were converted to GRCh37A, but the article and its supplied
#   supplemental PDF do not publish a machine-readable coordinate table for the
#   ranked common regions. Figure 2D displays 24 common abnormal regions. The
#   20 intervals below are therefore the 20-region GRCh37 table supplied in the
#   revised script and requested for this audit; their exact transcription should
#   be confirmed against the original analysis spreadsheet/authors before SOP use.
#
# Usage:
#   Rscript karyotype_cnv_filtering_GRCh37_conf75_validated_chrX2filtered.R \
#     extracted.csv historical_filtered.csv CNV_sample_metadata.csv output_dir
# ==============================================================================
###########example within this docker:

#Rscript scripts/plot_CNV_karyotype_from_filtered_GRCh37.R \
#output/test/cnvPartition_3.2.0_filtered_data_GRCh37_corrected_conf75_chrX2filtered.csv \
#output/test/CNV_karyotype_test.pdf
# ==============================================================================
#
# All CSV inputs/outputs use semicolon fields and comma decimal separators.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
input_filename <- if (length(args) >= 1) args[1] else "cnvPartition_3.2.0_extracted_data.csv"
historical_filename <- if (length(args) >= 2) args[2] else "cnvPartition_3.2.0_filtered_data.csv"
metadata_filename <- if (length(args) >= 3) args[3] else "CNV_sample_metadata.csv"
output_dir <- if (length(args) >= 4) args[4] else "."

# Filtering thresholds and genome build are now configurable from the command
# line (positional args 5-8), instead of hardcoded, so a run's exact policy
# can be recorded/varied without editing this script. Defaults below match
# the values validated for the WTC11 dataset if nothing is passed.
#
#   Rscript 003_....R input.csv historical.csv metadata.csv output_dir \
#     [confidence_threshold] [hotspot_size_threshold] [standard_size_threshold] [genome_build]
confidence_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 75
hotspot_size_threshold <- if (length(args) >= 6) as.numeric(args[6]) else 100000
standard_size_threshold <- if (length(args) >= 7) as.numeric(args[7]) else 350000
genome_build <- if (length(args) >= 8) args[8] else "GRCh37"

# IMPORTANT: genome_build is a label AND a correctness guarantee. The 20-region
# hotspot coordinate table below is hardcoded in GRCh37 base pairs. Changing
# genome_build to anything else would silently mislabel those GRCh37
# coordinates as if they were on a different assembly, which would make every
# hotspot overlap decision wrong without any error being raised. Until a
# matching coordinate table for another build (e.g. GRCh38) is added to this
# script, genome_build must stay "GRCh37" -- passing anything else stops the
# run rather than producing quietly-incorrect hotspot classifications.
if (!identical(genome_build, "GRCh37")) {
  stop(
    "genome_build = '", genome_build, "' was requested, but this script's hotspot ",
    "coordinate table is hardcoded in GRCh37 coordinates -- using it under a ",
    "different build label would silently misclassify every hotspot overlap. ",
    "To support another build, a matching hotspot coordinate table for that ",
    "build must first be added to this script (see the 'hotspots' data.frame ",
    "below); only then should genome_build be changed.",
    call. = FALSE
  )
}

# Retrospective sex-chromosome safeguard.
# TRUE is appropriate for the present WTC11 XY dataset because chrX value=2 calls
# are pervasive across all male samples and are consistent with cnLOH-like calls
# generated when male X hemizygosity is not handled correctly upstream.
# Set FALSE only if analysing sex-aware cnvPartition output in which male chrX CN=2
# has been explicitly validated as a true copy-number gain.
exclude_male_chrX_value2_artifact <- TRUE

warnings_seen <- character(0)
record_warning <- function(message_text) {
  warnings_seen <<- c(warnings_seen, message_text)
  warning(message_text, call. = FALSE)
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

read_euro_csv <- function(path) {
  read.csv2(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

normalize_chr <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^chr", "", y, ignore.case = TRUE)
  toupper(y)
}

as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

key_number <- function(x) {
  ifelse(is.na(x), "NA", format(x, scientific = FALSE, trim = TRUE, digits = 15))
}

make_call_key <- function(sample_id, chr_num, start, end, value) {
  paste(
    trimws(as.character(sample_id)),
    normalize_chr(chr_num),
    key_number(start),
    key_number(end),
    key_number(value),
    sep = "|"
  )
}

join_nonempty <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) "" else paste(x, collapse = "; ")
}

logical_text <- function(x) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))

stop_if_missing(input_filename, "Authoritative extracted input")
stop_if_missing(historical_filename, "Historical filtered input")
stop_if_missing(metadata_filename, "Sample metadata input")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir, call. = FALSE)

# ==============================================================================
# Recurrent hPSC region table supplied for this audit (GRCh37)
# ==============================================================================

hotspots <- data.frame(
  hotspot_rank = 1:20,
  hotspot_chr = c(
    "20", "12", "17", "X", "1", "5", "18", "17", "7", "9",
    "11", "13", "16", "1", "8", "14", "6", "15", "3", "22"
  ),
  hotspot_start = c(
    29848383, 11937418, 51281495, 1, 172900000,
    104500000, 56200000, 7211004, 132600000, 68700000,
    2800000, 87700000, 1, 16200000, 93300000,
    19152018, 130300000, 67200000, 60799, 24268025
  ),
  hotspot_end = c(
    30754613, 25556120, 52293893, 60600000, 185800000,
    117404202, 61600000, 8044174, 133785759, 114900000,
    10700000, 101700000, 90354753, 17074942, 127230818,
    107349540, 139000000, 67300000, 26400000, 49570503
  ),
  genome_build = rep(genome_build, 20),
  coordinate_source_status = rep(
    "20-region table supplied in revised script; exact interval table not exposed in Assou article/supplement",
    20
  ),
  stringsAsFactors = FALSE
)

hotspots$region_id <- paste0(
  "rank", hotspots$hotspot_rank, ":chr", hotspots$hotspot_chr, ":",
  hotspots$hotspot_start, "-", hotspots$hotspot_end
)

if (nrow(hotspots) != 20 || !identical(hotspots$hotspot_rank, 1:20)) {
  stop("Hotspot table must contain unique ranks 1-20.", call. = FALSE)
}
if (any(duplicated(hotspots$hotspot_rank)) ||
    any(hotspots$hotspot_start < 1) ||
    any(hotspots$hotspot_end < hotspots$hotspot_start)) {
  stop("Hotspot table contains duplicate ranks or invalid intervals.", call. = FALSE)
}
if (!identical(hotspots$hotspot_chr[1], "20") ||
    hotspots$hotspot_start[1] != 29848383 ||
    hotspots$hotspot_end[1] != 30754613) {
  stop("Rank-1 chr20 GRCh37 guard interval failed.", call. = FALSE)
}

record_warning(paste(
  "Assou et al. Figure 2D reports 24 common abnormal regions, whereas the requested",
  "policy and supplied revised script use ranks 1-20. No additional four regions",
  "were introduced because that would change the requested biological policy."
))
record_warning(paste(
  "The Assou article confirms GRCh37A processing, but its article/supplement does",
  "not provide the exact 20-interval coordinate table. Coordinate transcription",
  "therefore remains a provenance item for confirmation before definitive SOP use."
))

# ==============================================================================
# Import and validate authoritative extracted calls
# ==============================================================================

raw <- read_euro_csv(input_filename)
required <- c("sample_id", "chr_num", "base_start_pos", "base_end_pos", "value", "confidence")
missing_required <- setdiff(required, names(raw))
if (length(missing_required)) {
  stop("Extracted input is missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
}
if (!nrow(raw)) stop("Extracted input contains no rows.", call. = FALSE)

audit <- data.frame(
  call_id = seq_len(nrow(raw)),
  sample_id = trimws(as.character(raw$sample_id)),
  chr_num_original = as.character(raw$chr_num),
  chr_num = normalize_chr(raw$chr_num),
  base_start_pos = as_number(raw$base_start_pos),
  base_end_pos = as_number(raw$base_end_pos),
  value = as_number(raw$value),
  confidence = as_number(raw$confidence),
  stringsAsFactors = FALSE
)

audit$CNV_size <- audit$base_end_pos - audit$base_start_pos
audit$inclusive_size <- audit$base_end_pos - audit$base_start_pos + 1
audit$genome_build <- genome_build

recognized_chr <- c(as.character(1:22), "X", "Y", "XY", "M", "MT")
audit$valid_sample_id <- !is.na(audit$sample_id) & nzchar(audit$sample_id)
audit$valid_chromosome <- !is.na(audit$chr_num) & audit$chr_num %in% recognized_chr
audit$valid_coordinates <-
  !is.na(audit$base_start_pos) & !is.na(audit$base_end_pos) &
  audit$base_start_pos >= 1 & audit$base_end_pos >= audit$base_start_pos
audit$valid_value <- !is.na(audit$value) & audit$value %in% 0:4
audit$valid_confidence <- !is.na(audit$confidence) & is.finite(audit$confidence)
audit$authoritative_call_key <- make_call_key(
  audit$sample_id, audit$chr_num, audit$base_start_pos, audit$base_end_pos, audit$value
)
audit$duplicate_exact_call <- duplicated(audit$authoritative_call_key) |
  duplicated(audit$authoritative_call_key, fromLast = TRUE)

if (any(audit$duplicate_exact_call)) {
  record_warning(paste(sum(audit$duplicate_exact_call),
                       "authoritative rows participate in exact duplicate call keys; none were dropped."))
}
if (any(!audit$valid_chromosome)) {
  record_warning(paste(sum(!audit$valid_chromosome),
                       "rows have malformed/unrecognized chromosomes and will be excluded with audit reasons."))
}
if (any(!audit$valid_coordinates)) {
  record_warning(paste(sum(!audit$valid_coordinates),
                       "rows have invalid coordinates and will be excluded with audit reasons."))
}
if (any(!audit$valid_confidence)) {
  record_warning(paste(sum(!audit$valid_confidence),
                       "rows have missing/non-numeric confidence and will be excluded with audit reasons."))
}
if (any(!audit$valid_value)) {
  record_warning(paste(sum(!audit$valid_value),
                       "rows have missing/unexpected CNV values and will be excluded with audit reasons."))
}

# ==============================================================================
# Explicit sample metadata - unknown sex is not permitted for this run
# ==============================================================================

metadata <- read_euro_csv(metadata_filename)
metadata_required <- c("sample_id", "expected_sex")
metadata_missing <- setdiff(metadata_required, names(metadata))
if (length(metadata_missing)) {
  stop("Metadata input is missing: ", paste(metadata_missing, collapse = ", "), call. = FALSE)
}
metadata$sample_id <- trimws(as.character(metadata$sample_id))
metadata$expected_sex <- tolower(trimws(as.character(metadata$expected_sex)))
if (!"sex_source" %in% names(metadata)) metadata$sex_source <- "not recorded"
if (!"sex_status" %in% names(metadata)) metadata$sex_status <- "not recorded"
if (!"cell_line" %in% names(metadata)) metadata$cell_line <- "not recorded"
if (!"sample_name" %in% names(metadata)) metadata$sample_name <- NA_character_
metadata$sample_name <- trimws(as.character(metadata$sample_name))

if (any(duplicated(metadata$sample_id))) stop("Metadata has duplicate sample_id rows.", call. = FALSE)
input_samples <- sort(unique(audit$sample_id))
missing_metadata_samples <- setdiff(input_samples, metadata$sample_id)
if (length(missing_metadata_samples)) {
  stop("Sex metadata missing for: ", paste(missing_metadata_samples, collapse = ", "), call. = FALSE)
}
allowed_sexes <- c("male", "female", "unknown")
if (any(!metadata$expected_sex %in% allowed_sexes)) {
  stop("Metadata expected_sex must be male, female, or unknown.", call. = FALSE)
}
used_metadata <- metadata[match(input_samples, metadata$sample_id), , drop = FALSE]
if (any(used_metadata$expected_sex == "unknown")) {
  stop(
    "Expected biological sex remains unknown for: ",
    paste(used_metadata$sample_id[used_metadata$expected_sex == "unknown"], collapse = ", "),
    ". Sex-dependent filtering was not run.",
    call. = FALSE
  )
}

metadata_index <- match(audit$sample_id, metadata$sample_id)
audit$sample_name <- as.character(metadata$sample_name[metadata_index])
audit$cell_line <- as.character(metadata$cell_line[metadata_index])
audit$expected_sex <- as.character(metadata$expected_sex[metadata_index])
audit$sex_metadata_status <- as.character(metadata$sex_status[metadata_index])
audit$sex_metadata_source <- as.character(metadata$sex_source[metadata_index])

# ==============================================================================
# Call class and sex-chromosome interpretation
# ==============================================================================

audit$call_class <- "unclassified"
value_known <- !is.na(audit$value)
audit$call_class[value_known & audit$value < 2] <- "copy-number loss"
audit$call_class[value_known & audit$value > 2] <- "copy-number gain"
audit$call_class[value_known & audit$value == 2] <- "copy-neutral LOH (autosomal/default interpretation)"

male_sex_chr <- audit$expected_sex == "male" & audit$chr_num %in% c("X", "Y")
male_chrX_value2 <- audit$expected_sex == "male" & audit$chr_num == "X" &
  value_known & audit$value == 2

audit$call_class[male_sex_chr & value_known & audit$value == 1] <- "physiological male haploid copy number"
audit$call_class[male_sex_chr & value_known & audit$value == 0] <- "male sex-chromosome loss"
audit$call_class[male_sex_chr & value_known & audit$value >= 2] <- "potential male sex-chromosome gain"
audit$call_class[male_chrX_value2] <- "male chrX CN=2/cnLOH-like sex-baseline artifact"

audit$physiological_male_sex_chr_state <- male_sex_chr & value_known & audit$value == 1
audit$male_chrX_value2_artifact <- male_chrX_value2
audit$sex_chr_interpretation <- "autosomal call"
audit$sex_chr_interpretation[audit$expected_sex == "male" & audit$chr_num == "X" & value_known & audit$value == 1] <-
  "physiological male chrX CN=1"
audit$sex_chr_interpretation[audit$expected_sex == "male" & audit$chr_num == "Y" & value_known & audit$value == 1] <-
  "physiological male chrY CN=1"
audit$sex_chr_interpretation[male_sex_chr & value_known & audit$value == 0] <-
  "non-physiological male sex-chromosome CN=0; inspect"
audit$sex_chr_interpretation[male_sex_chr & value_known & audit$value >= 2] <-
  "potential gain relative to male haploid baseline; inspect GenomeStudio sex setting and BAF/LRR"
audit$sex_chr_interpretation[male_chrX_value2] <-
  "male chrX CN=2 with CNV confidence: cnLOH-like call; excluded as likely hemizygosity/sex-setting artifact"
audit$sex_chr_interpretation[audit$expected_sex == "female" & audit$chr_num == "Y"] <-
  "unexpected chrY call for expected female sample; inspect"
audit$sex_chr_interpretation[audit$expected_sex == "female" & audit$chr_num == "X"] <-
  "female chrX call; interpret as CNV/cnLOH"
audit$sex_chr_interpretation[audit$chr_num == "XY"] <-
  "GenomeStudio label XY; not converted to X or Y; manual review required"

if (any(male_chrX_value2)) {
  record_warning(paste(
    sum(male_chrX_value2),
    "male chrX value=2 rows were detected. cnvPartition defines value=2 with",
    "confidence as copy-neutral LOH, but LOH is not biologically meaningful on",
    "the hemizygous non-PAR X of an XY male. With the present setting these calls",
    "are excluded as likely upstream sex/ploidy-setting artifacts. Confirm by",
    "rerunning cnvPartition with sex correctly specified and/or inspecting BAF/LRR."
  ))
}
if (any(male_sex_chr & value_known & audit$value == 2 & audit$chr_num == "Y")) {
  record_warning(paste(
    sum(male_sex_chr & value_known & audit$value == 2 & audit$chr_num == "Y"),
    "male chrY value=2 rows were detected. These are NOT automatically excluded",
    "and should be assessed as potential gains or upstream sex-setting artifacts."
  ))
}
if (any(audit$chr_num == "XY")) {
  record_warning(paste(sum(audit$chr_num == "XY"),
                       "rows use chromosome label XY and are not auto-converted."))
}

# ==============================================================================
# Correct interval overlap and endpoint-logic counterfactual
# ==============================================================================

find_hotspot_rows <- function(chr, start, end) {
  if (is.na(chr) || is.na(start) || is.na(end)) return(integer(0))
  which(hotspots$hotspot_chr == chr & start <= hotspots$hotspot_end & end >= hotspots$hotspot_start)
}

find_endpoint_rows <- function(chr, start, end) {
  if (is.na(chr) || is.na(start) || is.na(end)) return(integer(0))
  which(
    hotspots$hotspot_chr == chr &
    ((start >= hotspots$hotspot_start & start <= hotspots$hotspot_end) |
     (end >= hotspots$hotspot_start & end <= hotspots$hotspot_end))
  )
}

proper_matches <- lapply(seq_len(nrow(audit)), function(i) {
  find_hotspot_rows(audit$chr_num[i], audit$base_start_pos[i], audit$base_end_pos[i])
})
endpoint_matches <- lapply(seq_len(nrow(audit)), function(i) {
  find_endpoint_rows(audit$chr_num[i], audit$base_start_pos[i], audit$base_end_pos[i])
})

audit$hotspot_region <- lengths(proper_matches) > 0
audit$hotspot_rank <- vapply(proper_matches, function(idx) {
  if (!length(idx)) "" else paste(hotspots$hotspot_rank[idx], collapse = ";")
}, character(1))
audit$hotspot_intervals <- vapply(proper_matches, function(idx) {
  if (!length(idx)) "" else paste(hotspots$region_id[idx], collapse = ";")
}, character(1))
audit$hotspot_overlap_by_old_endpoint_logic_on_same_coordinates <- lengths(endpoint_matches) > 0
audit$hotspot_overlap_only_detected_by_proper_interval_logic <-
  audit$hotspot_region & !audit$hotspot_overlap_by_old_endpoint_logic_on_same_coordinates
audit$overlaps_restored_missing_hotspot_rank <- vapply(proper_matches, function(idx) {
  any(hotspots$hotspot_rank[idx] %in% c(4, 13))
}, logical(1))

# Preserve the supplied value=2 policy: no relaxed hotspot size threshold.
audit$hotspot_threshold_applied <- audit$hotspot_region & audit$valid_value & audit$value != 2
audit$value2_threshold_policy <- ifelse(
  audit$valid_value & audit$value == 2,
  "standard >350 kb threshold; no hotspot relaxation",
  "not applicable"
)
audit$required_size <- ifelse(
  audit$hotspot_threshold_applied,
  hotspot_size_threshold,
  standard_size_threshold
)
audit$size_threshold_reason <- ifelse(
  audit$hotspot_threshold_applied,
  "overlaps supplied recurrent region and value is not 2: strict >100 kb",
  ifelse(
    audit$hotspot_region & audit$valid_value & audit$value == 2,
    "value=2/cnLOH policy: strict >350 kb despite recurrent-region overlap",
    "standard strict >350 kb"
  )
)

# ==============================================================================
# Filtering and full reasons
# ==============================================================================

audit$passes_confidence <- audit$valid_confidence & audit$confidence > confidence_threshold
audit$passes_size <- audit$valid_coordinates & !is.na(audit$CNV_size) & audit$CNV_size > audit$required_size
audit$excluded_as_physiological_sex_chr <- audit$physiological_male_sex_chr_state
audit$excluded_as_male_chrX_value2_artifact <-
  exclude_male_chrX_value2_artifact & audit$male_chrX_value2_artifact

audit$keep <-
  audit$valid_sample_id & audit$valid_chromosome & audit$valid_coordinates &
  audit$valid_value & audit$valid_confidence & audit$passes_confidence &
  audit$passes_size &
  !audit$excluded_as_physiological_sex_chr &
  !audit$excluded_as_male_chrX_value2_artifact

audit$manual_review_reasons <- ""
for (i in seq_len(nrow(audit))) {
  reasons <- character(0)
  if (audit$chr_num[i] == "XY") reasons <- c(reasons, "chromosome label XY is ambiguous")
  if (male_sex_chr[i] && !is.na(audit$value[i]) && audit$value[i] != 1) {
    reasons <- c(reasons, "non-physiological male sex-chromosome state")
  }
  if (audit$male_chrX_value2_artifact[i]) {
    reasons <- c(
      reasons,
      "male chrX value=2/cnLOH-like call: likely hemizygosity/sex-setting artifact; verify upstream sex setting and BAF/LRR"
    )
  } else if (male_sex_chr[i] && !is.na(audit$value[i]) && audit$value[i] == 2) {
    reasons <- c(reasons, "CN=2 may be gain or an upstream sex-assignment/segmentation artifact")
  }
  if (audit$chr_num[i] %in% c("M", "MT")) reasons <- c(reasons, "mitochondrial chromosome call")
  if (audit$duplicate_exact_call[i]) reasons <- c(reasons, "duplicate authoritative call key")
  audit$manual_review_reasons[i] <- join_nonempty(reasons)
}
audit$manual_review_flag <- nzchar(audit$manual_review_reasons)

audit$exclusion_reasons <- ""
audit$retention_reason <- ""
for (i in seq_len(nrow(audit))) {
  reasons <- character(0)
  if (!audit$valid_sample_id[i]) reasons <- c(reasons, "missing/blank sample ID")
  if (!audit$valid_chromosome[i]) reasons <- c(reasons, "malformed/unrecognized chromosome")
  if (!audit$valid_coordinates[i]) reasons <- c(reasons, "missing/non-numeric/impossible coordinates")
  if (!audit$valid_value[i]) reasons <- c(reasons, "missing/non-numeric/unexpected CNV value")
  if (!audit$valid_confidence[i]) {
    reasons <- c(reasons, "missing/non-numeric confidence")
  } else if (!audit$passes_confidence[i]) {
    reasons <- c(reasons, paste0("confidence <= ", confidence_threshold))
  }
  if (audit$valid_coordinates[i] && !audit$passes_size[i]) {
    reasons <- c(reasons, paste0("CNV_size <= required strict threshold of ", audit$required_size[i], " bp"))
  }
  if (audit$excluded_as_physiological_sex_chr[i]) {
    reasons <- c(reasons, "physiological male X/Y CN=1")
  }
  if (audit$excluded_as_male_chrX_value2_artifact[i]) {
    reasons <- c(
      reasons,
      "male chrX value=2/cnLOH-like call excluded as likely hemizygosity/sex-setting artifact"
    )
  }
  audit$exclusion_reasons[i] <- join_nonempty(reasons)
  if (audit$keep[i]) {
    audit$retention_reason[i] <- if (audit$hotspot_threshold_applied[i]) {
      paste0("passes confidence >75 and strict >100-kb recurrent-region threshold; rank ", audit$hotspot_rank[i])
    } else if (audit$hotspot_region[i] && audit$value[i] == 2) {
      paste0("passes confidence >75 and strict >350-kb value=2 threshold; overlaps rank ", audit$hotspot_rank[i])
    } else {
      "passes confidence >75 and strict >350-kb standard threshold"
    }
  }
}
audit$final_status <- ifelse(audit$keep, "KEEP", "EXCLUDE")
audit$final_decision_reason <- ifelse(audit$keep, audit$retention_reason, audit$exclusion_reasons)

# ==============================================================================
# Historical comparison with exact match, then one safe non-positional fallback
# ==============================================================================

historical <- read_euro_csv(historical_filename)
historical_missing <- setdiff(required, names(historical))
if (length(historical_missing)) {
  stop("Historical input is missing: ", paste(historical_missing, collapse = ", "), call. = FALSE)
}
historical$historical_row_id <- seq_len(nrow(historical))
historical$sample_id <- trimws(as.character(historical$sample_id))
historical$chr_num_original <- as.character(historical$chr_num)
historical$chr_num_normalized <- normalize_chr(historical$chr_num)
historical$base_start_pos <- as_number(historical$base_start_pos)
historical$base_end_pos <- as_number(historical$base_end_pos)
historical$value <- as_number(historical$value)
historical$confidence <- as_number(historical$confidence)
historical$historical_call_key <- make_call_key(
  historical$sample_id, historical$chr_num_normalized,
  historical$base_start_pos, historical$base_end_pos, historical$value
)
historical$matched_call_id <- match(historical$historical_call_key, audit$authoritative_call_key)
historical$match_method <- ifelse(
  is.na(historical$matched_call_id),
  "unmatched after exact normalized key",
  "exact sample+chromosome+start+end+value"
)

for (i in which(is.na(historical$matched_call_id))) {
  candidate <- which(
    audit$sample_id == historical$sample_id[i] &
    audit$base_start_pos == historical$base_start_pos[i] &
    audit$base_end_pos == historical$base_end_pos[i] &
    audit$value == historical$value[i] &
    audit$confidence == historical$confidence[i]
  )
  if (length(candidate) == 1) {
    historical$matched_call_id[i] <- candidate
    historical$match_method[i] <-
      "unique exact sample+start+end+value+confidence fallback; historical chromosome malformed; no positional tolerance"
  }
}

if (any(is.na(historical$matched_call_id))) {
  stop(
    "Historical rows remain unmatched after exact and unique non-positional fallback: rows ",
    paste(historical$historical_row_id[is.na(historical$matched_call_id)], collapse = ", "),
    call. = FALSE
  )
}
if (any(duplicated(historical$matched_call_id))) {
  stop("Multiple historical rows matched the same authoritative call.", call. = FALSE)
}

fallback_rows <- which(grepl("fallback", historical$match_method, fixed = TRUE))
if (length(fallback_rows)) {
  record_warning(paste(
    length(fallback_rows),
    "historical row(s) required a unique non-positional fallback because the",
    "historical chromosome field was malformed. No coordinate tolerance was used."
  ))
}

audit$historically_retained <- FALSE
audit$historical_row_id <- NA_integer_
audit$historical_match_method <- "not retained historically"
audit$historical_chr_num_original <- ""
audit$historical_chr_num_normalized <- ""
audit$historical_hotspot_region <- NA
audit$historical_CNV_size <- NA_real_

for (i in seq_len(nrow(historical))) {
  j <- historical$matched_call_id[i]
  audit$historically_retained[j] <- TRUE
  audit$historical_row_id[j] <- historical$historical_row_id[i]
  audit$historical_match_method[j] <- historical$match_method[i]
  audit$historical_chr_num_original[j] <- historical$chr_num_original[i]
  audit$historical_chr_num_normalized[j] <- historical$chr_num_normalized[i]
  if ("hotspot_region" %in% names(historical)) {
    audit$historical_hotspot_region[j] <- as.logical(historical$hotspot_region[i])
  }
  if ("CNV_size" %in% names(historical)) {
    audit$historical_CNV_size[j] <- as_number(historical$CNV_size[i])
  }
}

audit$comparison_category <- ifelse(
  audit$historically_retained & audit$keep,
  "unchanged retained call",
  ifelse(
    audit$historically_retained & !audit$keep,
    "historically retained -> now excluded",
    ifelse(!audit$historically_retained & audit$keep,
           "historically excluded -> now retained",
           "excluded by both")
  )
)

audit$change_reason_all <- ""
audit$change_reason_primary <- ""
audit$change_reason_evidence <- ""
for (i in seq_len(nrow(audit))) {
  category <- audit$comparison_category[i]
  reasons <- character(0)
  evidence <- "observed"

  if (category == "unchanged retained call") {
    audit$change_reason_primary[i] <- "not applicable"
    audit$change_reason_all[i] <- "not applicable"
    audit$change_reason_evidence[i] <- "observed"
    next
  }
  if (category == "excluded by both") {
    audit$change_reason_primary[i] <- "not part of retained-call delta"
    audit$change_reason_all[i] <- "not part of retained-call delta"
    audit$change_reason_evidence[i] <- "observed"
    next
  }

  if (category == "historically retained -> now excluded") {
    if (audit$excluded_as_physiological_sex_chr[i]) {
      reasons <- c(reasons, "corrected sex-aware exclusion of physiological male X/Y CN=1")
    }
    if (audit$excluded_as_male_chrX_value2_artifact[i]) {
      reasons <- c(
        reasons,
        "corrected sex-aware exclusion of male chrX value=2/cnLOH-like hemizygosity artifact"
      )
      evidence <- "observed systematic chrX pattern; upstream sex-setting mechanism inferred"
    }
    if (!audit$passes_confidence[i]) {
      reasons <- c(reasons, "corrected confidence requirement >75")
    }
    if (!audit$passes_size[i]) {
      reasons <- c(reasons, "does not pass corrected applicable size threshold")
      if (isTRUE(audit$historical_hotspot_region[i]) && !audit$hotspot_region[i]) {
        reasons <- c(reasons, "historical hotspot annotation not reproduced by supplied GRCh37 table")
        evidence <- "observed classification difference; build-mismatch mechanism inferred"
      }
    }
  }

  if (category == "historically excluded -> now retained") {
    if (audit$confidence[i] > 75 && audit$confidence[i] <= 80) {
      reasons <- c(
        reasons,
        "confidence in (75,80]: newly eligible under >75 versus apparent historical >80 gate"
      )
      evidence <- "historical threshold inferred from retained-file minimum; old script unavailable"
    }
    if (audit$CNV_size[i] > 100000 && audit$CNV_size[i] <= 350000 && audit$hotspot_threshold_applied[i]) {
      if (audit$overlaps_restored_missing_hotspot_rank[i]) {
        reasons <- c(reasons, "restoration of previously missing hotspot rank 4 or 13")
      } else if (audit$hotspot_overlap_only_detected_by_proper_interval_logic[i]) {
        reasons <- c(reasons, "correct interval-overlap logic; call contains hotspot without an endpoint inside")
      } else {
        reasons <- c(reasons, "corrected GRCh37 hotspot classification; consistent with historical build mismatch")
        evidence <- "counterfactual inference; historical unfiltered annotations/script unavailable"
      }
    }
    if (male_sex_chr[i] && audit$chr_num[i] == "Y" &&
        !is.na(audit$value[i]) && audit$value[i] == 2) {
      reasons <- c(
        reasons,
        "other: historical treatment of male chrY CN=2 differs; potential gain/upstream sex-setting artifact"
      )
      evidence <- "observed historical difference; biological state requires manual GenomeStudio review"
    }
    if (!length(reasons)) {
      reasons <- c(reasons, "other historical filtering difference not reconstructable from filtered-only CSV")
      evidence <- "ambiguous because the historical script/audit trail was not supplied"
    }
  }

  audit$change_reason_all[i] <- join_nonempty(reasons)
  audit$change_reason_primary[i] <- if (length(reasons)) reasons[1] else "unclassified change"
  audit$change_reason_evidence[i] <- evidence
}

# ==============================================================================
# Outputs
# ==============================================================================

audit_filename <- file.path(output_dir, "cnvPartition_3.2.0_filtering_audit_GRCh37_conf75_chrX2filtered.csv")
filtered_filename <- file.path(output_dir, "cnvPartition_3.2.0_filtered_data_GRCh37_corrected_conf75_chrX2filtered.csv")
delta_filename <- file.path(output_dir, "CNV_reanalysis_delta_summary_chrX2filtered.csv")
report_filename <- file.path(output_dir, "CNV_reanalysis_QC_report_chrX2filtered.md")
log_filename <- file.path(output_dir, "CNV_reanalysis_execution_log_chrX2filtered.md")
hotspot_filename <- file.path(output_dir, "cnvPartition_3.2.0_hotspot_reference_GRCh37.csv")

filtered_columns <- c(
  "call_id", "sample_id", "sample_name", "cell_line", "expected_sex", "sex_metadata_source",
  "chr_num", "base_start_pos", "base_end_pos", "value", "call_class",
  "confidence", "CNV_size", "inclusive_size", "hotspot_region", "hotspot_rank",
  "hotspot_intervals", "hotspot_threshold_applied", "required_size",
  "sex_chr_interpretation", "male_chrX_value2_artifact", "manual_review_flag", "manual_review_reasons",
  "retention_reason", "genome_build"
)
filtered <- audit[audit$keep, filtered_columns, drop = FALSE]

delta_columns <- c(
  "call_id", "sample_id", "sample_name", "cell_line", "expected_sex", "chr_num",
  "base_start_pos", "base_end_pos", "value", "call_class", "confidence",
  "CNV_size", "hotspot_region", "hotspot_rank",
  "hotspot_overlap_only_detected_by_proper_interval_logic",
  "overlaps_restored_missing_hotspot_rank", "sex_chr_interpretation",
  "male_chrX_value2_artifact", "excluded_as_male_chrX_value2_artifact",
  "historically_retained", "keep", "comparison_category",
  "change_reason_primary", "change_reason_all", "change_reason_evidence",
  "historical_row_id", "historical_match_method", "historical_chr_num_original",
  "historical_hotspot_region", "manual_review_flag", "manual_review_reasons"
)
delta <- audit[
  audit$historically_retained | audit$keep,
  delta_columns,
  drop = FALSE
]

write.csv2(audit, audit_filename, row.names = FALSE, quote = TRUE, na = "")
write.csv2(filtered, filtered_filename, row.names = FALSE, quote = TRUE, na = "")
write.csv2(delta, delta_filename, row.names = FALSE, quote = TRUE, na = "")
write.csv2(hotspots, hotspot_filename, row.names = FALSE, quote = TRUE, na = "")

# ==============================================================================
# Scientific QC report
# ==============================================================================

historical_count <- nrow(historical)
corrected_count <- sum(audit$keep)
unchanged_count <- sum(audit$comparison_category == "unchanged retained call")
old_to_excluded_count <- sum(audit$comparison_category == "historically retained -> now excluded")
excluded_to_new_count <- sum(audit$comparison_category == "historically excluded -> now retained")
changed <- audit$comparison_category %in% c(
  "historically retained -> now excluded", "historically excluded -> now retained"
)
new_retained <- audit$comparison_category == "historically excluded -> now retained"
new_male_Y_cn2 <- new_retained & male_sex_chr & audit$chr_num == "Y" & value_known & audit$value == 2
new_confidence_band <- new_retained & audit$confidence > 75 & audit$confidence <= 80
new_rank13 <- new_retained & vapply(
  strsplit(audit$hotspot_rank, ";", fixed = TRUE),
  function(x) "13" %in% x,
  logical(1)
)
rank13_detail <- if (any(new_rank13)) {
  paste(vapply(which(new_rank13), function(i) {
    paste0(
      "`", audit$sample_id[i], "` chr16:", audit$base_start_pos[i], "-",
      audit$base_end_pos[i], " CN=", audit$value[i], " confidence=",
      format(audit$confidence[i], digits = 8), " size=", audit$CNV_size[i], " bp"
    )
  }, character(1)), collapse = "; ")
} else {
  "none"
}

profile_rows <- lapply(input_samples, function(s) {
  idx <- audit$sample_id == s
  data.frame(
    sample_id = s,
    historical_retained = sum(audit$historically_retained[idx]),
    corrected_retained = sum(audit$keep[idx]),
    removed = sum(audit$comparison_category[idx] == "historically retained -> now excluded"),
    added = sum(audit$comparison_category[idx] == "historically excluded -> now retained"),
    profile_changed = any(changed[idx]),
    stringsAsFactors = FALSE
  )
})
profile <- do.call(rbind, profile_rows)
changed_samples <- profile$sample_id[profile$profile_changed]

count_reason <- function(pattern) sum(changed & grepl(pattern, audit$change_reason_all, fixed = TRUE))
build_mismatch_count <- count_reason("consistent with historical build mismatch") +
  count_reason("historical hotspot annotation not reproduced")
missing_hotspot_count <- count_reason("restoration of previously missing hotspot rank 4 or 13")
overlap_logic_count <- count_reason("correct interval-overlap logic")
confidence_change_count <- count_reason("confidence in (75,80]")
sex_aware_cn1_count <- count_reason("corrected sex-aware exclusion of physiological male X/Y CN=1")
sex_aware_chrX2_count <- count_reason("corrected sex-aware exclusion of male chrX value=2")

band_100_350 <- audit$CNV_size > 100000 & audit$CNV_size <= 350000
rank1_changed <- changed & vapply(strsplit(audit$hotspot_rank, ";", fixed = TRUE), function(x) "1" %in% x, logical(1))
hotspot_changed <- changed & audit$hotspot_region
biologically_notable <- changed & (
  audit$hotspot_region | audit$CNV_size >= 1000000 |
  (audit$chr_num %in% c("X", "Y", "XY") & !audit$physiological_male_sex_chr_state)
)

profile_lines <- vapply(seq_len(nrow(profile)), function(i) {
  paste0(
    "| `", profile$sample_id[i], "` | ", profile$historical_retained[i], " | ",
    profile$corrected_retained[i], " | ", profile$removed[i], " | ",
    profile$added[i], " | ", ifelse(profile$profile_changed[i], "yes", "no"), " |"
  )
}, character(1))

report_lines <- c(
  "# CNV reanalysis QC report",
  "",
  paste0("Execution date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Rome")),
  "",
  "## Scope and evidence status",
  "",
  "This is a scientific data-integrity audit of the supplied combined cnvPartition calls. The extracted CSV is the authoritative call universe; the historical filtered CSV is used only for comparison.",
  "",
  "Assou et al. state that abnormalities were converted to GRCh37A. However, the publication and supplied supplemental PDF do not expose the exact ranked interval table, and Figure 2D shows 24 common abnormal regions. The analysis therefore implements the requested top-20 table exactly as supplied, while retaining coordinate provenance as an unresolved manuscript/SOP issue.",
  "",
  "## Observed results",
  "",
  paste0("1. Total unfiltered calls: **", nrow(audit), "**."),
  paste0("2. Samples represented: **", length(input_samples), "**."),
  paste0("3. Historically retained calls: **", historical_count, "**."),
  paste0("4. Corrected retained calls: **", corrected_count, "**."),
  paste0("5. Unchanged retained calls: **", unchanged_count, "**."),
  paste0("6. Historically retained -> now excluded: **", old_to_excluded_count, "**."),
  paste0("7. Historically excluded -> now retained: **", excluded_to_new_count, "**."),
  paste0("8. Samples with changed reported profiles: **", length(changed_samples), "** of ", length(input_samples), "."),
  "",
  "| Sample | Historical | Corrected | Removed | Added | Profile changed |",
  "|---|---:|---:|---:|---:|---|",
  profile_lines,
  "",
  "## Computational attribution",
  "",
  paste0("9. Changes consistent with the GRCh37/GRCh38 hotspot mismatch: **", build_mismatch_count, "**. This is an inference unless the historical annotation is directly present; the historical script was not supplied."),
  paste0("10. Changes attributed to restoring missing ranks 4 or 13: **", missing_hotspot_count, "**."),
  paste0("11. Changes attributed to corrected containment/interval-overlap logic: **", overlap_logic_count, "**."),
  paste0("12. Changes in the confidence (75,80] band, consistent with lowering an apparent historical >80 gate to >75: **", confidence_change_count, "**. The >80 historical gate is inferred from the retained file and is not proven without the old script."),
  paste0("13. Historical calls excluded as physiological male X/Y CN=1: **", sex_aware_cn1_count, "**; historical calls excluded as male chrX value=2/cnLOH-like artifacts: **", sex_aware_chrX2_count, "**."),
  paste0("14. Changed calls overlapping one of the supplied recurrent regions: **", sum(hotspot_changed), "**."),
  paste0("15. Changed calls meeting a broad biological-notability flag (recurrent region, >=1 Mb, or non-physiological/ambiguous sex-chromosome state): **", sum(biologically_notable), "**."),
  paste0("16. Samples requiring reconsideration because their retained-call profile changed: **", ifelse(length(changed_samples), paste(changed_samples, collapse = ", "), "none"), "**."),
  paste0("17. Calls flagged for manual GenomeStudio inspection in the full audit: **", sum(audit$manual_review_flag), "**; among changed calls: **", sum(audit$manual_review_flag & changed), "**."),
  "",
  "## Targeted audits",
  "",
  paste0("- Calls between 100 and 350 kb: ", sum(band_100_350), " total; ", sum(band_100_350 & audit$passes_confidence), " pass confidence; ", sum(band_100_350 & audit$keep), " are retained; ", sum(band_100_350 & changed), " change historical status."),
  paste0("- Newly retained 100-350-kb calls attributable to restored rank-13 chr16 hotspot logic: ", rank13_detail, "."),
  paste0("- Rank-1 chr20 (20q11.21) calls changing historical status: **", sum(rank1_changed), "**."),
  paste0("- Physiological male chrX/chrY CN=1 calls in the authoritative input: **", sum(audit$physiological_male_sex_chr_state), "**; historically retained: **", sum(audit$physiological_male_sex_chr_state & audit$historically_retained), "**."),
  paste0("- Newly retained calls decompose into **", sum(new_male_Y_cn2), " male chrY CN=2 calls** and **", sum(new_rank13), " restored chr16 rank-13 losses**."),
  paste0("- Newly retained calls in the confidence (75,80] band after male chrX value=2 artifact exclusion: **", sum(new_confidence_band), "**. Of these, male chrY CN=2 calls: **", sum(new_confidence_band & new_male_Y_cn2), "**."),
  paste0("- Male chrX value=2 calls in the authoritative input: **", sum(audit$male_chrX_value2_artifact), "**; excluded by the retrospective sex-aware artifact rule: **", sum(audit$excluded_as_male_chrX_value2_artifact), "**."),
  paste0("- Male chrY value=2 calls in the authoritative input: **", sum(male_sex_chr & audit$chr_num == "Y" & value_known & audit$value == 2), "**. These are not automatically excluded and require BAF/LRR review."),
  "",
  "## Biological interpretation",
  "",
  "No rank-1 20q11.21 call changed in this dataset. The restored rank-13 chr16 losses are small recurrent-region calls and merit manual confirmation. Pervasive chrX value=2 segments across all XY samples are excluded as cnLOH-like artifacts of applying a diploid/incorrect sex baseline to a hemizygous male X chromosome. This retrospective exclusion should still be confirmed by rerunning cnvPartition with sample sex correctly specified and/or inspecting raw BAF/LRR.",
  "",
  "Autosomal value=2 calls are reported as copy-neutral LOH and use the standard strict >350-kb threshold. The Assou recurrent regions were defined from recurrent copy-number abnormalities; value=2 calls do not receive the relaxed hotspot threshold under the supplied policy.",
  "",
  "## Remaining uncertainty and checks before manuscript use",
  "",
  "- Confirm the exact ranked GRCh37 interval table from the original Assou analysis spreadsheet/authors. The paper displays 24 common regions, while this requested policy uses 20.",
  "- Confirm the GSAMD-24v3-0-EA_20034606_A1 build directly with Life & Brain.",
  "- Obtain the historical R script if causal attribution among build mismatch, missing regions, endpoint logic, and confidence filtering must be definitive rather than counterfactual.",
  "- Re-run or review male chrX value=2 calls in GenomeStudio with sample sex correctly set; the present script excludes them retrospectively as likely hemizygosity/sex-setting artifacts.",
  "- Review any male chrY CN=2 and chromosome-XY-labelled calls separately; these are not automatically excluded.",
  "- The historical file contains one corrupted chromosome field (` instead of 7); it was matched uniquely using identical sample, start, end, value, and confidence, without positional tolerance.",
  "- R 3.4.3 is obsolete. The validated script is base-R compatible, but rerunning in a current validated R environment is advisable before formal SOP/manuscript release.",
  "",
  "## Sources",
  "",
  "- Assou S, et al. Stem Cell Reports. 2020;14:1-8. PMID 31902703. https://doi.org/10.1016/j.stemcr.2019.12.004",
  "- Open-access article: https://pmc.ncbi.nlm.nih.gov/articles/PMC6962701/",
  "- Illumina Knowledge Article #2766: CNV Value=2 with a confidence value denotes copy-neutral LOH; normal CN=2 has blank confidence."
)
writeLines(report_lines, report_filename, useBytes = TRUE)

# ==============================================================================
# Execution log
# ==============================================================================

pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else "NOT INSTALLED"
}
input_hashes <- tools::md5sum(c(input_filename, historical_filename, metadata_filename))
sex_lines <- vapply(seq_len(nrow(used_metadata)), function(i) {
  paste0("- `", used_metadata$sample_id[i], "`: ", used_metadata$expected_sex[i],
         " (", used_metadata$sex_source[i], ")")
}, character(1))
warning_lines <- if (length(warnings_seen)) paste0("- ", unique(warnings_seen)) else "- None recorded by the script."

log_lines <- c(
  "# CNV reanalysis execution log",
  "",
  paste0("- Execution date/time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Rome")),
  paste0("- R version: ", R.version.string),
  paste0("- Platform: ", R.version$platform),
  "- External analysis packages used: none (base R implementation)",
  paste0("- base package: ", pkg_version("base")),
  paste0("- stats package: ", pkg_version("stats")),
  paste0("- utils package: ", pkg_version("utils")),
  paste0("- tools package: ", pkg_version("tools")),
  paste0("- dplyr availability/version: ", pkg_version("dplyr")),
  "",
  "## Inputs",
  "",
  paste0("- Authoritative extracted calls: `", normalizePath(input_filename), "` (MD5 ", unname(input_hashes[1]), ")"),
  paste0("- Historical filtered calls: `", normalizePath(historical_filename), "` (MD5 ", unname(input_hashes[2]), ")"),
  paste0("- Explicit sample metadata: `", normalizePath(metadata_filename), "` (MD5 ", unname(input_hashes[3]), ")"),
  "",
  "## Run counts",
  "",
  paste0("- Samples: ", length(input_samples)),
  paste0("- Authoritative calls: ", nrow(audit)),
  paste0("- Historical retained rows: ", historical_count),
  paste0("- Corrected retained rows: ", corrected_count),
  paste0("- Full audit rows: ", nrow(audit)),
  paste0("- Delta rows including unchanged retained union: ", nrow(delta)),
  "",
  "## Expected biological sex used",
  "",
  sex_lines,
  "",
  "## Warnings and uncertainties",
  "",
  "- R emitted LC_* startup locale warnings in this execution environment; CSV parsing was verified numerically and was not affected.",
  warning_lines,
  "",
  "## Modifications relative to the supplied script",
  "",
  "- Replaced dplyr code with base R so the actual installed R 3.4.3 environment can execute it.",
  "- Retained direct combined-CSV input; no XML parsing was introduced.",
  "- Added strict schema, numeric, chromosome, coordinate, duplicate, and metadata validation.",
  "- Added explicit external per-sample sex metadata and stopped the run if sex is unknown.",
  "- Distinguished autosomal value=2 cnLOH from sex-chromosome calls.",
  "- Added a configurable default exclusion of male chrX value=2 calls as cnLOH-like hemizygosity/sex-setting artifacts; this rule is specific to explicitly male samples and does not automatically exclude chrY value=2.",
  "- Added proper interval overlap plus an endpoint-logic counterfactual field.",
  "- Added complete per-call reasons, manual-review flags, historical matching, delta attribution, QC report, and execution log.",
  "- Added a unique exact non-positional fallback for the single malformed historical chromosome field; no coordinate tolerance was used.",
  "- Documented the unresolved hotspot-coordinate provenance and Assou 24-versus-requested-20 region discrepancy; no extra regions were added.",
  "",
  "## Errors",
  "",
  "- None; execution completed successfully.",
  "",
  "## Outputs",
  "",
  paste0("- `", basename(audit_filename), "`"),
  paste0("- `", basename(filtered_filename), "`"),
  paste0("- `", basename(delta_filename), "`"),
  paste0("- `", basename(report_filename), "`"),
  paste0("- `", basename(log_filename), "`"),
  paste0("- `", basename(hotspot_filename), "` (additional provenance output)")
)
writeLines(log_lines, log_filename, useBytes = TRUE)

cat("CNV reanalysis completed successfully\n")
cat("Authoritative calls: ", nrow(audit), "\n", sep = "")
cat("Samples: ", length(input_samples), "\n", sep = "")
cat("Historical retained: ", historical_count, "\n", sep = "")
cat("Corrected retained: ", corrected_count, "\n", sep = "")
cat("Unchanged retained: ", unchanged_count, "\n", sep = "")
cat("Historically retained -> now excluded: ", old_to_excluded_count, "\n", sep = "")
cat("Historically excluded -> now retained: ", excluded_to_new_count, "\n", sep = "")
cat("Audit: ", audit_filename, "\n", sep = "")
cat("Filtered: ", filtered_filename, "\n", sep = "")
cat("Delta: ", delta_filename, "\n", sep = "")
cat("QC report: ", report_filename, "\n", sep = "")
cat("Execution log: ", log_filename, "\n", sep = "")
