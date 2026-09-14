#!/usr/bin/env Rscript

# ==============================================================================
# Build CNV_sample_metadata.csv from one or more Illumina SampleSheets
# ==============================================================================
#
# Purpose
#   Produces the third input required by
#   karyotype_cnv_filtering_GRCh37_conf75_validated_chrX2filtered.R:
#   a sample_id -> expected_sex lookup table (that script hard-stops if any
#   sample it sees is missing from this table, or marked "unknown").
#
# Why two candidate rows per physical sample
#   GenomeStudio bookmark exports label each sample either as
#     "<Sample_ID> [<n>]"                              (when Sample_ID was set), or
#     "<SentrixBarcode_A>_<SentrixPosition_A> [<n>]"    (fallback, when it was not)
#   and there is no way to know in advance, from the SampleSheet alone, which
#   convention GenomeStudio actually used for a given sample -- your own
#   hand-built CNV_sample_metadata.csv mixes both within the same project.
#   This script therefore emits BOTH candidate sample_id strings per physical
#   sample. Only the one that actually matches a bookmark's sample_id is ever
#   used by the downstream filtering script; the unused row is simply ignored.
#
# Provenance of the "[<n>]" suffix  -- READ BEFORE TRUSTING ON A NEW BATCH
#   Reverse-engineered by comparing your hand-built CNV_sample_metadata.csv
#   against samplesheet_2024287ILL_GSAUTU_N8.csv (chip 208602150014, 8
#   samples): "<n>" equals the alphabetic ROW of Sample_Well (A=1, B=2, ...
#   H=8), NOT the row of SentrixPosition_A -- confirmed exactly on all 4
#   overlapping samples (IT97...D100022 well E12 -> [5], IT70...D100023 well
#   F12 -> [6], IT43...D100024 well G12 -> [7], IT16...D100025 well H12 ->
#   [8]). This is an empirical inference from one chip, not a documented
#   GenomeStudio rule. Spot-check it against a couple more bookmark exports
#   (chip 210276010024 in particular) before trusting it project-wide.
#
# Cell line -- single value, OR one value per sample
#   The cell line is NOT read from the samplesheet (it isn't a samplesheet
#   field) and is deliberately NOT hardcoded, since a silently-baked default
#   would happily mislabel the next batch. --cell_line accepts two forms:
#
#     1) A single bare value -> applied to every sample in this run:
#          --cell_line="WTC11 hiPSC"
#
#     2) A comma-separated list of Sample_ID=CellLine pairs (Sample_ID must
#        match the samplesheet's Sample_ID column exactly, case-sensitive) ->
#        different lines per sample. Add a "default=" entry to cover every
#        sample NOT explicitly listed; without it, unlisted samples get
#        cell_line = "not recorded" and a warning naming them:
#          --cell_line="IT16GSAUTUD100025=WTC11 hiPSC,IT97GSAUTUD100022=HUES9,default=WTC11 hiPSC"
#
#   Whichever form (and whichever entry within it) produced each row's value
#   is recorded per row in cell_line_source, the same way sex_source records
#   how expected_sex was derived. If --cell_line is omitted entirely, the
#   script prompts once for a single value applied to all samples; for
#   per-sample assignment you must pass --cell_line explicitly.
#
# Restricting to a subset of samples (optional)
#   --samples restricts the run to only the listed Sample_ID(s), instead of
#   every sample found in the samplesheet(s). Two forms:
#     --samples="IT16GSAUTUD100025,IT97GSAUTUD100022"      (inline list), or
#     --samples=my_sample_list.txt                          (one Sample_ID per line)
#   Requested IDs not found in any samplesheet are warned about, not fatal.
#   Omit --samples to process every sample in the samplesheet(s) (default).
#
# Resolving the Sample_ID-vs-Barcode_Position ambiguity against the real XML
# (optional, but the only way to know for sure)
#   This script cannot tell, from the SampleSheet alone, which of the two
#   candidate sample_id forms GenomeStudio actually used for a given sample --
#   that is decided when the bookmark XML is exported, not by this script.
#   karyotpe_cnv_filtering.R (script 1) reads that XML and writes out the real
#   sample_id values it found, in its "<xmlname>_extracted_data.csv" output.
#   Pass that file here with --extracted=... and this script will, per
#   physical sample, keep only whichever of the two candidate forms actually
#   appears in it (instead of emitting both), and warn about any sample for
#   which NEITHER form appears (kept as both, flagged for manual review) --
#   this is the definitive check, since it looks at the real bookmark data
#   rather than guessing from the SampleSheet:
#     --extracted="KS_CNV01_..._extracted_data.csv"
#   Without --extracted, both candidate forms are always kept per sample (the
#   safe default from before) and sample_id_verification says so.
#
# Usage
#   Rscript build_CNV_sample_metadata.R samplesheet1.csv [samplesheet2.csv ...] output.csv \
#     [--cell_line="WTC11 hiPSC"] [--samples="ID1,ID2,..."] [--extracted="..._extracted_data.csv"]
#
# Input samplesheet format
#   Standard Illumina SampleSheet: free-form header lines, a "[Data]" marker,
#   then a header row and one row per sample. Required columns in the [Data]
#   block: Sample_ID, Sample_Well, SentrixBarcode_A, SentrixPosition_A, Gender.
#   Sample_Name is optional but used when present (see Output format).
#
# Output format
#   Semicolon-separated, comma-decimal ("European" CSV via write.csv2), quoted
#   -- the format karyotype_cnv_filtering_GRCh37_conf75_validated_chrX2filtered.R
#   reads with read.csv2(). Columns: sample_id, sample_name, expected_sex,
#   sex_source, cell_line, cell_line_source.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

# Pull one "--flag=value" out of an args vector; returns the value (or NA if
# absent) and the args vector with that flag removed.
extract_flag <- function(args, flag_prefix) {
  pattern <- paste0("^", flag_prefix)
  hit <- grep(pattern, args, value = TRUE)
  value <- if (length(hit)) sub(pattern, "", hit[1]) else NA_character_
  list(value = value, remaining = args[!grepl(pattern, args)])
}

step_cell_line <- extract_flag(args, "--cell_line=")
cell_line_arg  <- step_cell_line$value
step_samples   <- extract_flag(step_cell_line$remaining, "--samples=")
samples_arg    <- step_samples$value
step_extracted <- extract_flag(step_samples$remaining, "--extracted=")
extracted_arg  <- step_extracted$value
positional_args <- step_extracted$remaining

if (length(positional_args) < 2) {
  stop(
    "Usage: Rscript build_CNV_sample_metadata.R samplesheet1.csv [samplesheet2.csv ...] output.csv ",
    "[--cell_line=\"...\"] [--samples=\"...\"] [--extracted=\"..._extracted_data.csv\"]\n",
    "  If --cell_line is omitted, you will be prompted for a single value applied to all samples.",
    call. = FALSE
  )
}
output_filename <- positional_args[length(positional_args)]
samplesheet_files <- positional_args[-length(positional_args)]

# ------------------------------------------------------------------------------
# Resolve --cell_line into either a single value or a per-sample map
# ------------------------------------------------------------------------------

parse_cell_line_arg <- function(x) {
  if (!grepl("=", x, fixed = TRUE)) {
    return(list(mode = "single", value = trimws(x)))
  }
  pairs <- strsplit(x, ",", fixed = TRUE)[[1]]
  map <- list()
  for (p in pairs) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2 || !nzchar(trimws(kv[1]))) {
      stop("Malformed --cell_line entry (expected Sample_ID=CellLine): '", trimws(p), "'", call. = FALSE)
    }
    key <- trimws(kv[1])
    val <- trimws(paste(kv[-1], collapse = "="))
    map[[key]] <- val
  }
  list(mode = "map", value = map)
}

resolve_cell_line <- function(sample_id, parsed, source_method) {
  if (parsed$mode == "single") {
    return(list(value = parsed$value, source = source_method))
  }
  m <- parsed$value
  if (!is.null(m[[sample_id]])) {
    return(list(
      value = m[[sample_id]],
      source = paste0(source_method, "; per-sample entry for '", sample_id, "'")
    ))
  }
  default_key <- names(m)[tolower(names(m)) == "default"]
  if (length(default_key)) {
    return(list(
      value = m[[default_key[1]]],
      source = paste0(source_method, "; default= fallback (no per-sample entry for '", sample_id, "')")
    ))
  }
  list(
    value = NA_character_,
    source = paste0(source_method, "; NO per-sample entry for '", sample_id, "' and no default= fallback given")
  )
}

if (is.na(cell_line_arg) || !nzchar(trimws(cell_line_arg))) {
  cat("Enter the cell line to apply to ALL samples in this run (e.g. 'WTC11 hiPSC').\n")
  cat("(For different cell lines per sample, re-run with e.g.\n")
  cat('  --cell_line="ID1=Line1,ID2=Line2,default=Line3" instead of answering here.)\n')
  cat("Cell line: ")
  cell_line_arg <- readLines(con = "stdin", n = 1)
  cell_line_source_method <- "entered interactively when the script was run"
} else {
  cell_line_source_method <- "passed via --cell_line on the command line"
}
if (!nzchar(trimws(cell_line_arg))) {
  stop("No cell line was provided (empty input) -- re-run and supply one.", call. = FALSE)
}
cell_line_source_method <- paste0(cell_line_source_method, " (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")")
cell_line_parsed <- parse_cell_line_arg(trimws(cell_line_arg))

default_sex_note <- "derived from samplesheet Gender column"

# ------------------------------------------------------------------------------
# Read the [Data] block of one Illumina SampleSheet
# ------------------------------------------------------------------------------

read_samplesheet_data_block <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  data_marker <- grep("^\\[Data\\]", raw_lines)
  if (length(data_marker) != 1) {
    stop("Could not find a single '[Data]' section in ", path, call. = FALSE)
  }

  header <- trimws(strsplit(raw_lines[data_marker + 1], ";", fixed = TRUE)[[1]])

  data_lines <- raw_lines[(data_marker + 2):length(raw_lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (!length(data_lines)) stop("No sample rows found after [Data] in ", path, call. = FALSE)

  split_rows <- strsplit(data_lines, ";", fixed = TRUE)
  # pad/truncate every row to header length so a ragged trailing ';' doesn't break things
  split_rows <- lapply(split_rows, function(r) {
    length(r) <- length(header)
    r
  })
  mat <- do.call(rbind, split_rows)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- header
  df$source_samplesheet <- basename(path)
  df
}

required_cols <- c("Sample_ID", "Sample_Well", "SentrixBarcode_A", "SentrixPosition_A", "Gender")

all_rows <- do.call(rbind, lapply(samplesheet_files, function(f) {
  df <- read_samplesheet_data_block(f)
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) {
    stop("Samplesheet ", f, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!"Sample_Name" %in% names(df)) df$Sample_Name <- NA_character_
  df
}))

all_rows$Sample_ID <- trimws(all_rows$Sample_ID)
all_rows$Sample_Name <- trimws(all_rows$Sample_Name)
all_rows$Sample_Well <- trimws(all_rows$Sample_Well)
all_rows$SentrixBarcode_A <- trimws(all_rows$SentrixBarcode_A)
all_rows$SentrixPosition_A <- trimws(all_rows$SentrixPosition_A)

if (any(!nzchar(all_rows$Sample_ID))) {
  stop("Blank Sample_ID in one or more samplesheet rows -- fix the samplesheet before proceeding.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Optional: restrict to a subset of samples via --samples
# ------------------------------------------------------------------------------

if (!is.na(samples_arg) && nzchar(trimws(samples_arg))) {
  if (file.exists(samples_arg)) {
    requested_ids <- trimws(readLines(samples_arg, warn = FALSE))
  } else {
    requested_ids <- trimws(strsplit(samples_arg, ",", fixed = TRUE)[[1]])
  }
  requested_ids <- unique(requested_ids[nzchar(requested_ids)])
  if (!length(requested_ids)) stop("--samples was given but contained no usable Sample_ID.", call. = FALSE)

  not_found <- setdiff(requested_ids, all_rows$Sample_ID)
  if (length(not_found)) {
    warning(
      "--samples requested Sample_ID(s) not found in any samplesheet (ignored): ",
      paste(not_found, collapse = ", "),
      call. = FALSE
    )
  }

  before_n <- nrow(all_rows)
  all_rows <- all_rows[all_rows$Sample_ID %in% requested_ids, , drop = FALSE]
  if (!nrow(all_rows)) {
    stop("--samples matched none of the samples found in the samplesheet(s).", call. = FALSE)
  }
  cat("--samples filter applied: ", nrow(all_rows), " of ", before_n, " samplesheet rows kept.\n", sep = "")
}

# ------------------------------------------------------------------------------
# Derive the GenomeStudio "[n]" index, expected_sex, and per-row cell_line
# ------------------------------------------------------------------------------

well_row_index <- function(well) {
  ifelse(nzchar(well), match(toupper(substr(well, 1, 1)), LETTERS), NA_integer_)
}

sex_map <- function(gender) {
  g <- tolower(trimws(gender))
  ifelse(g == "male", "male", ifelse(g == "female", "female", "unknown"))
}

idx <- well_row_index(all_rows$Sample_Well)
expected_sex <- sex_map(all_rows$Gender)
barcode_pos <- paste0(all_rows$SentrixBarcode_A, "_", all_rows$SentrixPosition_A)

unrecognized_gender <- unique(all_rows$Gender[expected_sex == "unknown" & nzchar(trimws(all_rows$Gender))])
if (length(unrecognized_gender)) {
  warning(
    "Unrecognized Gender value(s), mapped to 'unknown' (the downstream filtering ",
    "script stops on any 'unknown' expected_sex): ", paste(unrecognized_gender, collapse = ", "),
    call. = FALSE
  )
}
if (any(is.na(idx))) {
  warning(
    sum(is.na(idx)), " row(s) have an unreadable Sample_Well; no '[n]' suffix could be built for ",
    "them, so only the bare (no-suffix) sample_id form is emitted for those rows.",
    call. = FALSE
  )
}

suffix <- ifelse(is.na(idx), "", paste0(" [", idx, "]"))

cell_line_resolved <- lapply(all_rows$Sample_ID, resolve_cell_line,
                              parsed = cell_line_parsed, source_method = cell_line_source_method)
cell_line_values <- vapply(cell_line_resolved, function(x) {
  if (is.na(x$value)) "not recorded" else x$value
}, character(1))
cell_line_sources <- vapply(cell_line_resolved, function(x) x$source, character(1))

missing_cell_line <- unique(all_rows$Sample_ID[vapply(cell_line_resolved, function(x) is.na(x$value), logical(1))])
if (length(missing_cell_line)) {
  warning(
    "No --cell_line entry (and no default= fallback) for: ", paste(missing_cell_line, collapse = ", "),
    " -- cell_line set to 'not recorded' for these. Add a per-sample entry or a default= fallback ",
    "if this is not intended.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Optional: resolve Sample_ID-vs-Barcode_Position against the real bookmark
# sample_id values, read from script 1's "..._extracted_data.csv" output
# ------------------------------------------------------------------------------

form_a_id <- paste0(all_rows$Sample_ID, suffix)
form_b_id <- paste0(barcode_pos, suffix)

has_extracted <- !is.na(extracted_arg) && nzchar(trimws(extracted_arg))
if (has_extracted) {
  if (!file.exists(extracted_arg)) stop("--extracted file not found: ", extracted_arg, call. = FALSE)
  extracted <- read.csv2(extracted_arg, header = TRUE, stringsAsFactors = FALSE,
                          check.names = FALSE, na.strings = c("", "NA"))
  if (!"sample_id" %in% names(extracted)) {
    stop(
      "--extracted file has no 'sample_id' column -- is this really script 1's ",
      "*_extracted_data.csv (and was it written with write.csv2, not write.csv)?",
      call. = FALSE
    )
  }
  extracted_ids <- unique(trimws(as.character(extracted$sample_id)))

  match_a <- form_a_id %in% extracted_ids
  match_b <- form_b_id %in% extracted_ids
  resolved <- match_a | match_b

  unresolved <- unique(all_rows$Sample_ID[!resolved])
  if (length(unresolved)) {
    warning(
      "Neither candidate sample_id form was found in --extracted for: ",
      paste(unresolved, collapse = ", "),
      " -- BOTH forms kept for these (unresolved); check them by hand against the real bookmark XML.",
      call. = FALSE
    )
  }
  both_matched <- unique(all_rows$Sample_ID[match_a & match_b])
  if (length(both_matched)) {
    warning(
      "Both candidate sample_id forms were found in --extracted for: ",
      paste(both_matched, collapse = ", "),
      " -- unexpected (usually means the same well/position string coincides with a ",
      "Sample_ID elsewhere); both forms kept.",
      call. = FALSE
    )
  }
  unmatched_extracted <- setdiff(extracted_ids, c(form_a_id[match_a], form_b_id[match_b]))
  if (length(unmatched_extracted)) {
    cat(
      length(unmatched_extracted),
      " sample_id value(s) in --extracted did not match any candidate from this samplesheet ",
      "(likely a different chip/batch; not an error): ",
      paste(head(unmatched_extracted, 10), collapse = ", "),
      if (length(unmatched_extracted) > 10) ", ..." else "",
      "\n", sep = ""
    )
  }

  keep_a <- match_a | !resolved
  keep_b <- match_b | !resolved
  verification_a <- ifelse(match_a, "verified against --extracted",
                            ifelse(!resolved, "UNRESOLVED - not found in --extracted; kept as fallback", "not applicable"))
  verification_b <- ifelse(match_b, "verified against --extracted",
                            ifelse(!resolved, "UNRESOLVED - not found in --extracted; kept as fallback", "not applicable"))
} else {
  keep_a <- rep(TRUE, nrow(all_rows))
  keep_b <- rep(TRUE, nrow(all_rows))
  verification_a <- rep("not checked (no --extracted given); both forms kept", nrow(all_rows))
  verification_b <- verification_a
}

# ------------------------------------------------------------------------------
# Build candidate sample_id form(s) -- both by default, narrowed to the
# verified one per sample when --extracted resolved it
# ------------------------------------------------------------------------------

form_a_df <- data.frame(
  sample_id             = form_a_id,
  sample_name           = all_rows$Sample_Name,
  expected_sex          = expected_sex,
  sex_source            = paste0(default_sex_note, " (", all_rows$source_samplesheet, "); Sample_ID form"),
  cell_line             = cell_line_values,
  cell_line_source      = cell_line_sources,
  sample_id_verification = verification_a
)[keep_a, , drop = FALSE]

form_b_df <- data.frame(
  sample_id             = form_b_id,
  sample_name           = all_rows$Sample_Name,
  expected_sex          = expected_sex,
  sex_source            = paste0(default_sex_note, " (", all_rows$source_samplesheet, "); Barcode_Position form"),
  cell_line             = cell_line_values,
  cell_line_source      = cell_line_sources,
  sample_id_verification = verification_b
)[keep_b, , drop = FALSE]

metadata <- rbind(form_a_df, form_b_df)

# Guard against two physical samples accidentally colliding on the same
# sample_id string with conflicting sex calls (should not happen, but silently
# dropping one would be worse than stopping).
dup_ids <- unique(metadata$sample_id[duplicated(metadata$sample_id)])
if (length(dup_ids)) {
  conflict_ids <- dup_ids[vapply(dup_ids, function(id) {
    length(unique(metadata$expected_sex[metadata$sample_id == id])) > 1
  }, logical(1))]
  if (length(conflict_ids)) {
    stop(
      "sample_id collision with conflicting expected_sex for: ",
      paste(conflict_ids, collapse = ", "),
      " -- resolve in the samplesheet(s) before proceeding.",
      call. = FALSE
    )
  }
}
metadata <- metadata[!duplicated(metadata$sample_id), ]
metadata <- metadata[order(metadata$sample_id), ]

write.csv2(metadata, output_filename, row.names = FALSE, quote = TRUE, na = "")

cat("Samplesheets read: ", length(samplesheet_files), "\n", sep = "")
cat("Physical samples: ", nrow(all_rows), "\n", sep = "")
cat("Cell line mode: ", cell_line_parsed$mode,
    if (cell_line_parsed$mode == "single") paste0(" ('", cell_line_parsed$value, "')") else
      paste0(" (", length(cell_line_parsed$value), " entries)"),
    "\n", sep = "")
if (has_extracted) {
  cat("sample_id verified against --extracted for ", sum(resolved), " of ", nrow(all_rows),
      " samples (", length(unresolved), " unresolved, both forms kept for those).\n", sep = "")
} else {
  cat("sample_id NOT verified against a bookmark XML/extracted file -- both candidate forms kept for every sample.\n")
}
cat("Metadata rows written: ", nrow(metadata), "\n", sep = "")
cat("Output: ", output_filename, "\n", sep = "")
