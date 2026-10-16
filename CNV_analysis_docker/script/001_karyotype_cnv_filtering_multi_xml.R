#!/usr/bin/env Rscript
# Multi-XML companion to corrected 001_karyotpe_cnv_filtering.R.
# Original interface preserved: --xml=FILE1,FILE2|list.txt
# --comment_threshold=75 --hotspot_size_threshold=100000
# --non_hotspot_size_threshold=350000 --output_dir=. --combined_name=combined
# Optional: --missing-confidence=error|exclude
# Reuses the corrected single-XML parser/filter; no duplicated regex or hotspot table.
args <- commandArgs(trailingOnly=TRUE)
allowed <- c('xml','comment_threshold','hotspot_size_threshold',
             'non_hotspot_size_threshold','output_dir','combined_name','missing-confidence')
flags <- list()
for (arg in args) {
  eq <- regexpr('=',arg,fixed=TRUE)[1]
  if (!grepl('^--',arg) || eq<4L) stop('Expected --flag=value: ',arg,call.=FALSE)
  key <- substring(arg,3L,eq-1L); val <- substring(arg,eq+1L)
  if (!key %in% allowed || !is.null(flags[[key]]) || !nzchar(trimws(val)))
    stop('Unknown/repeated/empty option: ',key,call.=FALSE)
  flags[[key]] <- val
}
if (is.null(flags$xml)) stop('Required: --xml=FILE1.xml,FILE2.xml or --xml=list.txt',call.=FALSE)
x <- flags$xml
if (file.exists(x)) {
  files <- if (grepl('\\.xml$',x,ignore.case=TRUE)) x else trimws(readLines(x,warn=FALSE))
} else files <- trimws(strsplit(x,',',fixed=TRUE)[[1]])
files <- unique(files[nzchar(files)])
if (!length(files) || any(!file.exists(files))) stop('Missing XML input(s): ',paste(files[!file.exists(files)],collapse=', '),call.=FALSE)
if (anyDuplicated(basename(files))) stop('Different inputs with identical basenames: source_xml would be ambiguous.',call.=FALSE)
get <- function(key,default) if(is.null(flags[[key]])) default else flags[[key]]
out <- get('output_dir',getwd()); combined <- get('combined_name','combined')
if (basename(combined)!=combined || combined %in% c('.','..')) stop('combined_name must be a filename prefix, not a path.',call.=FALSE)
self <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
if (length(self)!=1L) stop('Run this wrapper with Rscript.',call.=FALSE)
parser_file <- file.path(dirname(normalizePath(sub('^--file=','',self))), '001_karyotpe_cnv_filtering.R')
parser <- new.env(parent=globalenv()); sys.source(parser_file,envir=parser)
if (!is.function(parser$parse_confidence) || !is.function(parser$main)) stop('Companion 001 is not the corrected generalized version.',call.=FALSE)
# Each input is validated in temporary staging; combined outputs are written only
# after every file has passed. Candidate filtering is delegated to corrected001; final sex policy stays in003.
staging <- tempfile('cnv-multi-');dir.create(staging)
on.exit_cleanup <- function() unlink(staging,recursive=TRUE)
results <- tryCatch(lapply(files,function(file) {
  parser$main(c(file,get('comment_threshold','75'),get('hotspot_size_threshold','100000'),
    get('non_hotspot_size_threshold','350000'),staging,
    paste0('--missing-confidence=',get('missing-confidence','error'))))
  prefix <- tools::file_path_sans_ext(basename(file))
  read <- function(suffix) read.csv2(file.path(staging,paste0(prefix,suffix)),stringsAsFactors=FALSE,check.names=FALSE)
  add_source <- function(d) { d$source_xml <- rep(basename(file),nrow(d)); d }
  list(extracted=add_source(read('_extracted_data.csv')),
       filtered=add_source(read('_filtered_data.csv')),
       audit=add_source(read('_extraction_audit.csv')))
}),finally=on.exit_cleanup())
extracted <- do.call(rbind,lapply(results,function(x)x$extracted))
filtered <- do.call(rbind,lapply(results,function(x)x$filtered))
audit <- do.call(rbind,lapply(results,function(x)x$audit))
key <- do.call(paste,c(extracted[c('sample_id','chr_num','base_start_pos','base_end_pos','value')],sep='|'))
audit$duplicate_call_key <- duplicated(key) | duplicated(key,fromLast=TRUE)
if (any(audit$duplicate_call_key)) warning(sum(audit$duplicate_call_key),' rows have duplicate call keys; all copies retained, see combined audit.',call.=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write <- function(d,suffix) write.csv2(d,file.path(out,paste0(combined,suffix)),row.names=FALSE,quote=TRUE,na='')
write(extracted,'_extracted_data.csv');write(filtered,'_filtered_data.csv');write(audit,'_extraction_audit.csv')
summary <- data.frame(source_xml=basename(files),
  extracted_calls=vapply(results,function(x)nrow(x$extracted),integer(1)),
  missing_confidence=vapply(results,function(x)sum(is.na(x$extracted$confidence)),integer(1)),
  step001_retained=vapply(results,function(x)nrow(x$filtered),integer(1)))
write(summary,'_extraction_summary.csv')
cat('Corrected multi-XML script:',normalizePath(sub('^--file=','',self)),'\n')
cat('Corrected parser:',parser_file,'\n');print(summary)
cat('Combined rows:',nrow(extracted),'| Missing confidence:',sum(is.na(extracted$confidence)),'\n')
