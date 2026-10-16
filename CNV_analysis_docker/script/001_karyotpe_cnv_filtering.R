#!/usr/bin/env Rscript
# Generalized GenomeStudio bookmark extraction; GRCh37 recurrent regions aligned with script003.
# Usage (same positional arguments as the existing Rmd):
# Rscript 001_karyotpe_cnv_filtering.R XML [confidence=75] [hotspot_bp=100000]
#   [standard_bp=350000] [output_dir=.] [--missing-confidence=error|exclude]
# Requires xml2. All outputs use write.csv2 (semicolon fields, comma decimals).
# Missing confidence: default error; explicit 'exclude' keeps it in extracted
# calls/audit as NA and excludes it from the INTERMEDIATE filtered CSV.
# Uses Assou et al. Table 1 (2020), all 20 GRCh37 regions and inclusive interval
# overlap. No hotspot size relaxation for value=2. This remains an intermediate
# candidate filter: ALL sex-chromosome decisions and the final output stay in
# script003, using sample sex from script002. No chromosome is excluded here.
# Coordinates from the XML are copied unchanged; no liftover is performed.

fail <- function(...) stop(..., call.=FALSE)
read_args <- function(args) {
  flag <- grep('^--missing-confidence=', args, value=TRUE)
  if (length(flag)>1L) fail('Repeated --missing-confidence flag.')
  mode <- if (length(flag)) sub('^--missing-confidence=', '', flag) else 'error'
  if (!mode %in% c('error','exclude')) fail('--missing-confidence must be error or exclude.')
  pos <- args[!grepl('^--missing-confidence=',args)]
  if (any(grepl('^--',pos)) || !length(pos) || length(pos)>5L)
    fail('Usage: Rscript 001_karyotpe_cnv_filtering.R XML [confidence] [hotspot_bp] [standard_bp] [output_dir] [--missing-confidence=error|exclude]')
  numeric_arg <- function(i, default) {
    raw <- if (length(pos)>=i) pos[i] else as.character(default)
    # system2() converts numeric 100000 to '1e+05'; accept R's scientific notation.
    threshold_pattern <- '^[+-]?(?:[0-9]+(?:[.,][0-9]+)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?$'
    if (!grepl(threshold_pattern,raw,perl=TRUE)) fail('Invalid numeric threshold: ',raw)
    value <- as.numeric(chartr(',','.',raw))
    if (!is.finite(value) || value<0) fail('Invalid threshold: ',raw)
    value
  }
  list(xml=pos[1], confidence=numeric_arg(2,75), hotspot=numeric_arg(3,100000),
       standard=numeric_arg(4,350000), out=if(length(pos)>=5L) pos[5] else getwd(), missing=mode)
}

parse_confidence <- function(text, context, missing='error') {
  label <- 'CNV[[:space:]]+Confidence[[:space:]]*:'
  hits <- gregexpr(label,text,perl=TRUE,ignore.case=TRUE)[[1]]
  n <- if (hits[1]<0L) 0L else length(hits)
  if (n>1L) fail('Multiple CNV Confidence labels: ',context)
  if (!n) {
    if (missing=='error') fail('Missing CNV Confidence: ',context)
    return(list(value=NA_real_, token='', status='missing_label'))
  }
  rest <- trimws(substring(text,hits[1]+attr(hits,'match.length')[1]))
  if (!nzchar(rest)) {
    if (missing=='error') fail('Empty CNV Confidence: ',context)
    return(list(value=NA_real_,token='',status='empty_value'))
  }
  # Parse a COMPLETE first token. A suffix such as 'LOH Region' is allowed.
  token <- strsplit(rest,'[[:space:]]+',perl=TRUE)[[1]][1]
  number_pattern <- '^[+-]?(?:[0-9]+(?:[.,][0-9]+)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?$'
  if (!grepl(number_pattern,token,perl=TRUE)) fail('Malformed CNV Confidence ',sQuote(token),': ',context)
  value <- as.numeric(chartr(',','.',token))
  if (!is.finite(value) || value<0) fail('Invalid CNV Confidence: ',context)
  list(value=value,token=token,status='parsed')
}

main <- function(args=commandArgs(trailingOnly=TRUE)) {
  config <- read_args(args)
  if (!requireNamespace('xml2',quietly=TRUE)) fail('Required R package: xml2.')
  if (!file.exists(config$xml)) fail('XML not found: ',config$xml)
  doc <- xml2::read_xml(config$xml)
  nodes <- xml2::xml_find_all(doc,"//*[local-name()='bookmark']")
  if (!length(nodes)) fail('No bookmark elements in XML: ',config$xml)
  field <- function(name) vapply(seq_along(nodes),function(i) {
    found <- xml2::xml_find_all(nodes[[i]],paste0("./*[local-name()='",name,"']"))
    if (length(found)!=1L) fail('Expected exactly one ',name,' in bookmark ',i)
    trimws(xml2::xml_text(found))
  },character(1))
  # A absent comment is treated as missing confidence under the explicit policy.
  comments <- vapply(seq_along(nodes),function(i) {
    found <- xml2::xml_find_all(nodes[[i]],"./*[local-name()='comment']")
    if (length(found)>1L) fail('Multiple comments in bookmark ',i)
    if (!length(found)) '' else trimws(xml2::xml_text(found))
  },character(1))
  ids <- field('sample_id'); chrom <- field('chr_num')
  if (any(!nzchar(ids))) fail('Empty sample_id in XML.')
  number <- function(raw, name) {
    valid <- grepl('^[+-]?[0-9]+([.,][0-9]+)?$',raw)
    if (any(!valid)) fail('Invalid ',name,' in bookmark(s): ',paste(which(!valid),collapse=', '))
    x <- as.numeric(chartr(',','.',raw))
    if (any(!is.finite(x))) fail('Non-finite ',name)
    x
  }
  start <- number(field('base_start_pos'),'base_start_pos')
  end <- number(field('base_end_pos'),'base_end_pos')
  value <- number(field('value'),'value')
  invalid <- start<1 | end<start | start!=floor(start) | end!=floor(end)
  if (any(invalid)) fail('Invalid genomic coordinates in bookmark(s): ',paste(which(invalid),collapse=', '))
  if (any(!value %in% 0:4)) fail('Unexpected cnvPartition value (expected 0..4).')
  if (any(!chrom %in% c(as.character(1:22),'X','Y','XY','M','MT')))
    fail('Unrecognized chromosome label; not automatically reinterpreted: ',paste(unique(chrom[!chrom %in% c(as.character(1:22),'X','Y','XY','M','MT')]),collapse=', '))
  parsed <- lapply(seq_along(nodes),function(i) parse_confidence(comments[i],
    paste0('bookmark ',i,' / ',ids[i],' / chr',chrom[i],':',start[i],'-',end[i]),config$missing))
  confidence <- vapply(parsed,function(x) x$value,numeric(1))
  extracted <- data.frame(sample_id=ids,chr_num=chrom,base_start_pos=start,
    base_end_pos=end,value=value,confidence=confidence,stringsAsFactors=FALSE)

  # Assou et al. 2020, Table 1, page 4; GRCh37A stated in Methods, page 7.
  # All 20 intervals checked against the publication and supplied script003.
  hotspot_locations <- data.frame(
    hot_chr=c("20","12","17","X","1","5","18","17","7","9",
              "11","13","16","1","8","14","6","15","3","22"),
    hot_start=c(29848383,11937418,51281495,1,172900000,
                104500000,56200000,7211004,132600000,68700000,
                2800000,87700000,1,16200000,93300000,
                19152018,130300000,67200000,60799,24268025),
    hot_end=c(30754613,25556120,52293893,60600000,185800000,
              117404202,61600000,8044174,133785759,114900000,
              10700000,101700000,90354753,17074942,127230818,
              107349540,139000000,67300000,26400000,49570503))
  size <- end-start
  # Overlap includes containment and touching interval boundaries.
  in_hotspot <- vapply(seq_along(nodes),function(i) any(
    hotspot_locations$hot_chr==chrom[i] &
    start[i]<=hotspot_locations$hot_end & end[i]>=hotspot_locations$hot_start),logical(1))
  # Preserve the supplied value=2 policy: standard size threshold even in hotspots.
  hotspot_threshold_applied <- in_hotspot & value!=2
  size_limit <- ifelse(hotspot_threshold_applied,config$hotspot,config$standard)
  pass_confidence <- !is.na(confidence) & confidence>config$confidence
  # No sex-based decision before sample metadata are available.
  # Script003 applies the supplied boss policy to the COMPLETE extracted table:
  # male X/Y value=1 excluded; male X value=2 excluded; other Y calls evaluated.
  keep <- pass_confidence & size>size_limit
  # Preserve original filtered schema/order; no other rows are removed on extraction.
  filtered <- extracted[keep,,drop=FALSE]
  filtered$CNV_size <- size[keep]
  filtered$hotspot_region <- in_hotspot[keep]
  audit <- data.frame(bookmark_index=seq_along(nodes),sample_id=ids,chr_num=chrom,
    base_start_pos=start,base_end_pos=end,comment_raw=comments,
    confidence_token=vapply(parsed,function(x)x$token,character(1)),confidence=confidence,
    confidence_status=vapply(parsed,function(x)x$status,character(1)),
    passes_confidence=pass_confidence,
    CNV_size=size,hotspot_region=in_hotspot,hotspot_threshold_applied=hotspot_threshold_applied,required_size=size_limit,
    passes_size=size>size_limit,keep_step001=keep,stringsAsFactors=FALSE)
  dir.create(config$out,recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(config$out)) fail('Cannot create output directory: ',config$out)
  prefix <- tools::file_path_sans_ext(basename(config$xml))
  save <- function(x,suffix) write.csv2(x,file.path(config$out,paste0(prefix,suffix)),row.names=FALSE,quote=TRUE,na='')
  save(extracted,'_extracted_data.csv'); save(filtered,'_filtered_data.csv'); save(audit,'_extraction_audit.csv')
  writeLines(c(paste('Dataset:',prefix),paste('Bookmarks:',nrow(extracted)),
    paste('XML sample IDs:',length(unique(ids))),paste('Missing confidence:',sum(is.na(confidence))),
    paste('Step001 retained:',nrow(filtered)),paste('Confidence >',config$confidence),
    paste('Size >',config$hotspot,'inside GRCh37 hotspots for value !=2; >',config$standard,'otherwise.'),
    'CNV_size=end-start; full inclusive interval overlap; 20 source-verified GRCh37 regions.',
    'Step001 contains confidence/size candidates only; no autonomous X/Y exclusion.',
    'Use script003 final CSV: sex-aware rules are applied there from script002 metadata.',
    'This filtered CSV is intermediate; final GRCh37/sex-aware policy belongs to script003.'),
    file.path(config$out,paste0(prefix,'_extraction_QC.txt')))
  if (anyNA(confidence)) warning(sum(is.na(confidence)),' missing confidence values retained in extraction and explicitly excluded from step001.',call.=FALSE)
  cat('Extracted:',nrow(extracted),'| XML IDs:',length(unique(ids)),'| Step001 retained:',nrow(filtered),'\n')
  invisible(extracted)
}
if (sys.nframe()==0L) main()
