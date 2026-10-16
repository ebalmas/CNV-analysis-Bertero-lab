#!/usr/bin/env Rscript
# Usage (original interface preserved):
# Rscript 002_build_CNV_sample_metadata.R SHEET1.csv [SHEET2.csv ...] OUTPUT.csv
#   [--cell_line=LABEL|Sample_ID=LABEL,default=LABEL]
#   [--samples=ID1,ID2|file.txt] [--extracted=EXTRACTED.csv | --xml=BOOKMARKS.xml]
# Output: write.csv2 metadata compatible with 003 + *_sample_checks.csv inventory.
# Exact XML/extracted IDs are authoritative; NEVER infer [n] from a well.
# Without --extracted/--xml only bare-ID inventory is emitted, explicitly unverified.
# No interactive prompt: omitted cell_line is 'not recorded', never guessed.
# Semicolon/comma/tab input delimiters and quoted cells are supported.

fail <- function(...) stop(...,call.=FALSE)
clean <- function(x) { x <- as.character(x); x[is.na(x)] <- ''; trimws(x) }
strip_index <- function(x) sub('[[:space:]]+\\[[0-9]+\\][[:space:]]*$', '', clean(x))
read_args <- function(args) {
  allowed <- c('cell_line','samples','extracted','xml')
  flags <- list(); pos <- character(0)
  for (arg in args) {
    if (grepl('^--',arg)) {
      eq <- regexpr('=',arg,fixed=TRUE)[1]
      if (eq<0L) fail('Expected --flag=value: ',arg)
      key <- substring(arg,3L,eq-1L); val <- substring(arg,eq+1L)
      if (!key %in% allowed) fail('Unknown flag: --',key)
      if (!is.null(flags[[key]])) fail('Repeated flag: --',key)
      if (!nzchar(trimws(val))) fail('Empty flag: --',key)
      flags[[key]] <- val
    } else pos <- c(pos,arg)
  }
  if (length(pos)<2L) fail('Usage: Rscript 002_build_CNV_sample_metadata.R SHEET.csv [SHEET2.csv ...] OUTPUT.csv [--extracted=CALLS.csv|--xml=BOOKMARKS.xml] [--cell_line=...] [--samples=...]')
  if (!is.null(flags$xml) && !is.null(flags$extracted)) fail('Use either --xml or --extracted, not both.')
  list(sheets=head(pos,-1L),out=tail(pos,1L),flags=flags)
}

read_delimited <- function(lines, required, context) {
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) fail('Empty table: ',context)
  # Decide from parsed column names, not delimiter counts inside quoted text.
  sep_candidates <- c(';',',','\t')
  good <- vapply(sep_candidates,function(sep) {
    header <- tryCatch(read.table(text=lines[1],sep=sep,header=FALSE,quote='"',
      comment.char='',colClasses='character',check.names=FALSE),error=function(e)NULL)
    !is.null(header) && all(required %in% clean(unlist(header,use.names=FALSE)))
  },logical(1))
  if (sum(good)!=1L) fail('Cannot determine one delimiter/required header for ',context,
    '. Required: ',paste(required,collapse=', '))
  connection <- textConnection(lines)
  counts <- tryCatch(count.fields(connection,sep=sep_candidates[good],quote='"',comment.char=''),
    error=function(e) fail('Malformed quoting in ',context,': ',conditionMessage(e)),
    finally=close(connection))
  if (length(counts) && any(!is.na(counts) & counts!=counts[1]))
    fail('Unequal field counts in ',context,'; no columns are padded or truncated.')
  result <- tryCatch(read.table(text=paste(lines,collapse='\n'),sep=sep_candidates[good],
    header=TRUE,quote='"',comment.char='',colClasses='character',check.names=FALSE,
    fill=FALSE,blank.lines.skip=TRUE,na.strings='',row.names=NULL),
    error=function(e)fail('Malformed delimited table ',context,': ',conditionMessage(e)))
  names(result) <- clean(names(result))
  if (any(!nzchar(names(result))) || anyDuplicated(names(result))) fail('Empty/duplicate header names: ',context)
  result[] <- lapply(result,clean)
  if (!all(required %in% names(result))) fail('Missing required columns: ',context)
  result
}

read_samplesheet <- function(path) {
  if (!file.exists(path)) fail('Samplesheet not found: ',path)
  lines <- readLines(path,warn=FALSE,encoding='UTF-8')
  if (length(lines)) lines[1] <- sub('^\ufeff','',lines[1])
  marker <- grep('^\\[Data\\]([;,\t].*)?$',trimws(lines))
  # Also accept a plain tabular samplesheet with no Illumina preamble.
  if (length(marker)>1L) fail('Multiple [Data] sections: ',path)
  if (length(marker)) {
    if (marker>=length(lines)) fail('Empty [Data] section: ',path)
    lines <- lines[seq.int(marker+1L,length(lines))]
    next_section <- grep('^\\[',trimws(lines))
    if (length(next_section)) lines <- head(lines,next_section[1]-1L)
  }
  df <- read_delimited(lines,c('Sample_ID','Gender'),path)
  if (!nrow(df)) fail('No samples in ',path)
  if (any(!nzchar(df$Sample_ID))) fail('Empty Sample_ID in ',path)
  for (name in c('Sample_Name','SentrixBarcode_A','SentrixPosition_A'))
    if (!name %in% names(df)) df[[name]] <- ''
  # Project-dependent additional columns are irrelevant to identity matching.
  df <- df[c('Sample_ID','Sample_Name','Gender','SentrixBarcode_A','SentrixPosition_A')]
  df$source_samplesheet <- normalizePath(path,mustWork=TRUE)
  df
}

parse_cell_line <- function(arg) {
  if (is.null(arg)) return(list(mode='absent',value=NULL))
  if (!grepl('=',arg,fixed=TRUE)) return(list(mode='single',value=trimws(arg)))
  pairs <- strsplit(arg,',',fixed=TRUE)[[1]]; values <- list()
  for (pair in pairs) {
    eq <- regexpr('=',pair,fixed=TRUE)[1]
    if (eq<2L) fail('Malformed --cell_line map: ',pair)
    key <- trimws(substr(pair,1,eq-1L)); value <- trimws(substring(pair,eq+1L))
    if (!nzchar(key) || !nzchar(value) || !is.null(values[[key]])) fail('Empty/duplicate --cell_line entry: ',key)
    values[[key]] <- value
  }
  if (sum(tolower(names(values))=='default')>1L) fail('Repeated default in --cell_line map.')
  list(mode='map',value=values)
}
resolve_line <- function(id, parsed) {
  if (parsed$mode=='absent') return(list(value='not recorded',source='not supplied; no biological identity inferred'))
  if (parsed$mode=='single') return(list(value=parsed$value,source='explicit --cell_line; common label'))
  if (!is.null(parsed$value[[id]])) return(list(value=parsed$value[[id]],source=paste('explicit --cell_line entry:',id)))
  default <- names(parsed$value)[tolower(names(parsed$value))=='default']
  if (length(default)) return(list(value=parsed$value[[default]],source='explicit --cell_line default'))
  list(value='not recorded',source=paste('no --cell_line entry for',id))
}

main <- function(args=commandArgs(trailingOnly=TRUE)) {
  config <- read_args(args)
  rows <- do.call(rbind,lapply(config$sheets,read_samplesheet))
  # Collapse identical repeated records across sheets; keep combined provenance.
  identity_cols <- setdiff(names(rows),'source_samplesheet')
  unique_rows <- rows[!duplicated(rows[identity_cols]),,drop=FALSE]
  for (i in seq_len(nrow(unique_rows))) {
    equal <- rep(TRUE,nrow(rows))
    for (col in identity_cols) equal <- equal & rows[[col]]==unique_rows[[col]][i]
    unique_rows$source_samplesheet[i] <- paste(unique(rows$source_samplesheet[equal]),collapse=' | ')
  }
  rows <- unique_rows; rownames(rows) <- NULL
  rows$barcode_position <- ifelse(nzchar(rows$SentrixBarcode_A) & nzchar(rows$SentrixPosition_A),
    paste0(rows$SentrixBarcode_A,'_',rows$SentrixPosition_A),'')
  # Preserve unknown sex explicitly; do not infer sex from a sample name or CN calls.
  sex <- tolower(rows$Gender)
  rows$expected_sex <- ifelse(sex %in% c('male','m','xy'),'male',ifelse(sex %in% c('female','f','xx'),'female','unknown'))
  selected <- rep(TRUE,nrow(rows))
  if (!is.null(config$flags$samples)) {
    arg <- config$flags$samples
    wanted <- if (file.exists(arg)) clean(readLines(arg,warn=FALSE)) else clean(strsplit(arg,',',fixed=TRUE)[[1]])
    wanted <- unique(wanted[nzchar(wanted)])
    unknown <- setdiff(wanted,rows$Sample_ID)
    if (!length(wanted) || length(unknown)) fail('Invalid --samples selection; unknown IDs: ',paste(unknown,collapse=', '))
    selected <- rows$Sample_ID %in% wanted
  }
  parsed_line <- parse_cell_line(config$flags$cell_line)
  if (parsed_line$mode=='map') {
    unknown_keys <- setdiff(names(parsed_line$value)[tolower(names(parsed_line$value))!='default'],rows$Sample_ID)
    if (length(unknown_keys)) fail('--cell_line map contains unknown Sample_ID: ',paste(unknown_keys,collapse=', '))
  }
  line <- lapply(rows$Sample_ID,resolve_line,parsed=parsed_line)
  rows$cell_line <- vapply(line,function(x)x$value,character(1))
  rows$cell_line_source <- vapply(line,function(x)x$source,character(1))
  observed <- !is.null(config$flags$xml) || !is.null(config$flags$extracted)
  if (!is.null(config$flags$xml)) {
    if (!requireNamespace('xml2',quietly=TRUE)) fail('Required for --xml: xml2.')
    doc <- xml2::read_xml(config$flags$xml)
    nodes <- xml2::xml_find_all(doc,"//*[local-name()='bookmark']")
    ids <- vapply(nodes,function(node) {
      hit <- xml2::xml_find_all(node,"./*[local-name()='sample_id']")
      if (length(hit)!=1L) fail('Missing/duplicate sample_id in XML bookmark.')
      trimws(xml2::xml_text(hit))
    },character(1))
  } else if (!is.null(config$flags$extracted)) {
    path <- config$flags$extracted
    if (!file.exists(path)) fail('--extracted not found: ',path)
    df <- read_delimited(readLines(path,warn=FALSE), 'sample_id',path)
    ids <- df$sample_id
  } else {
    ids <- rows$Sample_ID[selected]
    warning('No --extracted/--xml: bare Sample_ID inventory only; bracket indices are not guessed.',call.=FALSE)
  }
  if (!length(ids) || any(!nzchar(ids))) fail('No usable sample IDs in input.')
  unique_ids <- unique(ids); bases <- strip_index(unique_ids)
  owner <- integer(length(unique_ids)); method <- character(length(unique_ids))
  for (i in seq_along(unique_ids)) {
    identity_tokens <- unique(c(unique_ids[i],bases[i]))
    matches_id <- which(rows$Sample_ID %in% identity_tokens)
    matches_barcode <- which(nzchar(rows$barcode_position) & rows$barcode_position %in% identity_tokens)
    hits <- unique(c(matches_id,matches_barcode))
    if (length(hits)!=1L) fail('Missing or ambiguous samplesheet match for ',sQuote(unique_ids[i]),
      ': ',length(hits),' physical records. Provide complete sheets with unique identities.')
    owner[i] <- hits
    method[i] <- if (length(matches_id) && length(matches_barcode)) 'Sample_ID and barcode_position' else if (length(matches_id)) 'Sample_ID' else 'barcode_position'
  }
  chosen <- selected[owner]
  if (!any(chosen)) fail('No observed sample IDs match --samples selection.')
  count <- tabulate(match(ids,unique_ids),nbins=length(unique_ids))
  matched_rows <- rows[owner[chosen],,drop=FALSE]
  metadata <- data.frame(sample_id=unique_ids[chosen],sample_name=matched_rows$Sample_Name,
    expected_sex=matched_rows$expected_sex,
    sex_source=paste0('samplesheet Gender=',matched_rows$Gender,'; ',matched_rows$source_samplesheet),
    cell_line=matched_rows$cell_line,cell_line_source=matched_rows$cell_line_source,
    sample_id_verification=if(observed) 'verified exact observed ID; no index inference' else 'UNVERIFIED bare-ID inventory; no XML supplied',
    Sample_ID=matched_rows$Sample_ID,barcode_position=matched_rows$barcode_position,
    match_method=method[chosen],stringsAsFactors=FALSE)
  if (anyDuplicated(metadata$sample_id)) fail('Duplicate output metadata identity.')
  metadata <- metadata[order(metadata$sample_id),,drop=FALSE]
  coverage <- rows[c('Sample_ID','Sample_Name','barcode_position','expected_sex','source_samplesheet')]
  coverage$selected <- selected
  coverage$observed_sample_ids <- vapply(seq_len(nrow(rows)),function(i) if(observed) paste(unique_ids[owner==i],collapse=' | ') else '',character(1))
  coverage$xml_or_extracted_calls <- vapply(seq_len(nrow(rows)),function(i) if(observed) sum(count[owner==i]) else NA_integer_,integer(1))
  coverage$status <- ifelse(!selected,'excluded_by_samples',if(!observed) 'unverified_inventory' else ifelse(coverage$xml_or_extracted_calls>0,'matched','no_calls_in_supplied_XML_or_extracted'))
  if (any(metadata$expected_sex=='unknown')) warning('Unknown expected sex in output; script003 will stop until corrected in source metadata.',call.=FALSE)
  if (any(metadata$cell_line=='not recorded')) warning('Cell line not supplied for some/all samples; marked not recorded.',call.=FALSE)
  dir.create(dirname(config$out),recursive=TRUE,showWarnings=FALSE)
  prefix <- tools::file_path_sans_ext(config$out)
  write.csv2(metadata,config$out,row.names=FALSE,quote=TRUE,na='')
  write.csv2(coverage,paste0(prefix,'_sample_checks.csv'),row.names=FALSE,quote=TRUE,na='')
  writeLines(c(paste('Samplesheets:',length(config$sheets)),paste('Unique physical records:',nrow(rows)),
    paste('Selected physical records:',sum(selected)),paste('Metadata IDs written:',nrow(metadata)),
    paste('Observed IDs excluded by --samples:',sum(!chosen)),
    if(observed) paste('Selected physical records without observed calls:',sum(selected & coverage$xml_or_extracted_calls==0)) else 'UNVERIFIED: no observed IDs supplied.',
    paste('Unknown sex in emitted metadata:',sum(metadata$expected_sex=='unknown')),
    paste('Cell line not recorded:',sum(metadata$cell_line=='not recorded')),
    'No bracket index is inferred from Sample_Well, sheet order or chip position.'),paste0(prefix,'_QC.txt'))
  cat('Physical records:',nrow(rows),'| Metadata IDs:',nrow(metadata),'| Observed verification:',observed,'\n')
  invisible(metadata)
}
if (sys.nframe()==0L) main()
