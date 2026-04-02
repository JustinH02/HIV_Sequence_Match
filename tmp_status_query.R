suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(jsonlite))
con <- dbConnect(SQLite(), 'hiv_predictions.db')
on.exit(dbDisconnect(con), add = TRUE)
q <- dbGetQuery(con, "SELECT result_id, alleles_json FROM hlp_donor_pipeline_results WHERE final_status='timeout'")
if (nrow(q) == 0) {
  cat('TIMEOUT_ROWS=0\n')
  quit(status=0)
}
parse_alleles <- function(x) {
  a <- tryCatch(fromJSON(x), error=function(e) character())
  if (is.null(a)) a <- character()
  as.character(a)
}
all_alleles <- unlist(lapply(q$alleles_json, parse_alleles), use.names = FALSE)
all_alleles <- unique(all_alleles)
non_drb <- all_alleles[!grepl('^HLA-DRB', all_alleles, ignore.case = TRUE)]
rows_with_non_drb <- 0
for (i in seq_len(nrow(q))) {
  a <- parse_alleles(q$alleles_json[[i]])
  if (any(!grepl('^HLA-DRB', a, ignore.case = TRUE))) rows_with_non_drb <- rows_with_non_drb + 1
}
cat(sprintf('TIMEOUT_ROWS=%d\n', nrow(q)))
cat(sprintf('ROWS_WITH_NON_DRB=%d\n', rows_with_non_drb))
cat(sprintf('UNIQUE_TIMEOUT_ALLELES=%d\n', length(all_alleles)))
cat(sprintf('UNIQUE_NON_DRB_ALLELES=%d\n', length(non_drb)))
if (length(non_drb) > 0) {
  cat('NON_DRB_LIST=')
  cat(paste(non_drb, collapse=';'))
  cat('\n')
}
