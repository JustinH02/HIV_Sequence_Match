library(DBI)
library(RSQLite)
library(jsonlite)

con <- dbConnect(SQLite(), "hiv_predictions.db")

rows <- dbGetQuery(con, "SELECT hlp, donor, protein_region, error_message, alleles_json, started_at FROM hlp_donor_pipeline_results WHERE final_status = 'error' ORDER BY started_at")

cat("All error rows:\n")
for (i in seq_len(nrow(rows))) {
  cat("---\n")
  cat("HLP:", rows$hlp[i], " donor:", rows$donor[i], " region:", rows$protein_region[i], " started:", rows$started_at[i], "\n")
  cat("error_message:", rows$error_message[i], "\n")
  alleles <- tryCatch(fromJSON(rows$alleles_json[i]), error = function(e) character(0))
  cat("alleles:", paste(alleles, collapse = ", "), "\n")
}

dbDisconnect(con)
