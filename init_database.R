user_lib_env <- Sys.getenv("R_LIBS_USER")
win_user_lib_root <- file.path(path.expand("~"), "R", "win-library")
version_dirs <- if (dir.exists(win_user_lib_root)) {
  list.dirs(win_user_lib_root, recursive = FALSE, full.names = TRUE)
} else {
  character(0)
}

candidate_libs <- unique(c(user_lib_env, version_dirs))
candidate_libs <- candidate_libs[nzchar(candidate_libs) & dir.exists(candidate_libs)]

if (length(candidate_libs) > 0) {
  .libPaths(c(candidate_libs, .libPaths()))
}

if (!requireNamespace("DBI", quietly = TRUE)) {
  stop("Package 'DBI' is required. Install with: install.packages('DBI')")
}

if (!requireNamespace("RSQLite", quietly = TRUE)) {
  stop("Package 'RSQLite' is required. Install with: install.packages('RSQLite')")
}

library(DBI)
library(RSQLite)

args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[1] else "hiv_predictions.db"

con <- dbConnect(SQLite(), dbname = db_path)

add_column_if_missing <- function(connection, table_name, column_name, column_type) {
  cols <- dbGetQuery(connection, sprintf("PRAGMA table_info(%s)", table_name))
  if (!(column_name %in% cols$name)) {
    dbExecute(
      connection,
      sprintf("ALTER TABLE %s ADD COLUMN %s %s", table_name, column_name, column_type)
    )
  }
}

# Table 1: Donor (String), Alleles (String)
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS donor (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    donor TEXT NOT NULL,
    alleles TEXT NOT NULL
  )
")

# Table 2: HLP(String), Alleles(String), Peptides(String), start(int), end(int), protein region(String)
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS hlp (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    hlp TEXT NOT NULL,
    alleles TEXT NOT NULL,
    peptides TEXT NOT NULL,
    \"start\" INTEGER NOT NULL,
    \"end\" INTEGER NOT NULL,
    protein_region TEXT
  )
")

# Table 3: Provirus Donor(String), Alleles(String), Peptides(String), start(int), end(int), protein region(String)
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS provirus_donor (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    provirus_donor TEXT NOT NULL,
    alleles TEXT NOT NULL,
    peptides TEXT NOT NULL,
    \"start\" INTEGER NOT NULL,
    \"end\" INTEGER NOT NULL,
    protein_region TEXT
  )
")

# Table 4: HLP primary-keyed sequence regions (gag/pol/env)
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS hlp_sequence_regions (
    hlp TEXT PRIMARY KEY,
    gag TEXT NOT NULL,
    pol TEXT NOT NULL,
    env TEXT NOT NULL,
    pipeline_result_json TEXT,
    source_file TEXT,
    updated_at TEXT NOT NULL
  )
")

# Table 5: Provirus donor primary-keyed sequence regions (gag/pol)
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS provirus_donor_sequence_regions (
    donor TEXT PRIMARY KEY,
    gag TEXT NOT NULL,
    pol TEXT NOT NULL,
    pipeline_result_json TEXT,
    source_file TEXT,
    updated_at TEXT NOT NULL
  )
")

# Table 6: Batch API result storage for HLP x donor x region runs.
dbExecute(con, "
  CREATE TABLE IF NOT EXISTS hlp_donor_pipeline_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    hlp TEXT NOT NULL,
    donor TEXT NOT NULL,
    protein_region TEXT NOT NULL,
    alleles_json TEXT NOT NULL,
    request_payload_json TEXT NOT NULL,
    submission_http_status INTEGER,
    submission_response_json TEXT,
    result_id TEXT,
    final_status TEXT,
    result_json TEXT,
    error_message TEXT,
    started_at TEXT NOT NULL,
    completed_at TEXT
  )
")

# Backward-compatible migration for existing databases.
add_column_if_missing(con, "hlp", "protein_region", "TEXT")
add_column_if_missing(con, "provirus_donor", "protein_region", "TEXT")
add_column_if_missing(con, "hlp_sequence_regions", "pipeline_result_json", "TEXT")
add_column_if_missing(con, "provirus_donor_sequence_regions", "pipeline_result_json", "TEXT")

# Optional indexes for faster filtering.
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_donor_donor_alleles ON donor(donor, alleles)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_hlp_alleles ON hlp(hlp, alleles)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_provirus_donor_alleles ON provirus_donor(provirus_donor, alleles)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_sequence_regions_updated_at ON hlp_sequence_regions(updated_at)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_provirus_donor_sequence_regions_updated_at ON provirus_donor_sequence_regions(updated_at)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_donor_pipeline_lookup ON hlp_donor_pipeline_results(hlp, donor, protein_region)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_donor_pipeline_result_id ON hlp_donor_pipeline_results(result_id)")

dbDisconnect(con)

cat(sprintf("Database initialized at: %s\n", normalizePath(db_path, winslash = "/", mustWork = FALSE)))
