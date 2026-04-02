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

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
}

library(DBI)
library(RSQLite)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "provirusDonor"
db_path <- if (length(args) >= 2) args[[2]] else "hiv_predictions.db"

if (!file.exists(input_path)) {
  stop(sprintf("Input path not found: %s", input_path))
}

normalize_sequence_line <- function(line) {
  cleaned <- toupper(trimws(line))
  cleaned <- gsub("\\s+", "", cleaned, perl = TRUE)
  cleaned <- gsub("\\*", "", cleaned)
  cleaned <- gsub("[^A-Z]", "", cleaned)
  cleaned
}

extract_donor_id <- function(lines, file_path) {
  if (length(lines) > 0) {
    first <- trimws(lines[[1]])
    m <- regexec("^Donor\\s*([A-Z0-9]+)$", first, ignore.case = TRUE, perl = TRUE)
    captured <- regmatches(first, m)[[1]]
    if (length(captured) == 2) {
      return(toupper(captured[[2]]))
    }
  }

  file_name <- basename(file_path)
  m_file <- regexec("^Donor([A-Z0-9]+)\\.fasta$", file_name, ignore.case = TRUE, perl = TRUE)
  cap_file <- regmatches(file_name, m_file)[[1]]
  if (length(cap_file) == 2) {
    return(toupper(cap_file[[2]]))
  }

  stop(sprintf("Could not parse donor id from file: %s", file_path))
}

normalize_header <- function(header_line) {
  cleaned <- trimws(header_line)
  cleaned <- sub("^>", "", cleaned)
  cleaned <- sub("^&gt;", "", cleaned, ignore.case = TRUE)
  tolower(trimws(cleaned))
}

parse_provirus_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) < 2) {
    stop(sprintf("Input file is too short: %s", path))
  }

  donor <- extract_donor_id(lines, path)
  sections <- list(gag = character(0), pol = character(0))
  current_section <- NULL

  for (line in lines[-1]) {
    trimmed <- trimws(line)
    if (!nzchar(trimmed)) {
      next
    }

    if (startsWith(trimmed, ">") || startsWith(tolower(trimmed), "&gt;")) {
      header <- normalize_header(trimmed)
      if (grepl("gag", header, fixed = TRUE)) {
        current_section <- "gag"
      } else if (grepl("pol", header, fixed = TRUE)) {
        current_section <- "pol"
      } else {
        current_section <- NULL
      }
      next
    }

    if (is.null(current_section)) {
      next
    }

    seq_part <- normalize_sequence_line(trimmed)
    if (nzchar(seq_part)) {
      sections[[current_section]] <- c(sections[[current_section]], seq_part)
    }
  }

  gag <- paste0(sections$gag, collapse = "")
  pol <- paste0(sections$pol, collapse = "")

  if (!nzchar(gag) || !nzchar(pol)) {
    stop(sprintf("Missing gag or pol section in file: %s", path))
  }

  list(donor = donor, gag = gag, pol = pol)
}

resolve_input_files <- function(path_value) {
  if (dir.exists(path_value)) {
    files <- list.files(path_value, pattern = "\\.fasta$", full.names = TRUE)
    if (length(files) == 0) {
      stop(sprintf("No .fasta files found in directory: %s", path_value))
    }
    return(sort(files))
  }

  if (file.exists(path_value)) {
    return(normalizePath(path_value, winslash = "/", mustWork = FALSE))
  }

  stop(sprintf("Input path does not exist: %s", path_value))
}

input_files <- resolve_input_files(input_path)
records <- lapply(input_files, parse_provirus_file)

build_pipeline_json <- function(record) {
  jsonlite::toJSON(
    list(
      id = record$donor,
      type = "result",
      data = list(
        results = list(
          list(
            type = "input_sequence_table",
            table_columns = list(
              list(
                name = "sequence_number",
                display_name = "seq #",
                type = "int",
                source = "core",
                hidden = FALSE
              ),
              list(
                name = "sequence_name",
                display_name = "sequence name",
                type = "text",
                source = "core",
                hidden = FALSE
              ),
              list(
                name = "sequence",
                display_name = "sequence",
                type = "text",
                source = "core",
                hidden = FALSE
              )
            ),
            table_data = list(
              list(1, "Gag", record$gag),
              list(2, "Pol", record$pol)
            )
          )
        ),
        errors = list(),
        warnings = list()
      ),
      status = "done"
    ),
    auto_unbox = TRUE,
    null = "null"
  )
}

con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

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

upsert_sql <- "
  INSERT INTO provirus_donor_sequence_regions (donor, gag, pol, pipeline_result_json, source_file, updated_at)
  VALUES (?, ?, ?, ?, ?, datetime('now'))
  ON CONFLICT(donor) DO UPDATE SET
    gag = excluded.gag,
    pol = excluded.pol,
    pipeline_result_json = excluded.pipeline_result_json,
    source_file = excluded.source_file,
    updated_at = excluded.updated_at
"

for (i in seq_along(records)) {
  rec <- records[[i]]
  file_path <- normalizePath(input_files[[i]], winslash = "/", mustWork = FALSE)
  pipeline_result_json <- build_pipeline_json(rec)
  dbExecute(con, upsert_sql, params = list(rec$donor, rec$gag, rec$pol, pipeline_result_json, file_path))
}

cat(sprintf(
  "Upserted %d provirus donor sequence record(s) into %s\n",
  length(records),
  normalizePath(db_path, winslash = "/", mustWork = FALSE)
))

for (rec in records) {
  cat(sprintf("  Donor %s: gag=%d aa, pol=%d aa\n", rec$donor, nchar(rec$gag), nchar(rec$pol)))
}
