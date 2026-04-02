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
input_path <- if (length(args) >= 1) args[[1]] else "HLP4713.fasta"
db_path <- if (length(args) >= 2) args[[2]] else "hiv_predictions.db"

if (!file.exists(input_path)) {
  stop(sprintf("Input HLP file not found: %s", input_path))
}

extract_hlp_id <- function(first_line) {
  cleaned <- trimws(first_line)
  m <- regexec("^HLP\\s*([A-Z0-9]+)$", toupper(cleaned), perl = TRUE)
  captured <- regmatches(toupper(cleaned), m)[[1]]

  if (length(captured) == 2) {
    return(captured[[2]])
  }

  m_alt <- regexec("^HLP[_\\- ]?([A-Z0-9]+)$", toupper(cleaned), perl = TRUE)
  captured_alt <- regmatches(toupper(cleaned), m_alt)[[1]]
  if (length(captured_alt) == 2) {
    return(captured_alt[[2]])
  }

  stop("Could not parse HLP id from first line. Expected format like: HLP 4713 or HLP UH7")
}

normalize_header <- function(header_line) {
  cleaned <- trimws(header_line)
  cleaned <- sub("^>", "", cleaned)
  cleaned <- sub("^&gt;", "", cleaned, ignore.case = TRUE)
  tolower(trimws(cleaned))
}

normalize_sequence_line <- function(line) {
  cleaned <- toupper(trimws(line))
  cleaned <- gsub("\\s+", "", cleaned, perl = TRUE)
  # Remove terminal stop markers and any non-amino-acid placeholders.
  cleaned <- gsub("\\*", "", cleaned)
  cleaned <- gsub("[^A-Z]", "", cleaned)
  cleaned
}

parse_hlp_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) < 2) {
    stop("Input file is too short to contain HLP id and protein sections.")
  }

  hlp_id <- extract_hlp_id(lines[[1]])
  sections <- list(gag = character(0), pol = character(0), env = character(0))
  current_section <- NULL

  for (line in lines[-1]) {
    trimmed <- trimws(line)
    if (!nzchar(trimmed)) {
      next
    }

    if (startsWith(trimmed, ">") || startsWith(tolower(trimmed), "&gt;")) {
      hdr <- normalize_header(trimmed)
      if (hdr %in% c("gag", "pol", "env")) {
        current_section <- hdr
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
  env <- paste0(sections$env, collapse = "")

  if (!nzchar(gag) || !nzchar(pol) || !nzchar(env)) {
    stop("Missing one or more required sections (gag, pol, env) in input file.")
  }

  list(hlp = hlp_id, gag = gag, pol = pol, env = env)
}

record <- parse_hlp_file(input_path)

build_pipeline_json <- function(record) {
  jsonlite::toJSON(
    list(
      id = record$hlp,
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
              list(2, "Pol", record$pol),
              list(3, "Env", record$env)
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

pipeline_result_json <- build_pipeline_json(record)

con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

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

dbExecute(
  con,
  "
    INSERT INTO hlp_sequence_regions (hlp, gag, pol, env, pipeline_result_json, source_file, updated_at)
    VALUES (?, ?, ?, ?, ?, ?, datetime('now'))
    ON CONFLICT(hlp) DO UPDATE SET
      gag = excluded.gag,
      pol = excluded.pol,
      env = excluded.env,
      pipeline_result_json = excluded.pipeline_result_json,
      source_file = excluded.source_file,
      updated_at = excluded.updated_at
  ",
  params = list(
    record$hlp,
    record$gag,
    record$pol,
    record$env,
    pipeline_result_json,
    normalizePath(input_path, winslash = "/", mustWork = FALSE)
  )
)

cat(sprintf(
  "Upserted HLP %s into hlp_sequence_regions in %s (gag=%d aa, pol=%d aa, env=%d aa)\n",
  record$hlp,
  normalizePath(db_path, winslash = "/", mustWork = FALSE),
  nchar(record$gag),
  nchar(record$pol),
  nchar(record$env)
))