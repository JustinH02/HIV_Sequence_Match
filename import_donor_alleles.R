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
input_path <- if (length(args) >= 1) args[[1]] else "Donor.txt"
db_path <- if (length(args) >= 2) args[[2]] else "hiv_predictions.db"

if (!file.exists(input_path)) {
  stop(sprintf("Input donor file not found: %s", input_path))
}

format_allele_digits <- function(digits) {
  n <- nchar(digits)
  if (n < 4) {
    return(digits)
  }

  if (n == 4) {
    return(sprintf("%s:%s", substr(digits, 1, 2), substr(digits, 3, 4)))
  }

  sprintf("%s:%s", substr(digits, 1, n - 2), substr(digits, n - 1, n))
}

normalize_single_chain <- function(token) {
  parts <- strsplit(token, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(NA_character_)
  }

  locus <- parts[[1]]
  digits <- gsub("[^0-9]", "", parts[[2]])
  if (!nzchar(locus) || !nzchar(digits)) {
    return(NA_character_)
  }

  sprintf("HLA-%s*%s", locus, format_allele_digits(digits))
}

normalize_haplotype <- function(token) {
  m <- regexec("^HLA-([A-Z]+[0-9])([0-9]{4,5})-([A-Z]+[0-9])([0-9]{4,5})$", token, perl = TRUE)
  captured <- regmatches(token, m)[[1]]

  if (length(captured) != 5) {
    return(NA_character_)
  }

  locus_a <- captured[[2]]
  digits_a <- captured[[3]]
  locus_b <- captured[[4]]
  digits_b <- captured[[5]]

  sprintf(
    "HLA-%s*%s/%s*%s",
    locus_a,
    format_allele_digits(digits_a),
    locus_b,
    format_allele_digits(digits_b)
  )
}

normalize_allele <- function(token) {
  cleaned <- toupper(trimws(token))
  cleaned <- gsub("&QUOT;", "", cleaned, fixed = TRUE)
  cleaned <- gsub('"', "", cleaned, fixed = TRUE)
  cleaned <- gsub("[[:punct:]]+$", "", cleaned, perl = TRUE)
  cleaned <- gsub("^[[:punct:]]+", "", cleaned, perl = TRUE)
  cleaned <- gsub("_O", "_0", cleaned, fixed = TRUE)

  if (!nzchar(cleaned)) {
    return(NA_character_)
  }

  if (grepl("^HLA-[A-Z0-9]+\\*[0-9:]+(/[A-Z0-9]+\\*[0-9:]+)?$", cleaned)) {
    return(cleaned)
  }

  if (grepl("^HLA-[A-Z0-9]+[0-9]{4,5}-[A-Z0-9]+[0-9]{4,5}$", cleaned)) {
    return(normalize_haplotype(cleaned))
  }

  if (grepl("^[A-Z0-9]+_[0-9O]{4,5}$", cleaned)) {
    return(normalize_single_chain(cleaned))
  }

  NA_character_
}

extract_candidate_tokens <- function(line) {
  if (!nzchar(trimws(line))) {
    return(character(0))
  }

  line <- gsub("&quot;", '"', line, fixed = TRUE)
  normalized_line <- gsub(",", " ", line, fixed = TRUE)
  normalized_line <- gsub("\t", " ", normalized_line, fixed = TRUE)

  pieces <- unlist(strsplit(normalized_line, "\\s+", perl = TRUE), use.names = FALSE)
  pieces[nzchar(pieces)]
}

parse_donor_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  donor_map <- list()
  current_donor <- NULL

  for (line in lines) {
    cleaned_line <- trimws(line)

    # Accept both "Donor 110 -" and "0036 -" styles.
    header_match <- regexec("^(?:DONOR\\s+)?([0-9]+)\\s*-\\s*(.*)$", toupper(cleaned_line), perl = TRUE)
    captured <- regmatches(toupper(cleaned_line), header_match)[[1]]

    if (length(captured) > 0) {
      donor_id <- captured[[2]]
      trailing <- trimws(sub("^(?:Donor\\s+)?[0-9]+\\s*-\\s*", "", cleaned_line, perl = TRUE))

      current_donor <- donor_id
      if (is.null(donor_map[[current_donor]])) {
        donor_map[[current_donor]] <- character(0)
      }

      if (nzchar(trailing)) {
        for (token in extract_candidate_tokens(trailing)) {
          normalized <- normalize_allele(token)
          if (!is.na(normalized)) {
            donor_map[[current_donor]] <- c(donor_map[[current_donor]], normalized)
          }
        }
      }

      next
    }

    if (is.null(current_donor)) {
      next
    }

    tokens <- extract_candidate_tokens(cleaned_line)
    if (length(tokens) == 0) {
      next
    }

    for (token in tokens) {
      normalized <- normalize_allele(token)
      if (!is.na(normalized)) {
        donor_map[[current_donor]] <- c(donor_map[[current_donor]], normalized)
      }
    }
  }

  donors <- sort(names(donor_map))

  rows <- lapply(donors, function(donor_id) {
    alleles <- sort(unique(donor_map[[donor_id]]))
    data.frame(
      donor = donor_id,
      allele_count = length(alleles),
      alleles_json = toJSON(alleles, auto_unbox = FALSE),
      stringsAsFactors = FALSE
    )
  })

  if (length(rows) == 0) {
    return(data.frame(donor = character(0), allele_count = integer(0), alleles_json = character(0), stringsAsFactors = FALSE))
  }

  do.call(rbind, rows)
}

records <- parse_donor_file(input_path)

if (nrow(records) == 0) {
  stop("No donor records were parsed from the input file.")
}

build_donor_allele_rows <- function(parsed_records) {
  rows <- list()

  for (i in seq_len(nrow(parsed_records))) {
    donor_id <- parsed_records$donor[[i]]
    alleles <- tryCatch(fromJSON(parsed_records$alleles_json[[i]]), error = function(e) character(0))
    alleles <- sort(unique(trimws(as.character(alleles))))
    alleles <- alleles[nzchar(alleles)]

    if (length(alleles) == 0) {
      next
    }

    rows[[length(rows) + 1]] <- data.frame(
      donor = rep(donor_id, length(alleles)),
      allele = alleles,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame(donor = character(0), allele = character(0), stringsAsFactors = FALSE))
  }

  do.call(rbind, rows)
}

donor_allele_rows <- build_donor_allele_rows(records)

con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

donor_cols <- dbGetQuery(con, "PRAGMA table_info(donor)")
if (nrow(donor_cols) > 0 && !("allele" %in% donor_cols$name)) {
  backup_table_name <- sprintf("donor_schema_backup_%s", format(Sys.time(), "%Y%m%d%H%M%S"))
  cat(sprintf("Schema migration: replacing donor table. Old table backed up as %s\n", backup_table_name))
  dbExecute(con, sprintf("ALTER TABLE donor RENAME TO %s", backup_table_name))
}

dbExecute(con, "
  CREATE TABLE IF NOT EXISTS donor (
    donor TEXT NOT NULL,
    allele TEXT NOT NULL,
    source_file TEXT,
    updated_at TEXT NOT NULL,
    PRIMARY KEY (donor, allele)
  )
")

dbExecute(con, "
  CREATE TABLE IF NOT EXISTS donor_allele_list (
    donor TEXT PRIMARY KEY,
    alleles_json TEXT NOT NULL,
    allele_count INTEGER NOT NULL,
    source_file TEXT,
    updated_at TEXT NOT NULL
  )
")

dbExecute(con, "
  CREATE TABLE IF NOT EXISTS distinct_alleles (
    allele TEXT PRIMARY KEY,
    updated_at TEXT NOT NULL
  )
")

# Keep distinct_alleles synchronized with donor on all row changes.
dbExecute(con, "DROP TRIGGER IF EXISTS trg_distinct_alleles_after_insert")
dbExecute(con, "DROP TRIGGER IF EXISTS trg_distinct_alleles_after_update")
dbExecute(con, "DROP TRIGGER IF EXISTS trg_distinct_alleles_after_delete")

dbExecute(con, "
  CREATE TRIGGER trg_distinct_alleles_after_insert
  AFTER INSERT ON donor
  BEGIN
    DELETE FROM distinct_alleles;
    INSERT OR REPLACE INTO distinct_alleles (allele, updated_at)
    SELECT DISTINCT d.allele, datetime('now')
    FROM donor d
    WHERE d.allele IS NOT NULL AND trim(d.allele) <> '';
  END
")

dbExecute(con, "
  CREATE TRIGGER trg_distinct_alleles_after_update
  AFTER UPDATE ON donor
  BEGIN
    DELETE FROM distinct_alleles;
    INSERT OR REPLACE INTO distinct_alleles (allele, updated_at)
    SELECT DISTINCT d.allele, datetime('now')
    FROM donor d
    WHERE d.allele IS NOT NULL AND trim(d.allele) <> '';
  END
")

dbExecute(con, "
  CREATE TRIGGER trg_distinct_alleles_after_delete
  AFTER DELETE ON donor
  BEGIN
    DELETE FROM distinct_alleles;
    INSERT OR REPLACE INTO distinct_alleles (allele, updated_at)
    SELECT DISTINCT d.allele, datetime('now')
    FROM donor d
    WHERE d.allele IS NOT NULL AND trim(d.allele) <> '';
  END
")

dbExecute(con, "DELETE FROM donor")

if (nrow(donor_allele_rows) > 0) {
  donor_insert_sql <- "
    INSERT INTO donor (donor, allele, source_file, updated_at)
    VALUES (?, ?, ?, datetime('now'))
  "

  for (i in seq_len(nrow(donor_allele_rows))) {
    dbExecute(
      con,
      donor_insert_sql,
      params = list(
        donor_allele_rows$donor[[i]],
        donor_allele_rows$allele[[i]],
        normalizePath(input_path, winslash = "/", mustWork = FALSE)
      )
    )
  }
}

upsert_sql <- "
  INSERT INTO donor_allele_list (donor, alleles_json, allele_count, source_file, updated_at)
  VALUES (?, ?, ?, ?, datetime('now'))
  ON CONFLICT(donor) DO UPDATE SET
    alleles_json = excluded.alleles_json,
    allele_count = excluded.allele_count,
    source_file = excluded.source_file,
    updated_at = excluded.updated_at
"

for (i in seq_len(nrow(records))) {
  dbExecute(
    con,
    upsert_sql,
    params = list(
      records$donor[[i]],
      records$alleles_json[[i]],
      records$allele_count[[i]],
      normalizePath(input_path, winslash = "/", mustWork = FALSE)
    )
  )
}

# Initial (or full) refresh so the table is correct after each import run.
dbExecute(con, "DELETE FROM distinct_alleles")
dbExecute(con, "
  INSERT OR REPLACE INTO distinct_alleles (allele, updated_at)
  SELECT DISTINCT allele, datetime('now')
  FROM donor
  WHERE allele IS NOT NULL AND trim(allele) <> ''
")

cat(sprintf("Imported %d donors into donor_allele_list in %s\n", nrow(records), normalizePath(db_path, winslash = "/", mustWork = FALSE)))
cat(sprintf("Imported %d donor-allele rows into donor\n", nrow(donor_allele_rows)))
cat(sprintf("Synchronized %d distinct alleles into distinct_alleles\n", dbGetQuery(con, "SELECT COUNT(*) AS n FROM distinct_alleles")$n[[1]]))
