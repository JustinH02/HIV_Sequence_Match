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

required_packages <- c("DBI", "RSQLite", "jsonlite", "httr2", "future", "future.apply")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg))
  }
}

library(DBI)
library(RSQLite)
library(jsonlite)
library(httr2)
library(future)
library(future.apply)

args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[[1]] else "hiv_predictions.db"
max_hits <- if (length(args) >= 2) as.integer(args[[2]]) else 0L
if (is.na(max_hits) || max_hits < 0) {
  max_hits <- 0L
}
dry_run <- if (length(args) >= 3) tolower(args[[3]]) %in% c("1", "true", "yes", "y") else FALSE
workers <- if (length(args) >= 4) as.integer(args[[4]]) else 2L
if (is.na(workers) || workers < 1L) {
  workers <- 1L
}
resume_mode <- if (length(args) >= 5) tolower(args[[5]]) %in% c("1", "true", "yes", "y") else TRUE

api_url <- "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
results_base_url <- "https://api-nextgen-tools.iedb.org/api/v1/results"
poll_interval_seconds <- 10
max_poll_attempts <- 60
max_retries <- 3
retry_delay_seconds <- 30
job_timeout_seconds <- 360

pending_states <- c("pending", "queued", "running", "in_progress", "processing")
terminal_success_states <- c("completed", "complete", "succeeded", "success", "done")
terminal_failure_states <- c("failed", "error", "cancelled", "canceled")

now_utc <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

log_progress <- function(fmt, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  pid <- Sys.getpid()
  cat(sprintf("[%s][pid:%s] %s\n", timestamp, pid, sprintf(fmt, ...)))
  flush.console()
}

build_fasta_input <- function(region_name, sequence_text) {
  paste0(">", region_name, "\n", sequence_text)
}

find_results_id <- function(x) {
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      key_exact <- c("results_id", "result_id", "resultsId", "resultId")
      for (key in key_exact) {
        idx <- which(nms == key)
        if (length(idx) > 0) {
          val <- x[[idx[1]]]
          if (length(val) > 0 && !is.null(val)) {
            return(as.character(val[[1]]))
          }
        }
      }
    }

    for (item in x) {
      found <- find_results_id(item)
      if (!is.null(found) && nzchar(found)) {
        return(found)
      }
    }
  }

  NULL
}

find_status <- function(x) {
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      idx <- which(nms == "status")
      if (length(idx) > 0) {
        val <- x[[idx[1]]]
        if (length(val) > 0 && !is.null(val)) {
          return(as.character(val[[1]]))
        }
      }
    }

    for (item in x) {
      found <- find_status(item)
      if (!is.null(found) && nzchar(found)) {
        return(found)
      }
    }
  }

  NULL
}

get_status_response <- function(results_id, timeout_seconds = 30) {
  status_url <- sprintf("%s/%s?statusOnly=true", results_base_url, URLencode(results_id, reserved = TRUE))
  status_req <- request(status_url) |>
    req_method("GET") |>
    req_headers("Accept" = "application/json") |>
    req_timeout(timeout_seconds)

  status_resp <- req_perform(status_req)
  status_text <- resp_body_string(status_resp)
  parsed <- tryCatch(fromJSON(status_text, simplifyVector = FALSE), error = function(e) NULL)

  list(
    http_status = resp_status(status_resp),
    text = status_text,
    parsed = parsed,
    status = if (is.null(parsed)) NULL else find_status(parsed)
  )
}

get_full_results_response <- function(results_id, timeout_seconds = 30) {
  full_url <- sprintf("%s/%s", results_base_url, URLencode(results_id, reserved = TRUE))
  full_req <- request(full_url) |>
    req_method("GET") |>
    req_headers("Accept" = "application/json") |>
    req_timeout(timeout_seconds)

  full_resp <- req_perform(full_req)

  list(
    http_status = resp_status(full_resp),
    text = resp_body_string(full_resp)
  )
}

create_payload <- function(input_sequence_text, alleles_csv) {
  list(
    pipeline_id = "",
    pipeline_title = "",
    email = "",
    run_stage_range = list(1, 1),
    stages = list(
      list(
        stage_display_name = "T-Cell Prediction - Class II",
        stage_number = 1,
        stage_type = "prediction",
        tool_group = "mhcii",
        input_sequence_text = input_sequence_text,
        input_parameters = list(
          alleles = alleles_csv,
          peptide_shift = 5,
          peptide_length_range = list(15, 15),
          predictors = list(
            list(
              type = "binding",
              method = "netmhciipan_ba"
            )
          )
        ),
        table_state = list(
          list(table = "peptide_table", columns = list()),
          list(table = "consolidated_peptide_table", columns = list())
        )
      )
    )
  )
}

if (!file.exists(db_path)) {
  stop(sprintf("Database file not found: %s", db_path))
}

con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

dbExecute(con, "PRAGMA foreign_keys = ON")

dbExecute(con, "
  CREATE TABLE IF NOT EXISTS hlp_region_keys (
    hlp TEXT NOT NULL,
    protein_region TEXT NOT NULL,
    PRIMARY KEY (hlp, protein_region),
    FOREIGN KEY (hlp) REFERENCES hlp_sequence_regions(hlp)
  )
")

dbExecute(con, "
  INSERT OR IGNORE INTO hlp_region_keys (hlp, protein_region)
  SELECT hlp, 'gag' FROM hlp_sequence_regions
  UNION ALL
  SELECT hlp, 'pol' FROM hlp_sequence_regions
  UNION ALL
  SELECT hlp, 'env' FROM hlp_sequence_regions
")

dbExecute(con, "
  CREATE TABLE IF NOT EXISTS hlp_donor_pipeline_results (
    hlp TEXT NOT NULL,
    protein_region TEXT NOT NULL,
    allele TEXT NOT NULL,
    request_payload_json TEXT NOT NULL,
    submission_http_status INTEGER,
    submission_response_json TEXT,
    result_id TEXT,
    final_status TEXT,
    result_json TEXT,
    error_message TEXT,
    started_at TEXT NOT NULL,
    completed_at TEXT,
    FOREIGN KEY (hlp, protein_region) REFERENCES hlp_region_keys(hlp, protein_region),
    FOREIGN KEY (allele) REFERENCES distinct_alleles(allele)
  )
")

existing_result_cols <- dbGetQuery(con, "PRAGMA table_info(hlp_donor_pipeline_results)")
has_id_col <- "id" %in% existing_result_cols$name
has_donor_col <- "donor" %in% existing_result_cols$name
has_allele_col <- "allele" %in% existing_result_cols$name
has_alleles_json_col <- "alleles_json" %in% existing_result_cols$name

if (has_id_col || has_donor_col || !has_allele_col || has_alleles_json_col) {
  backup_table_name <- sprintf("hlp_donor_pipeline_results_schema_backup_%s", format(Sys.time(), "%Y%m%d%H%M%S"))
  cat(sprintf("Schema migration: replacing hlp_donor_pipeline_results to current schema. Old table backed up as %s\n", backup_table_name))
  dbExecute(con, sprintf("ALTER TABLE hlp_donor_pipeline_results RENAME TO %s", backup_table_name))
  dbExecute(con, "
    CREATE TABLE hlp_donor_pipeline_results (
      hlp TEXT NOT NULL,
      protein_region TEXT NOT NULL,
      allele TEXT NOT NULL,
      request_payload_json TEXT NOT NULL,
      submission_http_status INTEGER,
      submission_response_json TEXT,
      result_id TEXT,
      final_status TEXT,
      result_json TEXT,
      error_message TEXT,
      started_at TEXT NOT NULL,
      completed_at TEXT,
      FOREIGN KEY (hlp, protein_region) REFERENCES hlp_region_keys(hlp, protein_region),
      FOREIGN KEY (allele) REFERENCES distinct_alleles(allele)
    )
  ")
}

dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_donor_pipeline_lookup ON hlp_donor_pipeline_results(hlp, allele, protein_region)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hlp_donor_pipeline_result_id ON hlp_donor_pipeline_results(result_id)")
dbExecute(con, "CREATE UNIQUE INDEX IF NOT EXISTS idx_hlp_donor_pipeline_unique ON hlp_donor_pipeline_results(hlp, protein_region, allele)")

dbExecute(con, "
  DELETE FROM hlp_donor_pipeline_results
  WHERE NOT EXISTS (
    SELECT 1 FROM distinct_alleles d WHERE d.allele = hlp_donor_pipeline_results.allele
  )
")

dbExecute(con, "
  DELETE FROM hlp_donor_pipeline_results
  WHERE NOT EXISTS (
    SELECT 1
    FROM hlp_region_keys k
    WHERE k.hlp = hlp_donor_pipeline_results.hlp
      AND k.protein_region = hlp_donor_pipeline_results.protein_region
  )
")

dbExecute(con, "
  INSERT OR IGNORE INTO hlp_donor_pipeline_results (
    hlp, protein_region, allele, request_payload_json, final_status, started_at
  )
  SELECT
    k.hlp,
    k.protein_region,
    a.allele,
    '{}',
    'pending',
    strftime('%Y-%m-%dT%H:%M:%SZ', 'now')
  FROM hlp_region_keys k
  JOIN distinct_alleles a
    ON 1 = 1
")

hlp_rows <- dbGetQuery(con, "SELECT hlp, gag, pol, env FROM hlp_sequence_regions ORDER BY hlp")
distinct_allele_rows <- dbGetQuery(con, "SELECT allele FROM distinct_alleles ORDER BY allele")

if (nrow(hlp_rows) == 0) {
  stop("No rows found in hlp_sequence_regions.")
}

if (nrow(distinct_allele_rows) == 0) {
  stop("No rows found in distinct_alleles.")
}

work_items <- list()
for (i in seq_len(nrow(hlp_rows))) {
  hlp_id <- as.character(hlp_rows$hlp[[i]])
  region_map <- list(
    gag = as.character(hlp_rows$gag[[i]]),
    pol = as.character(hlp_rows$pol[[i]]),
    env = as.character(hlp_rows$env[[i]])
  )

  for (j in seq_len(nrow(distinct_allele_rows))) {
    allele_value <- trimws(as.character(distinct_allele_rows$allele[[j]]))
    if (!nzchar(allele_value)) {
      next
    }

    for (region_name in names(region_map)) {
      seq_text <- region_map[[region_name]]
      if (!nzchar(seq_text)) {
        next
      }

      work_items[[length(work_items) + 1]] <- list(
        hlp = hlp_id,
        protein_region = region_name,
        sequence_text = seq_text,
        allele = allele_value
      )
    }
  }
}

if (length(work_items) == 0) {
  stop("No valid work items were built from HLP and distinct_alleles tables.")
}

if (resume_mode) {
  success_like <- paste(sprintf("'%s'", tolower(terminal_success_states)), collapse = ",")
  completed_rows <- dbGetQuery(
    con,
    sprintf(
      "SELECT DISTINCT hlp, allele, protein_region FROM hlp_donor_pipeline_results WHERE lower(final_status) IN (%s)",
      success_like
    )
  )

  if (nrow(completed_rows) > 0) {
    done_keys <- paste(completed_rows$hlp, completed_rows$allele, completed_rows$protein_region, sep = "\r")
    work_keys <- vapply(work_items, function(x) paste(x$hlp, x$allele, x$protein_region, sep = "\r"), character(1))
    keep <- !(work_keys %in% done_keys)
    skipped <- sum(!keep)
    if (skipped > 0) {
      cat(sprintf("Resume mode: skipping %d already-successful item(s).\n", skipped))
    }
    work_items <- work_items[keep]
  }
}

if (max_hits > 0L) {
  work_items <- work_items[seq_len(min(length(work_items), max_hits))]
}

cat(sprintf("Prepared %d API hit(s) from hlp_sequence_regions x distinct_alleles x regions (gag/pol/env).\n", length(work_items)))
cat(sprintf("dry_run=%s\n", ifelse(dry_run, "TRUE", "FALSE")))
cat(sprintf("workers=%d\n", workers))

update_sql <- "
  UPDATE hlp_donor_pipeline_results
  SET
    request_payload_json = ?,
    submission_http_status = ?,
    submission_response_json = ?,
    result_id = ?,
    final_status = ?,
    result_json = ?,
    error_message = ?,
    started_at = ?,
    completed_at = ?
  WHERE hlp = ? AND protein_region = ? AND allele = ?
"

process_work_item <- function(item, idx, total_items, dry_run_mode) {
  started_at <- now_utc()
  started_at_time <- Sys.time()
  log_progress("[%d/%d] START HLP=%s allele=%s region=%s", idx, total_items, item$hlp, item$allele, item$protein_region)

  alleles_csv <- item$allele
  input_fasta <- build_fasta_input(toupper(item$protein_region), item$sequence_text)
  payload <- create_payload(input_fasta, alleles_csv)
  payload_json <- toJSON(payload, auto_unbox = TRUE, null = "null")

  if (dry_run_mode) {
    log_progress("[%d/%d] DRY_RUN prepared payload HLP=%s allele=%s region=%s", idx, total_items, item$hlp, item$allele, item$protein_region)
    return(list(
      idx = idx,
      total = total_items,
      hlp = item$hlp,
      protein_region = item$protein_region,
      allele = item$allele,
      request_payload_json = payload_json,
      submission_http_status = NA_integer_,
      submission_response_json = NA_character_,
      result_id = NA_character_,
      final_status = "dry_run",
      result_json = NA_character_,
      error_message = NA_character_,
      started_at = started_at,
      completed_at = now_utc()
    ))
  }

  submission_http_status <- NA_integer_
  submission_response_json <- NA_character_
  result_id <- NA_character_
  final_status <- NA_character_
  result_json <- NA_character_
  error_message <- NA_character_

  for (retry_attempt in seq_len(max_retries + 1)) {
    is_last_attempt <- retry_attempt == (max_retries + 1)

    elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at_time, units = "secs"))
    if (elapsed_seconds >= job_timeout_seconds) {
      final_status <- "timeout"
      error_message <- sprintf("Timed out after %d seconds waiting for pipeline response", job_timeout_seconds)
      log_progress("[%d/%d] TIMEOUT before attempt HLP=%s allele=%s region=%s", idx, total_items, item$hlp, item$allele, item$protein_region)
      break
    }

    if (retry_attempt > 1) {
      log_progress("[%d/%d] RETRY attempt=%d/%d HLP=%s allele=%s region=%s (waiting %ds)",
                   idx, total_items, retry_attempt - 1L, max_retries, item$hlp, item$allele, item$protein_region, retry_delay_seconds)
      Sys.sleep(retry_delay_seconds)
      # reset state for retry
      submission_http_status <- NA_integer_
      submission_response_json <- NA_character_
      result_id <- NA_character_
      final_status <- NA_character_
      result_json <- NA_character_
      error_message <- NA_character_
    }

    run_error <- NULL
    tryCatch({
      elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at_time, units = "secs"))
      remaining_timeout <- max(1, floor(job_timeout_seconds - elapsed_seconds))

      log_progress("[%d/%d] SUBMITTING request HLP=%s allele=%s region=%s", idx, total_items, item$hlp, item$allele, item$protein_region)
      req <- request(api_url) |>
        req_method("POST") |>
        req_headers(
          "Content-Type" = "application/json",
          "Accept" = "application/json"
        ) |>
        req_timeout(remaining_timeout) |>
        req_body_json(payload, auto_unbox = TRUE)

      resp <- req_perform(req)
      submission_http_status <- as.integer(resp_status(resp))
      submission_response_json <- resp_body_string(resp)
      log_progress("[%d/%d] SUBMITTED http_status=%s HLP=%s allele=%s region=%s", idx, total_items, submission_http_status, item$hlp, item$allele, item$protein_region)

      parsed_submit <- tryCatch(fromJSON(submission_response_json, simplifyVector = FALSE), error = function(e) NULL)
      result_id_candidate <- if (is.null(parsed_submit)) NULL else find_results_id(parsed_submit)
      if (is.null(result_id_candidate) || length(result_id_candidate) == 0 || !nzchar(as.character(result_id_candidate[[1]]))) {
        result_id <- NA_character_
        stop("Could not find results ID in submission response.")
      }
      result_id <- as.character(result_id_candidate[[1]])
      log_progress("[%d/%d] RESULT_ID=%s", idx, total_items, result_id)

      observed_status <- NULL
      for (attempt in seq_len(max_poll_attempts)) {
        elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at_time, units = "secs"))
        if (elapsed_seconds >= job_timeout_seconds) {
          observed_status <- "timeout"
          error_message <- sprintf("Timed out after %d seconds waiting for result status", job_timeout_seconds)
          log_progress("[%d/%d] POLL timeout reached result_id=%s", idx, total_items, result_id)
          break
        }

        remaining_timeout <- max(1, floor(job_timeout_seconds - elapsed_seconds))
        status_result <- get_status_response(result_id, timeout_seconds = remaining_timeout)
        status_value <- tolower(if (is.null(status_result$status)) "unknown" else status_result$status)

        if (attempt == 1L || attempt %% 6L == 0L || status_value %in% c(terminal_success_states, terminal_failure_states, "unknown")) {
          log_progress("[%d/%d] POLL attempt=%d/%d result_id=%s status=%s", idx, total_items, attempt, max_poll_attempts, result_id, status_value)
        }

        if (status_value %in% terminal_success_states) {
          observed_status <- status_value
          break
        }

        if (status_value %in% terminal_failure_states) {
          observed_status <- status_value
          break
        }

        if (!(status_value %in% pending_states) && status_value != "unknown") {
          observed_status <- status_value
          break
        }

        if (attempt < max_poll_attempts) {
          sleep_seconds <- min(poll_interval_seconds, remaining_timeout)
          if (sleep_seconds > 0) {
            Sys.sleep(sleep_seconds)
          }
        }
      }

      if (is.null(observed_status)) {
        observed_status <- "timeout"
      }

      final_status <- observed_status
      log_progress("[%d/%d] FINAL_STATUS=%s result_id=%s", idx, total_items, final_status, result_id)

      if (final_status %in% terminal_success_states) {
        elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at_time, units = "secs"))
        remaining_timeout <- max(1, floor(job_timeout_seconds - elapsed_seconds))
        full_result <- get_full_results_response(result_id, timeout_seconds = remaining_timeout)
        result_json <- full_result$text
        log_progress("[%d/%d] FULL_RESULT downloaded result_id=%s http_status=%s", idx, total_items, result_id, full_result$http_status)
      } else if (final_status %in% terminal_failure_states) {
        elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at_time, units = "secs"))
        remaining_timeout <- max(1, floor(job_timeout_seconds - elapsed_seconds))
        err_result <- tryCatch(get_full_results_response(result_id, timeout_seconds = remaining_timeout), error = function(e) {
          log_progress("[%d/%d] ERROR_FETCH_FAILED result_id=%s reason=%s", idx, total_items, result_id, conditionMessage(e))
          NULL
        })
        if (!is.null(err_result)) {
          result_json <- err_result$text
          log_progress("[%d/%d] ERROR_BODY http=%s result_id=%s body=%s", idx, total_items, err_result$http_status, result_id,
                       substr(err_result$text, 1, 300))
          parsed_err <- tryCatch(fromJSON(err_result$text, simplifyVector = FALSE), error = function(e) NULL)
          if (!is.null(parsed_err)) {
            msgs <- c(
              parsed_err$errors,
              parsed_err$error,
              parsed_err$message,
              parsed_err$data$errors,
              parsed_err$data$error
            )
            msgs <- msgs[!sapply(msgs, is.null)]
            if (length(msgs) > 0) {
              error_message <- paste(unlist(msgs), collapse = "; ")
            }
          }
          log_progress("[%d/%d] ERROR_DETAIL result_id=%s message=%s", idx, total_items, result_id, error_message)
        }
      }
    }, error = function(e) {
      run_error <<- conditionMessage(e)
    })

    if (!is.null(run_error)) {
      error_message <- run_error
      log_progress("[%d/%d] ERROR HLP=%s allele=%s region=%s message=%s", idx, total_items, item$hlp, item$allele, item$protein_region, run_error)
      if (is.na(final_status) || !nzchar(final_status)) {
        final_status <- "error"
      }
    }

    # break on success or if this was the last retry attempt
    if (!is.null(final_status) && !is.na(final_status) && final_status %in% terminal_success_states) {
      break
    }
    if (is_last_attempt) {
      break
    }
    # only retry on server-side failures, not on submission-level problems
    if (is.null(final_status) || is.na(final_status) || !(final_status %in% terminal_failure_states)) {
      break
    }
  }

  log_progress("[%d/%d] FINISH HLP=%s allele=%s region=%s final_status=%s", idx, total_items, item$hlp, item$allele, item$protein_region, final_status)

  list(
    idx = idx,
    total = total_items,
    hlp = item$hlp,
    protein_region = item$protein_region,
    allele = item$allele,
    request_payload_json = payload_json,
    submission_http_status = submission_http_status,
    submission_response_json = submission_response_json,
    result_id = result_id,
    final_status = final_status,
    result_json = result_json,
    error_message = error_message,
    started_at = started_at,
    completed_at = now_utc()
  )
}

total_items <- length(work_items)
run_started <- Sys.time()
processed_count <- 0L
ok_count <- 0L
error_count <- 0L
dry_run_count <- 0L

save_record <- function(record) {
  dbExecute(
    con,
    update_sql,
    params = list(
      record$request_payload_json,
      record$submission_http_status,
      record$submission_response_json,
      record$result_id,
      record$final_status,
      record$result_json,
      record$error_message,
      record$started_at,
      record$completed_at,
      record$hlp,
      record$protein_region,
      record$allele
    )
  )

  processed_count <<- processed_count + 1L
  if (identical(record$final_status, "dry_run")) {
    dry_run_count <<- dry_run_count + 1L
  } else if (!is.na(record$error_message) && nzchar(record$error_message)) {
    error_count <<- error_count + 1L
  } else {
    ok_count <<- ok_count + 1L
  }

  if (record$final_status == "dry_run") {
    cat(sprintf("[%d/%d] DRY RUN saved: HLP %s, allele %s, region %s\n", record$idx, record$total, record$hlp, record$allele, record$protein_region))
  } else if (!is.na(record$error_message) && nzchar(record$error_message)) {
    cat(sprintf("[%d/%d] Saved with error: HLP %s allele %s region %s -> %s\n", record$idx, record$total, record$hlp, record$allele, record$protein_region, record$error_message))
  } else {
    cat(sprintf("[%d/%d] Saved final_status=%s result_id=%s (HLP %s allele %s region %s)\n", record$idx, record$total, record$final_status, ifelse(is.na(record$result_id), "NA", record$result_id), record$hlp, record$allele, record$protein_region))
  }

  elapsed <- as.integer(difftime(Sys.time(), run_started, units = "secs"))
  cat(sprintf("GLOBAL progress: completed=%d/%d ok=%d errors=%d dry_run=%d elapsed=%ds\n", processed_count, total_items, ok_count, error_count, dry_run_count, elapsed))
}

if (workers > 1L) {
  cat(sprintf("Executing API calls in parallel with %d workers...\n", workers))
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)

  submit_one <- function(i) {
    future(
      process_work_item(work_items[[i]], i, total_items, dry_run),
      packages = c("jsonlite", "httr2")
    )
  }

  next_index <- 1L
  active <- list()

  while (next_index <= total_items && length(active) < workers) {
    cat(sprintf("Dispatching item %d of %d\n", next_index, total_items))
    active[[length(active) + 1L]] <- list(index = next_index, fut = submit_one(next_index))
    next_index <- next_index + 1L
  }

  while (length(active) > 0L) {
    ready <- vapply(active, function(x) resolved(x$fut), logical(1))

    if (!any(ready)) {
      Sys.sleep(0.25)
      next
    }

    ready_positions <- which(ready)

    for (pos in rev(ready_positions)) {
      item_index <- active[[pos]]$index
      record <- tryCatch(
        value(active[[pos]]$fut),
        error = function(e) {
          list(
            idx = item_index,
            total = total_items,
            hlp = work_items[[item_index]]$hlp,
            protein_region = work_items[[item_index]]$protein_region,
            allele = work_items[[item_index]]$allele,
            request_payload_json = "{}",
            submission_http_status = NA_integer_,
            submission_response_json = NA_character_,
            result_id = NA_character_,
            final_status = "error",
            result_json = NA_character_,
            error_message = sprintf("Future failed before result return: %s", conditionMessage(e)),
            started_at = now_utc(),
            completed_at = now_utc()
          )
        }
      )

      save_record(record)
      active[[pos]] <- NULL

      if (next_index <= total_items) {
        cat(sprintf("Dispatching item %d of %d\n", next_index, total_items))
        active[[length(active) + 1L]] <- list(index = next_index, fut = submit_one(next_index))
        next_index <- next_index + 1L
      }
    }
  }
} else {
  processed_items <- lapply(
    seq_along(work_items),
    function(i) process_work_item(work_items[[i]], i, total_items, dry_run)
  )

  for (record in processed_items) {
    save_record(record)
  }
}

summary_df <- dbGetQuery(con, "
  SELECT final_status, COUNT(*) AS n
  FROM hlp_donor_pipeline_results
  GROUP BY final_status
  ORDER BY n DESC
")

cat("\nRun summary by final_status:\n")
print(summary_df)
