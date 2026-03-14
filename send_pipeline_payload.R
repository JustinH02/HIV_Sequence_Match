if (!requireNamespace("httr2", quietly = TRUE)) {
  stop("Package 'httr2' is required. Install with: install.packages('httr2')")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
}

library(httr2)
library(jsonlite)

api_url <- "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
results_base_url <- "https://api-nextgen-tools.iedb.org/api/v1/results"
poll_interval_seconds <- 10
max_poll_attempts <- 30
output_dir <- "."

pretty_print_json <- function(text) {
  parsed <- tryCatch(fromJSON(text), error = function(e) NULL)
  if (is.null(parsed)) {
    cat(text, "\n")
  } else {
    cat(toJSON(parsed, pretty = TRUE, auto_unbox = TRUE), "\n")
  }
}

save_json_text <- function(json_text, file_path) {
  parsed <- tryCatch(fromJSON(json_text, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(parsed)) {
    writeLines(json_text, con = file_path, useBytes = TRUE)
  } else {
    formatted <- toJSON(parsed, pretty = TRUE, auto_unbox = TRUE)
    writeLines(formatted, con = file_path, useBytes = TRUE)
  }
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

  if (is.atomic(x) && length(x) == 1 && !is.na(x)) {
    return(NULL)
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

get_status_response <- function(results_id) {
  status_url <- sprintf("%s/%s?statusOnly=true", results_base_url, URLencode(results_id, reserved = TRUE))
  status_req <- request(status_url) |>
    req_method("GET") |>
    req_headers("Accept" = "application/json")

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

get_full_results_response <- function(results_id) {
  full_url <- sprintf("%s/%s", results_base_url, URLencode(results_id, reserved = TRUE))
  full_req <- request(full_url) |>
    req_method("GET") |>
    req_headers("Accept" = "application/json")

  full_resp <- req_perform(full_req)

  list(
    http_status = resp_status(full_resp),
    text = resp_body_string(full_resp)
  )
}

payload <- list(
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
      input_sequence_text = paste(
        ">Gag",
        "MGARASVLSGGQLDRWEKIRLRPGGKKKYQLKHVVWASRELERFAVNPGLLETSGGC",
        "KQILEQLQPSLQTGSEELKSLFNTVAVLYCVHQKIEVKDTKEALDKIEEERNKSK",
        "KMAQQAAAGTGNSSQVSQNYPIVQNLQGQMVHQAISPRTLNAWVKVVEEKAF",
        "SPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHP",
        "VHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNK",
        "IVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNA",
        "NPDCKTILKALGPAATLEEMMTACQGVGGPGHKAKILAEAMSQVTSSATIMMQR",
        "GNFKNQRPIKCFNCGKVGHLAKHCRAPRKRGCWKCGKEGHQMKDCTERQAN",
        "FLGKIWPSNKGRPGNFLQSRPQPTAPPAPPEESFRFGEETTTPPQEQNQIDKEL",
        "YPLTSLKSLFGNDPSSQ",
        sep = "\n"
      ),
      input_parameters = list(
        alleles = "HLA-DRB1*11:01",
        peptide_shift = 5,
        peptide_length_range = list(15, 15),
        predictors = list(
          list(
            type = "binding",
            method = "netmhciipan_el"
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

req <- request(api_url) |>
  req_method("POST") |>
  req_headers(
    "Content-Type" = "application/json",
    "Accept" = "application/json"
  ) |>
  req_body_json(payload, auto_unbox = TRUE)

resp <- req_perform(req)

cat("Status:", resp_status(resp), "\n\n")

resp_text <- resp_body_string(resp)
pretty_print_json(resp_text)

parsed_raw <- tryCatch(fromJSON(resp_text, simplifyVector = FALSE), error = function(e) NULL)
results_id <- if (is.null(parsed_raw)) NULL else find_results_id(parsed_raw)

if (is.null(results_id) || !nzchar(results_id)) {
  stop("Could not find results ID in first API response.")
}

cat("\nResults ID:", results_id, "\n\n")

pending_states <- c("pending", "queued", "running", "in_progress", "processing")
terminal_success_states <- c("completed", "complete", "succeeded", "success", "done")
terminal_failure_states <- c("failed", "error", "cancelled", "canceled")

final_status <- NULL
for (attempt in seq_len(max_poll_attempts)) {
  status_result <- get_status_response(results_id)
  status_value <- tolower(if (is.null(status_result$status)) "unknown" else status_result$status)

  cat(sprintf("Poll %d/%d - HTTP %s - status: %s\n", attempt, max_poll_attempts, status_result$http_status, status_value))
  pretty_print_json(status_result$text)
  cat("\n")

  if (status_value %in% terminal_success_states) {
    final_status <- status_value
    break
  }

  if (status_value %in% terminal_failure_states) {
    stop(sprintf("Pipeline result reached failure state: %s", status_value))
  }

  if (!(status_value %in% pending_states) && status_value != "unknown") {
    final_status <- status_value
    break
  }

  if (attempt < max_poll_attempts) {
    Sys.sleep(poll_interval_seconds)
  }
}

if (is.null(final_status)) {
  stop(sprintf("Polling timeout after %d attempts (~%d seconds).", max_poll_attempts, max_poll_attempts * poll_interval_seconds))
}

cat(sprintf("Final observed status: %s\n", final_status))

if (final_status %in% terminal_success_states) {
  cat("\nFetching full results payload...\n\n")
  full_result <- get_full_results_response(results_id)
  cat("Full Results Endpoint HTTP Status:", full_result$http_status, "\n\n")
  pretty_print_json(full_result$text)

  output_file <- file.path(output_dir, sprintf("results_%s.json", results_id))
  save_json_text(full_result$text, output_file)
  cat("\nSaved full results to:", normalizePath(output_file, winslash = "/", mustWork = FALSE), "\n")
}
