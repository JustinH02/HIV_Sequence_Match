library(shiny)

read_dna <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dna <- toupper(paste(lines, collapse = ""))
  dna <- gsub("[^ACGTN]", "", dna)

  if (nchar(dna) == 0) {
    stop(paste("No valid DNA bases (A/C/G/T/N) found in file:", path))
  }

  dna
}

extract_kmers <- function(sequence, k) {
  n <- nchar(sequence)
  if (n < k) {
    return(character(0))
  }

  starts <- seq_len(n - k + 1)
  substring(sequence, starts, starts + k - 1)
}

build_match_summary <- function(kmers1, kmers2) {
  matches <- sort(unique(intersect(kmers1, kmers2)))
  if (length(matches) == 0) {
    return(list(summary = data.frame(), positions1 = list(), positions2 = list()))
  }

  positions1 <- split(seq_along(kmers1), kmers1)
  positions2 <- split(seq_along(kmers2), kmers2)

  collapse_positions <- function(pos) {
    paste(pos, collapse = ", ")
  }

  summary_rows <- lapply(matches, function(seq_match) {
    pos1 <- positions1[[seq_match]]
    pos2 <- positions2[[seq_match]]

    data.frame(
      match = seq_match,
      occurrences_in_file1 = length(pos1),
      occurrences_in_file2 = length(pos2),
      total_position_pairs = as.double(length(pos1)) * as.double(length(pos2)),
      starts_in_file1 = collapse_positions(pos1),
      starts_in_file2 = collapse_positions(pos2),
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  summary <- summary[order(-summary$total_position_pairs, summary$match), ]

  list(summary = summary, positions1 = positions1, positions2 = positions2)
}

build_pairs_for_match <- function(seq_match, positions1, positions2, max_rows = 2000L) {
  pos1 <- positions1[[seq_match]]
  pos2 <- positions2[[seq_match]]

  if (is.null(pos1) || is.null(pos2)) {
    return(data.frame())
  }

  pairs <- expand.grid(start_in_file1 = pos1, start_in_file2 = pos2)
  pairs$match <- seq_match

  if (!is.na(max_rows) && nrow(pairs) > max_rows) {
    pairs <- pairs[seq_len(max_rows), ]
  }

  pairs
}

ui <- fluidPage(
  titlePanel("DNA Match Explorer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "DNA File 1", accept = c(".raw", ".txt", ".fasta", ".fa")),
      fileInput("file2", "DNA File 2", accept = c(".raw", ".txt", ".fasta", ".fa")),
      numericInput("k", "Sequence Length (k)", value = 9, min = 1, step = 1),
      numericInput("max_rows", "Max detail rows to display", value = 2000, min = 50, step = 50),
      actionButton("run", "Find Matches", class = "btn-primary")
    ),
    mainPanel(
      h4("Run Status"),
      verbatimTextOutput("status"),
      h4("Matching Sequence Summary"),
      tableOutput("summary_table"),
      h4("Position Pairs For Selected Match"),
      selectInput("selected_match", "Choose a matching sequence", choices = character(0)),
      tableOutput("pairs_table")
    )
  )
)

server <- function(input, output, session) {
  run_result <- eventReactive(input$run, {
    req(input$file1$datapath, input$file2$datapath)

    k <- as.integer(input$k)
    validate(
      need(!is.na(k) && k > 0, "k must be a positive integer.")
    )

    dna1 <- read_dna(input$file1$datapath)
    dna2 <- read_dna(input$file2$datapath)

    kmers1 <- extract_kmers(dna1, k)
    kmers2 <- extract_kmers(dna2, k)

    if (length(kmers1) == 0 || length(kmers2) == 0) {
      return(list(
        status = sprintf("No %d-letter sequences can be formed from one or both files.", k),
        summary = data.frame(),
        positions1 = list(),
        positions2 = list(),
        k = k
      ))
    }

    built <- build_match_summary(kmers1, kmers2)

    status <- if (nrow(built$summary) == 0) {
      sprintf("No matching %d-letter sequences found.", k)
    } else {
      sprintf(
        "Found %d unique matching %d-letter sequences.",
        nrow(built$summary),
        k
      )
    }

    list(
      status = status,
      summary = built$summary,
      positions1 = built$positions1,
      positions2 = built$positions2,
      k = k
    )
  })

  observeEvent(run_result(), {
    result <- run_result()
    choices <- result$summary$match
    updateSelectInput(
      session,
      "selected_match",
      choices = choices,
      selected = if (length(choices) > 0) choices[1] else character(0)
    )
  })

  output$status <- renderText({
    req(run_result())
    run_result()$status
  })

  output$summary_table <- renderTable({
    req(run_result())
    run_result()$summary
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$pairs_table <- renderTable({
    req(run_result())
    req(input$selected_match)

    result <- run_result()
    if (nrow(result$summary) == 0) {
      return(data.frame())
    }

    build_pairs_for_match(
      seq_match = input$selected_match,
      positions1 = result$positions1,
      positions2 = result$positions2,
      max_rows = as.integer(input$max_rows)
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shinyApp(ui = ui, server = server)
