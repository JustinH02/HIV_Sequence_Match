args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("Usage: Rscript main.R <dna_file_1> <dna_file_2> [k]\nExample: Rscript main.R A04321.raw B04321.raw 9")
}

file1 <- args[1]
file2 <- args[2]
k <- if (length(args) >= 3) as.integer(args[3]) else 9L

if (!file.exists(file1)) {
	stop(paste("Input file not found:", file1))
}

if (!file.exists(file2)) {
	stop(paste("Input file not found:", file2))
}

if (is.na(k) || k <= 0) {
	stop("k must be a positive integer.")
}

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

dna1 <- read_dna(file1)
dna2 <- read_dna(file2)

kmers1 <- extract_kmers(dna1, k)
kmers2 <- extract_kmers(dna2, k)

if (length(kmers1) == 0 || length(kmers2) == 0) {
	cat(sprintf("No %d-letter sequences can be formed from one or both files.\n", k))
	quit(save = "no", status = 0)
}

matches <- sort(unique(intersect(kmers1, kmers2)))

if (length(matches) == 0) {
	cat(sprintf("No matching %d-letter sequences found.\n", k))
} else {
	positions1 <- split(seq_along(kmers1), kmers1)
	positions2 <- split(seq_along(kmers2), kmers2)

	rows <- lapply(matches, function(seq_match) {
		pos1 <- positions1[[seq_match]]
		pos2 <- positions2[[seq_match]]

		# Build all index pairs where the same k-mer appears in both strands.
		grid <- expand.grid(start_in_file1 = pos1, start_in_file2 = pos2)
		data.frame(
			start_in_file1 = grid$start_in_file1,
			start_in_file2 = grid$start_in_file2,
			match = seq_match,
			stringsAsFactors = FALSE
		)
	})

	result <- do.call(rbind, rows)
	rownames(result) <- NULL

	cat(sprintf("Found %d matching %d-letter sequence occurrence pairs:\n", nrow(result), k))
	write.table(result, row.names = FALSE, quote = FALSE, sep = "\t")
}
