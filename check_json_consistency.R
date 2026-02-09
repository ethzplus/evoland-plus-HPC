# Check consistency of class values across JSON files in perc_area directory
# Also checks that percentage values don't exceed 100% in any single file
# Excludes lulc files as they have a different structure

library(jsonlite)
library(dplyr)

# Source config file to get directory path
config_file <- "src/config.yml"
config <- yaml::read_yaml(config_file)

# Define the directory path from config
perc_area_dir <- file.path(
  config$Summarisation$OutputDir,
  config$Summarisation$PercAreaDir
)

cat("Checking JSON files in directory:", perc_area_dir, "\n\n")

# Get all JSON files, excluding those with "lulc" in the name
json_files <- list.files(perc_area_dir, pattern = "*.json", full.names = TRUE)
json_files_filtered <- json_files[
  !grepl("lulc", basename(json_files), ignore.case = TRUE)
]

cat("Total JSON files found:", length(json_files), "\n")
cat("Non-lulc JSON files to analyze:", length(json_files_filtered), "\n\n")

if (length(json_files_filtered) == 0) {
  stop("No non-lulc JSON files found in the directory!")
}

# Function to extract class values and check percentage sums from a JSON file
extract_classes <- function(file_path) {
  tryCatch(
    {
      # Read JSON file
      data <- fromJSON(file_path, simplifyVector = FALSE)

      # Extract class values (the array values in the JSON)
      # The structure should be: percentage -> [class_label]
      classes <- character(0)
      percentages <- numeric(0)

      for (percentage in names(data)) {
        if (is.list(data[[percentage]]) && length(data[[percentage]]) > 0) {
          classes <- c(classes, data[[percentage]][[1]])
          # Convert percentage key to numeric value
          percentages <- c(percentages, as.numeric(percentage))
        }
      }

      # Sort classes for comparison
      sorted_classes <- sort(classes)

      # Calculate total percentage
      total_percentage <- sum(percentages)
      exceeds_100 <- total_percentage > 100.01 # Allow small floating point tolerance

      return(list(
        file = basename(file_path),
        classes = sorted_classes,
        n_classes = length(sorted_classes),
        total_percentage = total_percentage,
        exceeds_100 = exceeds_100,
        success = TRUE,
        error = NULL
      ))
    },
    error = function(e) {
      return(list(
        file = basename(file_path),
        classes = character(0),
        n_classes = 0,
        total_percentage = NA,
        exceeds_100 = NA,
        success = FALSE,
        error = e$message
      ))
    }
  )
}

# Process files in batches to avoid memory issues
batch_size <- 1000
n_batches <- ceiling(length(json_files_filtered) / batch_size)

cat(
  "Processing files in",
  n_batches,
  "batches of",
  batch_size,
  "files each...\n"
)

all_results <- list()
for (batch in 1:n_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(json_files_filtered))

  cat(
    "Processing batch",
    batch,
    "of",
    n_batches,
    "(files",
    start_idx,
    "to",
    end_idx,
    ")...\n"
  )

  batch_files <- json_files_filtered[start_idx:end_idx]
  batch_results <- lapply(batch_files, extract_classes)

  all_results <- c(all_results, batch_results)

  # Clean up memory
  gc()
}

cat("\nAnalysis complete! Processing results...\n\n")

# Check for failures
failed_files <- sapply(all_results, function(x) !x$success)
success_count <- sum(!failed_files)
failed_count <- sum(failed_files)

cat("=== PROCESSING SUMMARY ===\n")
cat("Successfully processed:", success_count, "files\n")
cat("Failed to process:", failed_count, "files\n\n")

# Check for files exceeding 100%
if (success_count > 0) {
  successful_results <- all_results[!failed_files]
  exceeds_100_files <- sapply(successful_results, function(x) x$exceeds_100)
  exceeds_count <- sum(exceeds_100_files, na.rm = TRUE)

  cat("=== PERCENTAGE SUM CHECK ===\n")
  if (exceeds_count == 0) {
    cat("✅ ALL FILES HAVE PERCENTAGES SUMMING TO ≤100%!\n\n")

    # Show distribution of percentage totals
    all_totals <- unlist(sapply(successful_results, function(x) {
      x$total_percentage
    }))
    cat("Distribution of percentage totals:\n")
    cat("  Min:", round(min(all_totals, na.rm = TRUE), 2), "%\n")
    cat("  Max:", round(max(all_totals, na.rm = TRUE), 2), "%\n")
    cat("  Mean:", round(mean(all_totals, na.rm = TRUE), 2), "%\n")
    cat("  Median:", round(median(all_totals, na.rm = TRUE), 2), "%\n\n")

    # Count files at exactly 100%
    exactly_100 <- sum(abs(all_totals - 100) < 0.01, na.rm = TRUE)
    cat(
      "Files with totals at exactly 100%:",
      exactly_100,
      "(",
      round(exactly_100 / success_count * 100, 1),
      "%)\n\n"
    )
  } else {
    cat("❌ FOUND", exceeds_count, "FILES WITH PERCENTAGES EXCEEDING 100%!\n\n")

    exceeds_info <- successful_results[exceeds_100_files]

    # Show first 20 problematic files
    cat("Files exceeding 100% (first 20):\n")
    for (i in 1:min(20, length(exceeds_info))) {
      info <- exceeds_info[[i]]
      cat(sprintf("  %-60s Total: %.2f%%\n", info$file, info$total_percentage))
    }
    if (exceeds_count > 20) {
      cat("  ... and", exceeds_count - 20, "more files\n")
    }
    cat("\n")

    # Distribution of percentage totals for exceeding files
    exceeds_totals <- unlist(sapply(exceeds_info, function(x) {
      x$total_percentage
    }))
    cat("Distribution of totals for files exceeding 100%:\n")
    cat("  Min:", round(min(exceeds_totals, na.rm = TRUE), 2), "%\n")
    cat("  Max:", round(max(exceeds_totals, na.rm = TRUE), 2), "%\n")
    cat("  Mean:", round(mean(exceeds_totals, na.rm = TRUE), 2), "%\n")
    cat("  Median:", round(median(exceeds_totals, na.rm = TRUE), 2), "%\n\n")

    # Check if specific ES types are affected
    cat("Files exceeding 100% by ES type:\n")
    es_types <- c("ndr", "pol", "ff", "wy", "car", "rec", "sdr", "hab")
    for (es_type in es_types) {
      es_exceeds <- grep(
        paste0("^", es_type, "-"),
        sapply(exceeds_info, function(x) x$file),
        value = FALSE
      )
      if (length(es_exceeds) > 0) {
        cat("  ", toupper(es_type), ":", length(es_exceeds), "files\n")
      }
    }
    cat("\n")
  }
}


if (failed_count > 0) {
  cat("Failed files (first 10):\n")
  failed_info <- all_results[failed_files]
  for (i in 1:min(10, length(failed_info))) {
    cat("  ", failed_info[[i]]$file, ": ", failed_info[[i]]$error, "\n")
  }
  if (failed_count > 10) {
    cat("  ... and", failed_count - 10, "more\n")
  }
  cat("\n")
}

# Analyze successful extractions
if (success_count == 0) {
  stop("No files were successfully processed!")
}

successful_results <- all_results[!failed_files]

# Group files by their class patterns
class_patterns <- list()

cat("Grouping files by class patterns...\n")
for (i in seq_along(successful_results)) {
  if (i %% 1000 == 0) {
    cat("  Processed", i, "of", length(successful_results), "files\n")
  }

  result <- successful_results[[i]]
  class_key <- paste(result$classes, collapse = "|")

  if (!(class_key %in% names(class_patterns))) {
    class_patterns[[class_key]] <- list(
      classes = result$classes,
      n_classes = result$n_classes,
      files = character(0)
    )
  }

  class_patterns[[class_key]]$files <- c(
    class_patterns[[class_key]]$files,
    result$file
  )
}

# Sort patterns by number of files (most common first)
pattern_counts <- sapply(class_patterns, function(x) length(x$files))
class_patterns <- class_patterns[order(pattern_counts, decreasing = TRUE)]

cat("\n=== CONSISTENCY ANALYSIS RESULTS ===\n\n")

if (length(class_patterns) == 1) {
  cat("✅ ALL FILES USE CONSISTENT CLASS VALUES!\n\n")
  cat("Standard class values used across all", success_count, "files:\n")
  standard_classes <- class_patterns[[1]]$classes
  for (i in seq_along(standard_classes)) {
    cat("  ", i, ". '", standard_classes[i], "'\n", sep = "")
  }
  cat("\nTotal classes per file:", class_patterns[[1]]$n_classes, "\n")
} else {
  cat("❌ INCONSISTENT CLASS VALUES FOUND!\n\n")
  cat("Found", length(class_patterns), "different class patterns:\n\n")

  # Summary table
  summary_df <- data.frame(
    Pattern = 1:length(class_patterns),
    Files = sapply(class_patterns, function(x) length(x$files)),
    Percentage = round(
      sapply(class_patterns, function(x) length(x$files)) / success_count * 100,
      1
    ),
    Classes = sapply(class_patterns, function(x) x$n_classes),
    stringsAsFactors = FALSE
  )

  print(summary_df)
  cat("\n")

  # Detailed breakdown of each pattern
  for (i in seq_along(class_patterns)) {
    pattern <- class_patterns[[i]]
    file_count <- length(pattern$files)
    percentage <- round(file_count / success_count * 100, 1)

    cat("--- Pattern", i, "(", file_count, "files,", percentage, "%) ---\n")
    cat("Number of classes:", pattern$n_classes, "\n")
    cat("Classes:\n")
    for (j in seq_along(pattern$classes)) {
      cat("  ", j, ". '", pattern$classes[j], "'\n", sep = "")
    }

    cat("\nSample files (first 5):\n")
    sample_files <- head(pattern$files, 5)
    for (file in sample_files) {
      cat("  - ", file, "\n", sep = "")
    }
    if (file_count > 5) {
      cat("  ... and", file_count - 5, "more files\n")
    }
    cat("\n")

    # Stop after showing first 5 patterns to avoid overwhelming output
    if (i >= 5 && length(class_patterns) > 5) {
      cat("... and", length(class_patterns) - 5, "more patterns\n\n")
      break
    }
  }

  # Identify files that don't match the most common pattern
  most_common_pattern <- class_patterns[[1]]
  non_standard_files <- character(0)

  for (i in 2:length(class_patterns)) {
    pattern <- class_patterns[[i]]
    non_standard_files <- c(non_standard_files, pattern$files)
  }

  cat("=== SUMMARY ===\n")
  cat("- Total files analyzed:", success_count, "\n")
  cat(
    "- Files following most common pattern:",
    length(most_common_pattern$files),
    "(",
    round(length(most_common_pattern$files) / success_count * 100, 1),
    "%)\n"
  )
  cat(
    "- Files with non-standard patterns:",
    length(non_standard_files),
    "(",
    round(length(non_standard_files) / success_count * 100, 1),
    "%)\n"
  )

  # Extract ES types from filenames for pattern analysis
  cat("\n=== PATTERN BY ES TYPE ===\n")
  es_types <- c("ndr", "pol", "ff", "wy", "car", "rec", "sdr", "hab")

  for (es_type in es_types) {
    es_files <- grep(
      paste0("^", es_type, "-"),
      sapply(successful_results, function(x) x$file),
      value = TRUE
    )
    if (length(es_files) > 0) {
      # Find which patterns these files belong to
      es_patterns <- character(0)
      for (i in seq_along(class_patterns)) {
        pattern_files <- class_patterns[[i]]$files
        es_in_pattern <- sum(es_files %in% pattern_files)
        if (es_in_pattern > 0) {
          es_patterns <- c(
            es_patterns,
            paste0("Pattern ", i, " (", es_in_pattern, " files)")
          )
        }
      }
      cat(
        "- ",
        toupper(es_type),
        ": ",
        length(es_files),
        " files -> ",
        paste(es_patterns, collapse = ", "),
        "\n"
      )
    }
  }
}

# Save detailed results
output_file <- "json_consistency_analysis_results.rds"
saveRDS(
  list(
    summary = list(
      total_files = length(json_files_filtered),
      successful = success_count,
      failed = failed_count,
      n_patterns = length(class_patterns),
      exceeds_100_count = if (success_count > 0) exceeds_count else NA
    ),
    class_patterns = class_patterns,
    failed_files = if (failed_count > 0) all_results[failed_files] else NULL,
    exceeds_100_files = if (success_count > 0 && exceeds_count > 0) {
      successful_results[exceeds_100_files]
    } else {
      NULL
    }
  ),
  output_file
)

cat("\nDetailed results saved to:", output_file, "\n")
cat("\nAnalysis completed!\n")
