# R script to comprehensively analyze class consistency in JSON files
# Excludes lulc files as they have a different structure

library(jsonlite)

# Define the directory path
perc_area_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/chart_data/perc_area"

# Get all JSON files, excluding those with "lulc" in the name
json_files <- list.files(perc_area_dir, pattern = "*.json", full.names = TRUE)
json_files <- json_files[!grepl("lulc", basename(json_files))]

cat("Found", length(json_files), "JSON files (excluding lulc files)\n\n")

# Function to extract class values from a JSON file
extract_classes <- function(file_path) {
  tryCatch(
    {
      data <- fromJSON(file_path, simplifyVector = FALSE)

      # Extract class values (the array values in the JSON)
      classes <- character(0)
      for (key in names(data)) {
        if (is.list(data[[key]]) && length(data[[key]]) > 0) {
          classes <- c(classes, data[[key]][[1]])
        }
      }

      # Sort classes for comparison
      sorted_classes <- sort(classes)

      return(list(
        file = basename(file_path),
        classes = sorted_classes,
        success = TRUE
      ))
    },
    error = function(e) {
      return(list(
        file = basename(file_path),
        classes = character(0),
        success = FALSE,
        error = e$message
      ))
    }
  )
}

# Process all files
cat("Processing all files...\n")
all_results <- lapply(json_files, extract_classes)

# Check for failures
failed_files <- sapply(all_results, function(x) !x$success)
success_count <- sum(!failed_files)

cat("Successfully processed:", success_count, "files\n")
cat("Failed to process:", sum(failed_files), "files\n\n")

if (sum(failed_files) > 0) {
  cat("Failed files:\n")
  failed_info <- all_results[failed_files]
  for (failed in failed_info[1:min(5, length(failed_info))]) {
    cat("  ", failed$file, "\n")
  }
  if (sum(failed_files) > 5) {
    cat("  ... and", sum(failed_files) - 5, "more\n")
  }
  cat("\n")
}

# Analyze successful extractions
successful_results <- all_results[!failed_files]

if (length(successful_results) > 0) {
  # Group files by their class patterns
  class_patterns <- list()

  for (result in successful_results) {
    class_key <- paste(result$classes, collapse = "|")

    if (!(class_key %in% names(class_patterns))) {
      class_patterns[[class_key]] <- list(
        classes = result$classes,
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

  cat("=== ANALYSIS RESULTS ===\n\n")

  if (length(class_patterns) == 1) {
    cat("✅ ALL FILES USE CONSISTENT CLASS VALUES!\n\n")
    cat("Standard class values used across all", success_count, "files:\n")
    standard_classes <- class_patterns[[1]]$classes
    for (i in seq_along(standard_classes)) {
      cat("  ", i, ". ", standard_classes[i], "\n", sep = "")
    }
  } else {
    cat("❌ INCONSISTENT CLASS VALUES FOUND!\n\n")
    cat("Found", length(class_patterns), "different class patterns:\n\n")

    for (i in seq_along(class_patterns)) {
      pattern <- class_patterns[[i]]
      file_count <- length(pattern$files)
      percentage <- round(file_count / success_count * 100, 1)

      cat("--- Pattern", i, "(", file_count, "files,", percentage, "%) ---\n")
      cat("Classes:\n")
      for (j in seq_along(pattern$classes)) {
        cat("  ", j, ". ", pattern$classes[j], "\n", sep = "")
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
    }

    # Identify the most common pattern as the "standard"
    most_common_pattern <- class_patterns[[1]]
    cat(
      "Most common pattern (Pattern 1) appears to be the intended standard.\n"
    )
    cat("Files NOT following the most common pattern:\n")

    non_standard_count <- 0
    for (i in 2:length(class_patterns)) {
      pattern <- class_patterns[[i]]
      non_standard_count <- non_standard_count + length(pattern$files)
      cat("\n Pattern", i, "files (", length(pattern$files), "):\n")
      for (file in pattern$files) {
        cat("  - ", file, "\n", sep = "")
      }
    }

    cat("\nSUMMARY:\n")
    cat("- Total files analyzed:", success_count, "\n")
    cat(
      "- Files following most common pattern:",
      length(most_common_pattern$files),
      "(",
      round(length(most_common_pattern$files) / success_count * 100, 1),
      "%)\n"
    )
    cat(
      "- Files with inconsistent patterns:",
      non_standard_count,
      "(",
      round(non_standard_count / success_count * 100, 1),
      "%)\n"
    )
  }
}

cat("\nAnalysis completed.\n")
