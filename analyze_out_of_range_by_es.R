# Analyze which ES types have files with "Bin_0_OutOfRange" class
# Using the json_check object loaded from the RDS file

# Extract files that contain "Bin_0_OutOfRange" class
files_with_out_of_range <- character(0)

# Loop through all class patterns to find those containing "Bin_0_OutOfRange"
for (pattern_key in names(json_check$class_patterns)) {
  pattern <- json_check$class_patterns[[pattern_key]]

  # Check if this pattern contains "Bin_0_OutOfRange"
  if ("Bin_0_OutOfRange" %in% pattern$classes) {
    files_with_out_of_range <- c(files_with_out_of_range, pattern$files)
  }
}

cat(
  "Total files with Bin_0_OutOfRange:",
  length(files_with_out_of_range),
  "\n\n"
)

# Extract ES types from filenames
# ES types are typically the first part of the filename before the first hyphen
es_types <- c("ndr", "pol", "ff", "wy", "car", "rec", "sdr", "hab")

# Initialize counters
es_analysis <- data.frame(
  ES_Type = character(0),
  Files_With_OutOfRange = integer(0),
  Total_Files = integer(0),
  Percentage = numeric(0),
  stringsAsFactors = FALSE
)

# Get all successful files for comparison
all_successful_files <- character(0)
for (pattern in json_check$class_patterns) {
  all_successful_files <- c(all_successful_files, pattern$files)
}

# Analyze each ES type
for (es_type in es_types) {
  # Find all files for this ES type
  es_files <- grep(
    paste0("^", es_type, "-"),
    all_successful_files,
    value = TRUE
  )

  # Find files with out-of-range values for this ES type
  es_files_with_out_of_range <- grep(
    paste0("^", es_type, "-"),
    files_with_out_of_range,
    value = TRUE
  )

  if (length(es_files) > 0) {
    percentage <- round(
      length(es_files_with_out_of_range) / length(es_files) * 100,
      1
    )

    es_analysis <- rbind(
      es_analysis,
      data.frame(
        ES_Type = toupper(es_type),
        Files_With_OutOfRange = length(es_files_with_out_of_range),
        Total_Files = length(es_files),
        Percentage = percentage,
        stringsAsFactors = FALSE
      )
    )
  }
}

# Sort by percentage (highest first)
es_analysis <- es_analysis[order(es_analysis$Percentage, decreasing = TRUE), ]

cat("=== ES ANALYSIS: Files with Bin_0_OutOfRange ===\n")
print(es_analysis)
cat("\n")

# Summary statistics
cat("=== SUMMARY ===\n")
cat("ES types with ALL files having out-of-range values:\n")
all_affected <- es_analysis[es_analysis$Percentage == 100, ]
if (nrow(all_affected) > 0) {
  for (i in 1:nrow(all_affected)) {
    cat(
      "  -",
      all_affected$ES_Type[i],
      "(",
      all_affected$Total_Files[i],
      "files )\n"
    )
  }
} else {
  cat("  None\n")
}

cat("\nES types with NO files having out-of-range values:\n")
no_affected <- es_analysis[es_analysis$Percentage == 0, ]
if (nrow(no_affected) > 0) {
  for (i in 1:nrow(no_affected)) {
    cat(
      "  -",
      no_affected$ES_Type[i],
      "(",
      no_affected$Total_Files[i],
      "files )\n"
    )
  }
} else {
  cat("  None\n")
}

cat("\nES types with SOME files having out-of-range values:\n")
some_affected <- es_analysis[
  es_analysis$Percentage > 0 & es_analysis$Percentage < 100,
]
if (nrow(some_affected) > 0) {
  for (i in 1:nrow(some_affected)) {
    cat(
      "  -",
      some_affected$ES_Type[i],
      ":",
      some_affected$Percentage[i],
      "% (",
      some_affected$Files_With_OutOfRange[i],
      "of",
      some_affected$Total_Files[i],
      "files )\n"
    )
  }
} else {
  cat("  None\n")
}

# Show sample files for ES types with highest percentage
cat("\n=== SAMPLE FILES (Top 3 ES types by percentage) ===\n")
top_es <- head(es_analysis[es_analysis$Percentage > 0, ], 3)

for (i in 1:nrow(top_es)) {
  es_type <- tolower(top_es$ES_Type[i])
  cat(
    "\n",
    top_es$ES_Type[i],
    " (",
    top_es$Percentage[i],
    "% affected):\n",
    sep = ""
  )

  es_out_of_range_files <- grep(
    paste0("^", es_type, "-"),
    files_with_out_of_range,
    value = TRUE
  )
  sample_files <- head(es_out_of_range_files, 5)

  for (file in sample_files) {
    cat("  - ", file, "\n")
  }

  if (length(es_out_of_range_files) > 5) {
    cat("  ... and", length(es_out_of_range_files) - 5, "more files\n")
  }
}

# Additional analysis: Check if certain scenarios/configurations are more prone to out-of-range
cat("\n=== SCENARIO ANALYSIS ===\n")
cat(
  "Checking if certain scenario components are more prone to out-of-range values...\n\n"
)

# Extract scenario components from filenames
# Typical format: es-year-rcp-scenario-intensity-intervention-canton-perc_area.json
scenario_components <- list(
  rcp = c("rcp26", "rcp45", "rcp85"),
  scenario = c(
    "ref_central",
    "ref_peri_urban",
    "ecolo_central",
    "ecolo_urban",
    "combined_urban"
  ),
  intensity = c("ref", "low", "high"),
  intervention = c("bau", "ei_cul", "ei_nat", "ei_soc", "gr_ex")
)

for (component_type in names(scenario_components)) {
  cat(toupper(component_type), "analysis:\n")

  for (component in scenario_components[[component_type]]) {
    # Count files with this component
    component_files <- grep(component, all_successful_files, value = TRUE)
    component_out_of_range <- grep(
      component,
      files_with_out_of_range,
      value = TRUE
    )

    if (length(component_files) > 0) {
      percentage <- round(
        length(component_out_of_range) / length(component_files) * 100,
        1
      )
      cat(
        "  ",
        component,
        ": ",
        percentage,
        "% (",
        length(component_out_of_range),
        " of ",
        length(component_files),
        " files)\n"
      )
    }
  }
  cat("\n")
}

cat("Analysis completed.\n")
