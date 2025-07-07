# List of your required packages
required_packages <- c(
  "lme4",        # For fitting linear mixed-effects models
  "pbkrtest",    # For Kenward-Roger approximation for mixed models
  "simpleboot",  # For simple bootstrap methods
  "lmeresampler",# For resampling methods for mixed models
  "foreach",     # For looping constructs for parallel execution
  "doParallel",  # For parallel backend to the foreach
  "doRNG",       # For parallel Reproducibility
  "dplyr",         # Data manipulation
  "pbkrtest",      # Parametric bootstrap and Kenward-Roger approximation
  "lmeresampler",  # Resampling methods for mixed models
  "emmeans",       # Estimated marginal means (least-squares means)
  "parameters",    # Model parameter extraction and summary
  "bayestestR",    # Bayesian model diagnostics and summaries
  "ggplot2",       # Data visualization
  "ggbeeswarm",    # Quasirandom plots (e.g., beeswarm-style jittering)
  "knitr",         # Report generation
  "readxl",   # For reading Excel files
  "janitor",   # For cleaning data
  "car",       # For regression diagnostics
  "sjPlot",    # For data visualization and tabulation
  "tictoc",    # For timing code execution
  "skimr",     # For summarizing data
  "here",
  "gridExtra",
  "tidyverse",
  "tinytex",
  "extrafont",
  "ggpubr",
  "ggsci",
  "freshr"
)

# Initialize an empty list to store citations
all_citations <- list()

# Loop through each package and get its citation
for (pkg in required_packages) {
  # Try to get the citation for the package
  # Some packages might not have a citation or might throw an error if not installed
  tryCatch({
    pkg_citation <- citation(pkg)
    all_citations[[pkg]] <- pkg_citation
    message(paste("Successfully retrieved citation for:", pkg))
  }, error = function(e) {
    warning(paste("Could not retrieve citation for", pkg, ":", e$message))
    message(paste("Consider installing", pkg, "if not already, or manually citing it."))
  })
}

# Define the output BibTeX file name
bib_file_name <- "R_package_citations.bib"

# Open the file in write mode
file_conn <- file(bib_file_name, "w")

# Write each citation to the BibTeX file
for (pkg_citation in all_citations) {
  # Convert the citation object to BibTeX format and write it to the file
  # The 'toBibtex' method is usually available for citation objects
  writeLines(toBibtex(pkg_citation), file_conn)
  writeLines("\n", file_conn) # Add a newline for separation
}

# Close the file connection
close(file_conn)

message(paste("All available package citations saved to:", bib_file_name))

# --- Important Note for R Itself ---
# Don't forget to cite R itself!
message("\nDon't forget to cite R itself:")
print(citation())
writeLines(toBibtex(citation()), file("R_citations_base.bib", "w"))
message("Base R citation saved to R_citations_base.bib")


cat("--- R Package Versions ---\n")
for (pkg in required_packages) {
  # Check if the package is installed first to avoid errors for missing packages
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("%s: %s\n", pkg, as.character(version)))
  } else {
    cat(sprintf("%s: Not installed\n", pkg))
  }
}
cat("--------------------------\n")
