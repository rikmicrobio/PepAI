# setup_environment.R
# Run this once to install all PepAI dependencies
# Usage: source("setup_environment.R")

cat("==============================================\n")
cat("  PepAI — Environment Setup\n")
cat("  InteractIQ — AI-first Drug Discovery\n")
cat("==============================================\n\n")

required_pkgs <- c(
  "shiny",
  "shinydashboard",
  "visNetwork",
  "ggplot2",
  "DT",
  "dplyr",
  "plotly",
  "stringr",
  "scales",
  "htmlwidgets"
)

cat("Checking required packages...\n")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]

if (length(missing_pkgs) == 0) {
  cat("All packages already installed.\n\n")
} else {
  cat(sprintf("Installing %d missing package(s): %s\n\n",
              length(missing_pkgs),
              paste(missing_pkgs, collapse = ", ")))
  install.packages(missing_pkgs,
                   repos    = "https://cloud.r-project.org/",
                   dependencies = TRUE)
}

# Verify all packages load
cat("Verifying package loading...\n")
errors <- character(0)
for (pkg in required_pkgs) {
  tryCatch(
    library(pkg, character.only = TRUE, quietly = TRUE),
    error = function(e) errors <<- c(errors, pkg)
  )
}

if (length(errors) == 0) {
  cat("All packages loaded successfully.\n\n")
  cat("==============================================\n")
  cat("  Setup complete! Run the app with:\n")
  cat("  shiny::runApp('app.R')\n")
  cat("==============================================\n")
} else {
  cat(sprintf("WARNING: Failed to load: %s\n", paste(errors, collapse = ", ")))
  cat("Try installing manually:\n")
  cat(sprintf("  install.packages(c('%s'))\n", paste(errors, collapse = "', '")))
}
