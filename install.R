# install.R

# Run this script whithin the cecast main directory.
# R < install.R --no-save  

# Step 1: Install required packages for cecast
required_packages <- c("argparse", "parallel", "RColorBrewer", "fields", "here", "jsonlite")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load the necessary packages for here
library("here", quietly=T)
library("jsonlite", quietly=T)


# Step 2: Check and create necessary directories
dirs <- c("data", "scripts")
for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    cat("Created directory:", dir, "\n")
  } else {
    cat("Directory already exists:", dir, "\n")
  }
}


# Step 3: Download external resources (if necessary)
main_url <- "https://github.com/cesaredef/cecast/"
counts <- c(
  "expected_counts_boot.RData", # previously cecast_exp_boot.RData
  "expected_counts_fine_grid.RData", 
  "expected_counts_fine_grid_archaics_afr_combined.RData",
  "expected_counts_fine_grid_archaics_var.RData",
  "some_human_counts.RData" # previously some_human_frequencies.RData
  )

local_files <- c(here("data", counts), here("scripts", "cecast_functions.R"))
url_files <- paste(main_url, c(rep("data/", length(counts)), "scripts/"), c(counts,"cecast_functions.R"),  sep="")

for (f in seq_along(local_files)) {
	if (!file.exists(local_files[f])) {
		download.file(url_files[f], local_files[f])
		cat("Downloaded file to:", local_files[f], "\n")
	} else {
		cat("File already exists:", local_files[f], "\n")
	}
}


# Step 4: Create or check configuration file 
config_file <- here("~/.config_cecast.json")
work_dir <- getwd()
# If ~/.config_cecast.json doesn't exist, create it with default values
if (!file.exists(config_file)) {
  default_config <- list(
    data_dir = paste(work_dir, "/data", sep=""),
    data_files = counts[c(1,2,5)],
	  scripts_dir = paste(work_dir, "/scripts", sep=""),
    functions_file = "cecast_functions.R"
  )
  write_json(default_config, config_file)
  cat("Created default config file:", config_file, "\n")
} else {
  cat("Config file already exists:", config_file, "\n")
}

