# 🍳 simChef culinary school

## Cooking up reliable simulations in R

- [Slides](https://tiffanymtang.github.io/simChef-culinary-school)
- [Repo](https://github.com/Yu-Group/simChef)
- [Docs](https://yu-group.github.io/simChef/)
- [Example simulation](https://github.com/PhilBoileau/simChef-case-study)

## Authors

- Tiffany Tang, UC Berkeley Statistics -- developer, maintainer
- James Duncan, UC Berkeley Biostatistics -- developer
- Corrine F. Elliott, UC Berkeley Statistics -- contributer
- Philippe Boileau, UC Berkeley Biostatistics -- contributer
- Bin Yu, UC Berkeley Statistics, EECS, and Biostatistics -- thesis adviser

## Getting started with the case study

### 1. Clone this repository

In a terminal, run:

```bash
git clone https://github.com/tiffanymtang/simChef-culinary-school.git
cd simChef-culinary-school
```

Open the cloned `simChef-culinary-school` folder in RStudio or your preferred R
development environment. Make sure the repository root is your working
directory.

### 2. Install the dependencies

Start R from the repository root, then restore the project environment:

```r
renv::restore()
```

This installs the package versions recorded in `renv.lock`.

### 3. Run the interaction-discovery case study

From the repository root, run:

```r
source("case-studies/interaction_discovery.R")
```

The script runs 20 simulation replicates for each sample size and can take a
while to finish. It saves the results under `case-studies/results/`.

### Download the completed results instead

To skip fitting the simulation, download the cached results from
[Zenodo](https://zenodo.org/records/23861655). Run the following code from the
repository root. It places the files in the same directory structure created by
the case-study script.

```r
record_id <- "23861655"
result_files <- c(
  "experiment.rds",
  "experiment_cached_params.rds",
  "fit_results.rds"
)

results_dir <- file.path(
  "case-studies", "results", "Example iRF Simulation",
  "Enhancer LSS with RF Noise", "Varying n"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Allow enough time to download the largest file.
options(timeout = max(600, getOption("timeout")))

for (filename in result_files) {
  download.file(
    url = sprintf(
      "https://zenodo.org/records/%s/files/%s?download=1",
      record_id, filename
    ),
    destfile = file.path(results_dir, filename),
    mode = "wb"
  )
}
```

Then run the following code, which will used these cached results to generate the simulation figures and documentation:

```r
source("case-studies/interaction_discovery.R")
```
