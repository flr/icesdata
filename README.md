# icesdata: ICES Age-Structured Stock Assessment Data and Utilities

A comprehensive R package for ICES age-structured stock assessment data, including utility functions for data processing, analysis, and visualization. This package consolidates functionality from multiple ICES-related packages into a single, maintainable codebase.

## Features

- **Data Loading**: Functions to load and process ICES stock assessment data from various formats
- **Analysis Tools**: Utilities for calculating benchmarks, time series, and management plan statistics
- **Visualization**: Plotting functions for stock status, fishing mortality, and reference points
- **Core Utilities**: FLR-compatible functions for stock assessment workflows

## Installation

```r
# Install from GitHub (if available)
devtools::install_github("your-username/icesdata")

# Or install locally
devtools::install("path/to/icesdata")
```

## Usage

### Loading Data

```r
library(icesdata)

# Load ICES stock data
stock_data <- loadIcesStockData("path/to/stock_data.csv")

# Load ICES RData files
loadIcesRdata("data/icesdata.RData")

# Process ICES data
processed_data <- processIcesData("data/data.RData", "data/info.RData")
```

### Analysis

```r
# Calculate benchmark reference points
benchmarks <- calculateBenchmarks(processed_data)

# Calculate time series
ts_data <- calculateTimeseries(processed_data, benchmarks)

# Calculate management plan statistics
mp_stats <- calculateMpStats(ts_data)
```

### Visualization

```r
# Plot stock status over time
plotStockStatus(stock_data)

# Plot fishing mortality
plotFishingMortality(stock_data)

# Plot management plan statistics
plotManagementPlanStats(mp_stats)

# Plot benchmark reference points
plotBenchmarks(benchmarks)
```

## Package Structure

```
icesdata/
├── R/                          # R source code
│   ├── data-loading.R         # Data loading functions
│   ├── analysis-functions.R   # Analysis utilities
│   ├── plotting-functions.R   # Visualization functions
│   └── icesdata-utils.R      # Core utility functions
├── data/                      # Package datasets
├── vignettes/                 # Package vignettes
├── tests/                     # Unit tests
└── man/                       # Documentation
```

## Dependencies

### Required
- R (>= 4.0.0)
- methods
- stats
- utils

### Suggested
- testthat (>= 3.0.0)
- knitr
- rmarkdown
- ggplot2
- dplyr
- tidyr
- plyr

## Dataset Storage

**Note**: The datasets in this package are too large for GitHub. They are stored externally and can be accessed through:

1. **Zenodo**: [DOI: TBD] - Long-term storage for research data
2. **ICES Data Portal**: Official ICES data repository
3. **Local Network**: For internal institutional use
4. **Cloud Storage**: AWS S3, Google Cloud, or Azure Blob Storage

To use the package with datasets, download them from the external source and place them in the `data/` directory.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## License

This package is licensed under the MIT License. See the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
ICES Working Group. (2024). icesdata: ICES Age-Structured Stock Assessment Data and Utilities. 
R package version 0.1.0.
```

## Contact

For questions or issues, please contact:
- **Maintainer**: ICES Working Group <ices@ices.dk>
- **Issues**: [GitHub Issues](https://github.com/your-username/icesdata/issues)

## Acknowledgments

This package consolidates functionality from:
- `bimResilience` package
- `BRebuildUtils` package
- Original `icesdata` utilities

Special thanks to the ICES community for providing the underlying data and methodologies.
