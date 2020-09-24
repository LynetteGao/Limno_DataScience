# odem.data

This R-package is still under development and functions are likely to change with future updates. `odem.data` is maintained by UW-Madison Center for Limnology. All data were obtained from online repositories, e.g. NTL-LTER and the USGS Water Quality Portal.

## Installation

```r
remotes::install_github("LynetteGao/odem.data")
```

## Usage

```R
# build odem input data for a specific lake using R
ii <- "analysis/data/143249470_Mendota/"
source("analysis/scripts/01_data_merge.R")
```

```bash
# or build odem input data for a specific lake using the command line
make analysis/data/143249470_(Mendota)/input.txt

# build odem input data for all lakes using the command line
make all
```

## Data

### Input data

 * temperature: `pball_*_temperatures.csv`
  
 * meteorology: `NLDAS_*.csv`
 
 * water quality data from LTER: `wq_data_*.csv`
 
 * nml file : `nhdhr_*.nml`

### Output data for odem

|variable                 |
|:-----------------|
|datetime          |
|thermocline_depth |
|temperature_epi   |
|temperature_hypo  |
|temperature_total |
|volume_total      |
|volume_epi        |
|volume_hypo       |
|area_thermocline  |
|area_surface      |
|upper_meta        |
|lower_meta        |
|year              |
|day_of_year       |
|max.d             |
|wind              |
|airtemp           |

## Postscript

This project is part of the Data Science Initiative at UW-Madison.
