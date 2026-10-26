# Rainfall Shocks, Child Anthropometrics, and Remittances (MS Thesis)

**Zachary Serby, Cornell University**

> **This project is fully reproducible. All data needed to run the analysis are included in this repository, and every file path is relative. Thus, the project can be downloaded to any computer and run as is, with no paths or settings to change. Running the single master file, `Scripts/00_Master.do`, rebuilds the cleaned survey data and every regression, table, figure, and the results PDF from the raw data.**

**Paper:** [Writing Sample.pdf](Writing%20Sample.pdf) (current draft of the results, also at [Paper/Results.pdf](Paper/Results.pdf))
**Code:** [Scripts/](Scripts) (start with [Scripts/00_Master.do](Scripts/00_Master.do))

---

## Overview

This thesis studies how the 2009 drought in Mexico affected child health, measured by WHO height-for-age (HAZ), weight-for-age (WAZ), and BMI-for-age (BMIZ) z-scores, and whether household remittances buffered the shock. It combines the Mexican Family Life Survey (2002, 2005, 2009) with municipality-level CHIRPS rainfall and estimates binary difference-in-differences, triple differences with 2002 remittances, event studies, and the continuous difference-in-differences of Callaway, Goodman-Bacon and Sant'Anna (2024), with subanalyses by child sex and by rural vs urban locality.

## How to run

1. Download or clone this repository.
2. Open `Scripts/00_Master.do` (double-clicking it sets Stata's working folder to `Scripts`; if Stata is already open, `cd` to the project or `Scripts` folder first). The project root is found automatically.
3. Run it. The full pipeline takes about 7 minutes and runs, in order: MxFLS cleaning (2002, 2005, 2009), the DID analysis data, summary statistics, DID/DDD and event studies, continuous DID (R), subanalyses by sex and by rural vs urban (Stata, then R), the regression tables, and the results document (R, compiled to `Paper/Results.pdf`).

Each step can be switched on or off with the `run_*` switches at the top of `00_Master.do` (all are on). The master calls R through `Rscript` and finds it automatically in the usual Mac locations; on another system, set `global Rscript` in `00_Master.do`.

### Rainfall data

The cleaned municipality rainfall panel (`Data/Cleaned/Rainfall/municipality_year_rainfall.dta`) is included, so the analysis runs without the raw rainfall file. The raw CHIRPS file (7.7 GB) is not included because it exceeds GitHub's file size limit. To rebuild the rainfall panel from scratch, download it from the UCSB Climate Hazards Center:

https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/netcdf/chirps-v2.0.monthly.nc

place it at `Data/Raw/Rainfall/chirps-v2.0.monthly.nc`, and set `run_R` to 1 in `00_Master.do`. If the cleaned rainfall file is ever missing, the master rebuilds it automatically from the raw file.

## Repository structure

```
Scripts/
  00_Master.do                      Runs everything in order
  MS Thesis.Rproj                   R project (R scripts take the project root as one folder up from Scripts)
  Cleaning/
    Cleaned_rainfall_municipality_drought.R   CHIRPS to a municipality-year rainfall panel and drought indicator
    Cleaning Script for MxFLS_2002_03.do      Cleans each MxFLS wave: WHO z-scores, remittances, controls
    Cleaning Script for MxFLS_2005_06.do
    Cleaning Script for MxFLS_2009_12.do      Also appends the three waves into the panel
    Build_DID_Analysis_Data.do                Merges the panel with rainfall, builds treatment, controls, samples
  Statistics/
    Summary Stats DID.do            Summary statistics, sample sizes, attrition tables
    Regression Stats.do             All regression tables in LaTeX
  Regressions/
    DID_Regs.do                     Binary DID, triple difference, event studies, robustness
    Contdid_regs.R                  Continuous DID, SPI drought classes, remittance buffer
    functions.R                     Functions used by the R scripts
    Subanalysis_Gender_*            Same analyses for boys and girls
    Subanalysis_Rural_*             Same analyses for rural and urban children
  Paper/
    Build_Results.R                 Writes Paper/Results.tex and compiles Paper/Results.pdf

Data/
  Raw/        MxFLS data and documentation (codebooks, questionnaires, user guides),
              INEGI municipal shapefiles, WHO z-score macros (never edited by any script)
  Cleaned/    Produced by the scripts: MxFLS wave files and panel, DID analysis data,
              WHO z-score files, municipality rainfall panel

Output/
  Tables/                 Summary statistics and regression tables (.tex), Main and Appendix
  Regression_Estimates/   Event study coefficients (.dta) and continuous DID estimates (.csv)
  Figures/                Trends, event studies, continuous DID figures (.png)
  Logs/                   Every Stata and R log

Paper/
  Results.tex, Results.pdf          All main and appendix tables and the main figures
  Proposal.pdf                      Thesis proposal

README.rtf                          Detailed project readme
Writing Sample.pdf                  Copy of Paper/Results.pdf
```

## Data sources

- **Mexican Family Life Survey (MxFLS)**: MxFLS-1 (2002), MxFLS-2 (2005-06), MxFLS-3 (2009-12). Data, codebooks, questionnaires and user guides are publicly available at no charge from [ennvih-mxfls.org](http://www.ennvih-mxfls.org).
- **CHIRPS v2.0** monthly precipitation, 0.05 degree grid, UCSB Climate Hazards Center: [chirps-v2.0.monthly.nc](https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/netcdf/chirps-v2.0.monthly.nc).
- **INEGI** Marco Geoestadistico 2010 municipal boundaries (Municipios_2010_5).
- **WHO** Child Growth Standards (igrowup) and 2007 Growth Reference (who2007) Stata macros.

## Software requirements

- **Stata** (tested with Stata 19 MP). reghdfe, ftools, require, estout and distinct install automatically.
- **R** 4.x (tested with R 4.4.2). terra, sf, data.table, readxl, haven, did (>= 2.5.0), BMisc, sandwich and tinytex install automatically. contdid 0.1.1 installs from GitHub (bcallaway11/contdid) and needs ptetools >= 1.0.1. fwildclusterboot 0.13.0 and summclust 0.7.2 install from the CRAN archive and need the Xcode command line tools and gfortran on a Mac.
- **TinyTeX** compiles Results.pdf; `Build_Results.R` installs it if no LaTeX is found.

Bootstrap results are reproducible: `set.seed(1234)` is set before every bootstrap in R.
