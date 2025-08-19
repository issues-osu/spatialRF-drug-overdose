# Spatial Random Forest Modeling of Drug-Related Mortality Rates Using Environmental and Socioeconomic Predictors

![Study Map](map.png)

This repository contains R code for modeling county-level **drug-related mortality rates** in the Atlantic states using **spatial random forest** and **geographically weighted random forest (GWRF)** approaches. The workflow integrates CDC mortality data, environmental justice indicators, and socioeconomic predictors to evaluate spatial heterogeneity in risk factors.

---

## Project Structure

### 1. Libraries and Setup
The analysis relies on spatial, statistical, and machine learning packages including:

- Spatial data management: `sf`, `tigris`, `sp`, `raster`
- Spatial statistics: `spdep`, `GWmodel`, `SpatialML`
- Machine learning: `h2o`, `vip`, `pdp`
- Visualization: `ggplot2`, `RColorBrewer`, `gridExtra`, `fmsb`

Working directory and cache options are set to ensure reproducibility.

---

### 2. Data Acquisition
- **Geographic data:** Downloaded using the `tigris` package for Atlantic states (NY, NJ, NC, SC, GA, DE, MD, PA, VA, WV).  
- **Mortality data:** County-level drug-related crude rates from CDC CSV extract.  
- **Socioeconomic/environmental data:** Indicators merged from an EJI dataset.  
- County geometries and data are merged by **FIPS code**.

---

### 3. Data Preparation
- Coordinate system transformed to **EPSG:5070 (Albers Conus)**.  
- County centroids (X, Y) computed for use in spatial models.  
- Data reshaped into wide and long formats for correlation analysis and modeling.  
- Predictors standardized across training, validation, and test splits (2018, 2020, 2022).

---

### 4. Random Forest Modeling (H2O)
- **Training:** H2O distributed random forest (`h2o.randomForest`) with hyperparameter tuning (`h2o.grid`).  
- **Validation/Test:** RMSE, MAE, and R² statistics computed.  
- **Visualization:** Observed vs. predicted plots for validation (2020–2021) and test (2018–2019).  

Outputs include:
- Cross-validation metrics  
- Feature importance (`vip`)  
- Partial dependence plots (`pdp`)  

---

### 5. Geographically Weighted Random Forest (GWRF)
- Local model fit with **adaptive kernel weighting**.  
- Outputs:
  - Local R² maps  
  - Local feature importance (IncMSE) for unemployment, uninsurance, disability, parks, mental health, mobility, pollution, TRI, and walkability.  
- Validation and test performance metrics (RMSE, MAE, R²) compared against global RF.

---

### 6. Visualization
- **Choropleths:** Predicted and observed mortality rates by county for training, validation, and test sets.  
- **Feature importance plots:** Permutation-based variable ranking.  
- **Partial dependence plots:** Top predictors individually and in interaction.  
- **Radar charts:** State-level summaries of socioeconomic and environmental indicators.

---

## Outputs
- **Maps:** Predicted vs observed drug-related mortality rates (`spplot` outputs).  
- **Variable importance:** Global and local feature rankings.  
- **Diagnostics:** RMSE, MAE, R² for both RF and GWRF.  
- **Figures:** Radar charts of state-level predictor profiles.

---

## 🔧 Requirements
- R version ≥ 4.0  
- H2O cluster initialized with sufficient memory (tested with `max_mem_size = "48g"`).  
- Required R packages:  
  `sf`, `tigris`, `spdep`, `SpatialML`, `h2o`, `ggplot2`, `gridExtra`, `vip`, `pdp`, `fmsb`.

---

## How to Run
1. Clone this repository.  
2. Place input files in the working directory:  
   - `cdc-drug-deaths-atlantic-all.csv`  
   - `eji-indicators.csv`  
   - `STATE_ATLANTIC.shp`  
3. Run the R script in sequence:  
   - Data download and preparation  
   - Random forest modeling  
   - GWRF modeling  
   - Visualization and summary statistics  

---

## Citation
This work is based on the tutorial provided here:
[**"Spatial Random Forest Modeling of Drug-Related Mortality Rates Using Environmental and Socioeconomic Predictors"**](https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html)

---
