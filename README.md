# VSALM-paper

This repository contains code for the analyses described in the manuscript "A spatial variance-smoothing area level model for small area estimation of demographic rates" by Peter A. Gao and Jon Wakefield [[arXiv to come]].

## Data

The following data sources are available freely or upon request online:

* Data from the [2018 Nigeria DHS survey](https://dhsprogram.com/publications/publication-fr359-dhs-final-reports.cfm) can be requested 
[here](https://dhsprogram.com/data/dataset/Nigeria_Standard-DHS_2018.cfm). Data from the [2015-16 Malawi DHS survey](https://dhsprogram.com/publications/publication-FR319-DHS-Final-Reports.cfm) can be requested 
[here](https://dhsprogram.com/data/dataset/Malawi_Standard-DHS_2015.cfm).

* Boundary data for Nigeria and Malawi's administrative divisons (`data/Nigeria/shapeFiles_gadm/` and `data/Malawi/shapeFiles_gadm/`) was obtained from the [Database of Global Administrative Areas (GADM)](https://gadm.org/).

* Age-specific population density estimates for Nigeria were obtained from [WorldPop](https://hub.worldpop.org/geodata/summary?id=16385).

* `nga_2018_ea_frame_strata_sample_size.csv` is a file containing information on number of sampled enumeration areas in each stratum with columns (see [Table A.3](https://dhsprogram.com/pubs/pdf/FR359/FR359.pdf)):
  - `DHS`: DHS  region name
  - `GADM`: GADM region name
  - `Urban`: Number of urban EAs sampled
  - `Rural`: Number of rural EAs sampled

The data should be organized in the following tree structure:

```bash
├── data
│   ├── Malawi
│   │   ├── MW_DHS
│   │   │   ├── MWAR7ADT
│   │   │   │   ├── MWAR7AFL.DCT
│   │   │   │   ├── MWAR7AFL.DO
│   │   │   │   ├── MWAR7AFL.DTA
│   │   │   │   └── MWAR7AFL.MAP
│   │   │   ├── MWGE7AFL
│   │   │   │   ├── DHS_README.txt
│   │   │   │   ├── GPS_Displacement_README.txt
│   │   │   │   ├── MWGE7AFL.CPG
│   │   │   │   ├── MWGE7AFL.dbf
│   │   │   │   ├── MWGE7AFL.prj
│   │   │   │   ├── MWGE7AFL.sbn
│   │   │   │   ├── MWGE7AFL.sbx
│   │   │   │   ├── MWGE7AFL.shp
│   │   │   │   ├── MWGE7AFL.shp.xml
│   │   │   │   └── MWGE7AFL.shx
│   │   │   └──  MWIR7HDT
│   │   │       ├── MWIR7HFL.DCT
│   │   │       ├── MWIR7HFL.DO
│   │   │       ├── MWIR7HFL.DTA
│   │   │       ├── MWIR7HFL.FRQ
│   │   │       ├── MWIR7HFL.FRW
│   │   │       └── MWIR7HFL.MAP
│   │   └── shapeFiles_gadm
│   │       ├── gadm36_MWI_0.dbf
│   │       ├── gadm36_MWI_0.prj
│   │       ├── gadm36_MWI_0.shp
│   │       ├── gadm36_MWI_0.shx
│   │       ├── gadm36_MWI_1.dbf
│   │       ├── gadm36_MWI_1.prj
│   │       ├── gadm36_MWI_1.shp
│   │       ├── gadm36_MWI_1.shx
│   │       ├── gadm36_MWI_2.dbf
│   │       ├── gadm36_MWI_2.prj
│   │       ├── gadm36_MWI_2.shp
│   │       ├── gadm36_MWI_2.shx
│   │       ├── gadm36_MWI_3.cpg
│   │       ├── gadm36_MWI_3.dbf
│   │       ├── gadm36_MWI_3.prj
│   │       ├── gadm36_MWI_3.shp
│   │       ├── gadm36_MWI_3.shx
│   │       └── license.txt
│   └── Nigeria
│       ├── NGGE7BFL
│       │   ├── NGGE7BFL.cpg
│       │   ├── NGGE7BFL.dbf
│       │   ├── NGGE7BFL.prj
│       │   ├── NGGE7BFL.sbn
│       │   ├── NGGE7BFL.sbx
│       │   ├── NGGE7BFL.shp
│       │   └── NGGE7BFL.shx
│       ├── NGKR7BFL.DTA
│       ├── 2015_accessibility_to_cities_v1.0.tif
│       ├── nga10povcons200.tif
│       ├── nga_2018_poppa.csv
│       ├── nga_f_0_2006.tif
│       ├── nga_f_1_2006.tif
│       ├── nga_f_1_2018.tif
│       ├── nga_m_0_2006.tif
│       ├── nga_m_1_2006.tif
│       ├── nga_m_1_2018.tif
│       ├── nga_ppp_2006_UNadj.tif
│       ├── nga_ppp_2018_UNadj.tif
│       └── shapeFiles_gadm
│           ├── gadm36_NGA_0.cpg
│           ├── gadm36_NGA_0.dbf
│           ├── gadm36_NGA_0.prj
│           ├── gadm36_NGA_0.shp
│           ├── gadm36_NGA_0.shx
│           ├── gadm36_NGA_1.dbf
│           ├── gadm36_NGA_1.prj
│           ├── gadm36_NGA_1.shp
│           ├── gadm36_NGA_1.shx
│           ├── gadm36_NGA_2.dbf
│           ├── gadm36_NGA_2.prj
│           ├── gadm36_NGA_2.shp
│           ├── gadm36_NGA_2.shx
│           └── license.txt
```

## Analysis

* The associated [`VSALM` package]

* `analysis/manuscript_figures.R` contains code for generating some figures used in the manuscript.

* `analysis/Nigeria` contains scripts for generating Admin-1 level estimates of 
measles vaccination rates for Nigeria based on the 2018 DHS survey.

  - `analysis/Nigeria/process-population-covariates.R` generates R objects containing population information for all areas.
  - `analysis/Nigeria/Nigeria-joint-model-MCV.R` runs analysis and generates figures for the Nigeria example.
  
* `analysis/Malawi` contains scripts for generating Admin-1 level estimates of 
HIV prevalence rates for Malawi based on the 2015-16 DHS survey.

* `analysis/simulations` contains code for reproducing the simulations described 
in the manuscript.

## Results

Outputs are collected in the `results` and `paper/figures` folders.

## Notes

The above scripts were run on a cluster computer using R version 4.1.2.

```
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8   
 [6] LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C        
[11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.1.1    readstata13_0.10.0 SUMMER_1.2.0       raster_3.5-11      rgdal_1.5-28       purrr_0.3.4       
 [7] survey_4.1-1       survival_3.2-13    INLA_21.11.22      foreach_1.5.1      Matrix_1.3-4       spdep_1.2-2       
[13] spData_2.0.1       sp_1.4-6           sf_1.0-5           ggplot2_3.3.5      dplyr_1.0.7        tidyr_1.1.4       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8         lattice_0.20-45    deldir_1.0-6       class_7.3-19       assertthat_0.2.1   utf8_1.2.2        
 [7] R6_2.5.1           e1071_1.7-9        pillar_1.6.5       rlang_1.0.0        data.table_1.14.2  splines_4.1.2     
[13] munsell_0.5.0      proxy_0.4-26       compiler_4.1.2     xfun_0.29          pkgconfig_2.0.3    mitools_2.4       
[19] tidyselect_1.1.1   tibble_3.1.6       gridExtra_2.3      codetools_0.2-18   fansi_1.0.2        viridisLite_0.4.0 
[25] crayon_1.4.2       withr_2.4.3        wk_0.5.0           gtable_0.3.0       lifecycle_1.0.1    DBI_1.1.2         
[31] magrittr_2.0.2     units_0.7-2        scales_1.1.1       KernSmooth_2.23-20 cli_3.1.1          viridis_0.6.2     
[37] ellipsis_0.3.2     generics_0.1.1     vctrs_0.3.8        boot_1.3-28        RColorBrewer_1.1-2 s2_1.0.7          
[43] iterators_1.0.13   tools_4.1.2        glue_1.6.1         shadowtext_0.1.1   colorspace_2.0-2   terra_1.4-22      
[49] classInt_0.4-3     knitr_1.37 

```
