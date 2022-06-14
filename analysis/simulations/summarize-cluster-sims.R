#### 0 LIBRARIES AND PATHS #####################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)

#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}

### SET TO LOCATION
sim_res_dir <- "../BALM-SAE/results/cluster/sims/area-level-models/"
summarize_res <- function(res) {
  res %>% 
    group_by(method) %>%
    summarize(rmse = sqrt(mean((median - pop_mean)^2)),
              abs_bias = mean(abs(median - pop_mean)),
              cov90 = mean(lower < pop_mean & upper > pop_mean),
              int_len = mean(upper - lower))
}

pub_order <- c(
  "Hajek",
  "iidMeanSmooth",
  "iidMeanSmoothLogit",
  "iidMeanSmoothUnmatched",
  "spatialMeanSmooth",
  "spatialMeanSmoothLogit",
  "spatialMeanSmoothUnmatched",
  "iidJointSmooth",
  "iidJointSmoothLogit",
  "iidJointSmoothUnmatched",
  "spatialJointSmooth",
  "spatialJointSmoothLogit",
  "spatialJointSmoothUnmatched"
)

selected_methods <-
  c(
    "Hajek",
    "iidMeanSmooth",
    "iidMeanSmoothUnmatched",
    "spatialMeanSmooth",
    "spatialMeanSmoothUnmatched",
    "iidJointSmooth",
    "iidJointSmoothUnmatched",
    "spatialJointSmooth",
    "spatialJointSmoothUnmatched"
  )

fmt_tbl <- function(res, methods = unique(res$method), 
                    rows = NULL) {
  res <- res %>%
    filter(method %in% methods) %>%
    arrange(match(method, pub_order)) %>%
    mutate(rmse = rmse * 100, abs_bias = abs_bias * 100, 
           int_len = int_len * 100, cov90 = cov90 * 100) %>%
    setNames(c("Method", 
               paste0("RMSE (x 100)"),
               paste0("MAE (x 100)"),
               "90% Cov.",
               "Int. Len. (x 100)")) %>%
    knitr::kable(digits = c(0, 2, 2, 0, 2), format = "latex", booktabs = T,
                 linesep = "")
  if (!is.null(rows)) {
    # res <- res %>%  kableExtra::row_spec(rows, hline_after = T) 
  }
  return(res)
  
}



#### 3 SPATIAL MODEL ###########################################################
spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "area-level-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "area-level-model_res_", x, ".rds"))
              }
            })) 
fmt_tbl(summarize_res(spa_res)) %>% writeLines("paper/figures/area-level-model-res.tex")
loprev_spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "loprev-area-level-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "loprev-area-level-model_res_", x, ".rds"))
              }
            })) 
fmt_tbl(summarize_res(loprev_spa_res)) %>% writeLines("paper/figures/loprev-area-level-model-res.tex")