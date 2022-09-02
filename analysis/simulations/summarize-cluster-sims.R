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
#setwd(paste0(home_dir, "VSALM-paper"))
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
pub_names <- rbind(
  c("Hajek", "Direct (Hájek)"),
  c("iidMeanSmooth", "MS"),
  c("iidMeanSmoothLogit", "Logit MS"),
  c("iidMeanSmoothUnmatched", "Unmatched MS"),
  c("spatialMeanSmooth", "Spatial MS"),
  c("spatialMeanSmoothLogit", "Spatial Logit MS"),
  c("spatialMeanSmoothUnmatched", "Spatial Unmatched MS"),
  c("iidJointSmooth", "JS"),
  c("iidJointSmoothLogit", "Logit JS"),
  c("iidJointSmoothUnmatched", "Unmatched JS"),
  c("spatialJointSmooth", "Spatial JS"),
  c("spatialJointSmoothLogit", "Spatial Logit JS"),
  c("spatialJointSmoothUnmatched", "Spatial Unmatched JS")
) %>%
  as.data.frame() %>%
  setNames(c("internal", "publication")) 
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
    mutate(method = 
             pub_names$publication[match(method, pub_names$internal)]) %>%
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
fmt_tbl(summarize_res(spa_res), methods = selected_methods) %>%
  writeLines("paper/figures/area-level-model-res.tex")
loprev_spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "loprev-area-level-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "loprev-area-level-model_res_", x, ".rds"))
              }
            })) 
fmt_tbl(summarize_res(loprev_spa_res), methods = selected_methods) %>%
  writeLines("paper/figures/loprev-area-level-model-res.tex")

sim_res_dir <- "results/cluster/sims/lgsample-area-level-models/"
spa_res <- 
  do.call(rbind,
          lapply(
            1:1000,
            function(x) {
              if (file.exists(paste0(sim_res_dir, "lgsample-area-level-model_res_", x, ".rds"))) {
                readRDS(paste0(sim_res_dir, "lgsample-area-level-model_res_", x, ".rds"))
              }
            })) 
fmt_tbl(summarize_res(spa_res)) %>% writeLines("results/figures/lgsample-area-level-model-res.tex")


summarize_res(loprev_spa_res) %>%
  left_join(summarize_res(spa_res) %>% 
              setNames(paste0("med_", colnames(.))), 
            by =  c("method" = "med_method")) %>%
  arrange(match(method, pub_order)) %>%
  filter(method %in% selected_methods) %>%
  mutate(method = 
           pub_names$publication[match(method, pub_names$internal)]) %>%
  mutate(rmse = rmse * 100, abs_bias = abs_bias * 100, 
         int_len = int_len * 100, cov90 = cov90 * 100,
         med_rmse = med_rmse * 100, med_abs_bias = med_abs_bias * 100, 
         med_int_len = med_int_len * 100, med_cov90 = med_cov90 * 100) %>%
  knitr::kable(digits = c(0, 2, 2, 0, 2, 2, 2, 0, 2), format = "latex", booktabs = T,
               linesep = "") %>%
  writeLines("paper/figures/combined_sim_results_table.tex")


pop_dat <- readRDS(paste0(sim_res_dir, "spatial_pop_dat.rds"))
cov_dat <- pop_dat %>%
  select(id_cluster, lon, lat, x1_ns, x1_ns_sig, x1_BYM2, x2_BYM2, x1_SPDE) %>%
  setNames(c("id_cluster", "lon", "lat", "X1", "X2", "X3", "X4", "X5")) %>%
  pivot_longer(cols = c("X1", "X2", "X3", "X4", "X5")) %>%
  rename(Covariate = name, Cov_val = value)
cov_sf <- st_as_sf(cov_dat, coords = c("lon", "lat"))
temp<-ggplot(cov_sf, aes(color = Cov_val)) + geom_sf(size = .08) +
  facet_wrap(~Covariate) + scale_color_viridis_c(option = "inferno", name = "Value")
ggsave(paste0("paper/figures/covariate_maps.png"), temp, dpi = 600, width = 10, height = 6) 
