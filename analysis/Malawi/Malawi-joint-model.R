#### 0 LIBRARIES AND PATHS ####################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(INLA)
library(survey)
library(purrr)
library(rgdal)
#library(geosphere)
library(raster)
library(patchwork)
library(rstan)
library(SUMMER)
library(readstata13)

library(VSALM)
#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/VSALM-paper/"))

res_dir <- "results/"
dir.create(file.path(res_dir), showWarnings = FALSE)
cluster_res_dir <- "results/cluster/"
dir.create(file.path(cluster_res_dir), showWarnings = FALSE)
sim_res_dir <- "results/cluster/sims/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)

figs_dir <- "results/figures/"
#### 0.1 COUNTRY INFO ##########################################################
country <- "Malawi"
gadm_abbrev <- "MWI"
pop_abbrev <- 'mwi'
country_dir <- paste0("data/", country, "/")
poly_path <- paste0(country_dir, "shapeFiles_gadm")

# dhsStata, which contains survey data

ir_file <- "data/Malawi/MW_DHS/MWIR7HDT/MWIR7HFL.DTA"
hiv_file <- "data/Malawi/MW_DHS/MWAR7ADT/MWAR7AFL.DTA"
# survey GPS
dhsFlat_file <- "data/Malawi/MW_DHS/MWGE7AFL"
#### 1 GENERATE POPULATION #####################################################
#### 1.1 LOAD MALAWI INFO #####################################################

poly_path <- paste0("data/", country, "/shapeFiles_gadm/")
#### 1.1.1 Load GADM polygons ####
poly_layer_adm0 <- paste('gadm36', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm36', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm36', gadm_abbrev,
                         '2', sep = "_")

poly_adm0 <- readOGR(dsn = poly_path, encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly_layer_adm0)) 
# use encoding to read special characters
poly_adm1 <- readOGR(dsn = poly_path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly_layer_adm1))

if(sum(grepl(paste('gadm36', gadm_abbrev,
                   '2', sep = "_"), list.files(poly_path))) != 0){
  poly_adm2 <- readOGR(dsn = poly_path,encoding = "UTF-8", use_iconv = TRUE,
                       layer = as.character(poly_layer_adm2))}

if (exists("poly_adm2")) {
  proj4string(poly_adm0) <- proj4string(poly_adm1)  <- proj4string(poly_adm2)
}else {
  proj4string(poly_adm0) <- proj4string(poly_adm1)
}
poly_adm2$NAME_1 <- as.character(poly_adm2$NAME_1)
poly_adm1$NAME_1 <- as.character(poly_adm1$NAME_1)
#### 1.1.2 Create adjacency matrices and tables ####
if(exists("poly_adm1")){
  admin1_mat <- poly2nb(SpatialPolygons(poly_adm1@polygons))
  admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)
  colnames(admin1_mat) <- 
    rownames(admin1_mat) <-
    paste0("admin1_", 1:dim(admin1_mat)[1])
  admin1_names <- data.frame(GADM = poly_adm1@data$NAME_1,
                             Internal = rownames(admin1_mat))
} else {
  message("There is no Admin1 polygon file.")
}
if(exists("poly_adm2")){
  admin2_mat <- poly2nb(SpatialPolygons(poly_adm2@polygons))
  admin2_mat <- nb2mat(admin2_mat, zero.policy = TRUE)
  colnames(admin2_mat) <- 
    rownames(admin2_mat) <- 
    paste0("admin2_", 1:dim(admin2_mat)[1])
  admin2_names <- data.frame(GADM = poly_adm2@data$NAME_2,
                             Internal = rownames(admin2_mat))
} else {
  message("There is no Admin2 polygon file.")
}

#### 1.2 DHS DATA ##############################################################
if (T) {
  svy_dat <- readRDS("data/Malawi/clean_DHS_data.rds")
} else {
  #### 1.2.1 Load EA locations ####
  ea_locs_path <- "data/Malawi/MW_DHS/MWGE7AFL/"
  ea_locs <- readOGR(dsn = path.expand(ea_locs_path),
                     layer = "MWGE7AFL")
  
  # Remove EAs with missing geo information
  missing_idx <- ea_locs$LATNUM == 0
  ea_locs <- ea_locs[!missing_idx, ]
  
  #### 1.2.2 Load Recode data #h##
  ir_dat <- read.dta13(paste0(ir_file))
  hiv_dat <- read.dta13(paste0(hiv_file))
  
  svy_dat <- data.frame(
    cluster = hiv_dat$hivclust, 
    hshold = hiv_dat$hivnumb,
    line = hiv_dat$hivline,
    wt = hiv_dat$hiv05/1000000,
    hiv = 1 * (hiv_dat$hiv06 == "positive")
  )
  join_dat <- data.frame(
    cluster = ir_dat$v001,
    hshold = ir_dat$v002,
    line = ir_dat$v003,
    stratum = ir_dat$v023
  )
  svy_dat <- svy_dat %>%
    merge(join_dat, by = c("cluster", "hshold", "line"))
  
  #### 1.2.3 Merge geographic info ####
  ea_dat <- data.frame(cluster = ea_locs$DHSCLUST,
                       urban = ea_locs$URBAN_RURA,
                       lon = ea_locs$LONGNUM,
                       lat = ea_locs$LATNUM)
  svy_dat = merge(svy_dat, ea_dat, by = "cluster")
  
  points_frame <- as.data.frame(svy_dat[,c("lon", "lat")])
  points_frame <- SpatialPoints(points_frame)
  
  # assign points to admin 2
  poly_over_adm2 <- SpatialPolygons(poly_adm2@polygons)
  proj4string(points_frame) <-
    proj4string(poly_over_adm2) <- 
    proj4string(poly_adm2)  <- 
    proj4string(poly_adm1)  
  admin2_key <- over(points_frame, poly_over_adm2)
  miss_frame_adm2 <- 
    matrix(unique(points_frame@coords[which(is.na(admin2_key)),]), ncol = 2)
  
  if(dim(miss_frame_adm2)[1] != 0){
    miss_poly_adm2 <- dist2Line( miss_frame_adm2, poly_over_adm2)
    for(i in 1:dim(miss_poly_adm2)[1]){
      long_ids <- 
        which(points_frame@coords[,c("lon")] %in% miss_frame_adm2[i,1])
      lat_ids <-
        which(points_frame@coords[,c("lat")] %in% miss_frame_adm2[i,2])
      ids <- intersect(long_ids, lat_ids)
      admin2_key[ids] <- rep(miss_poly_adm2[i, 'ID'], length(ids))
    }
  }
  svy_dat$admin2 <- admin2_key
  svy_dat$admin2_char <- paste0("admin2_", admin2_key)
  svy_dat$admin2_name <- as.character(poly_adm2@data$NAME_2)[admin2_key]
  
  # assign points to admin 2
  poly_over_adm1 <- SpatialPolygons(poly_adm1@polygons)
  proj4string(points_frame) <-
    proj4string(poly_over_adm1) <- 
    proj4string(poly_adm1) 
  admin1_key <- over(points_frame, poly_over_adm1)
  miss_frame_adm1 <- 
    matrix(unique(points_frame@coords[which(is.na(admin1_key)),]), ncol = 2)
  if(dim(miss_frame_adm1)[1] != 0){
    miss_poly_adm1 <- dist2Line(miss_frame_adm1, poly_over_adm1)
    for(i in 1:dim(miss_poly_adm1)[1]){
      long_ids <- 
        which(points_frame@coords[,c("lon")] %in% miss_frame_adm1[i,1])
      lat_ids <- 
        which(points_frame@coords[,c("lat")] %in% miss_frame_adm1[i,2])
      ids <- intersect(long_ids, lat_ids)
      admin1_key[ids] <- rep(miss_poly_adm1[i, 'ID'], length(ids))
    }
  }
  svy_dat$admin1 <- admin1_key
  svy_dat$admin1_char <- paste0("admin1_", admin1_key)
  svy_dat$admin1_name <- as.character(poly_adm1@data$NAME_1)[admin1_key]
  
  # add number of trials
  svy_dat$n_trials <- 1
  saveRDS(svy_dat,
          file = paste0("data/Malawi/clean_DHS_data.rds"))
}
#svy_dat$DHSadm1 <- stringr::str_sub(svy_dat$stratum, end = -7)
#### 2 DESIGN-BASED METHODS ####################################################

#### 2.1.1 Hajek ####
sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest=T, 
                        weights = ~wt, data=svy_dat)

dir_est <- svyby(~hiv, ~admin1_char, sample_des, svymean)
hajek_est <- dir_est %>%
  rename(mean = hiv) %>%
  mutate(median = mean, var = se ^ 2, domain = admin1_char) %>%
  dplyr::select(domain, mean, median, var) %>%
  mutate(lower = mean + qnorm((1-.9)/2) * sqrt(var),
         upper = mean + qnorm((1+.9)/2) * sqrt(var),
         method = "Hajek")  %>%
  arrange(match(domain, rownames(admin1_mat)))

nc_dat <- svy_dat %>%
  group_by(admin1_char) %>%
  summarize(na = n(),
            df = length(unique(cluster)) - 1,
            c_size = n() / length(unique(cluster))) %>%
  arrange(match(admin1_char, rownames(admin1_mat)))

admin1_mat <- admin1_mat[-10, -10]
hajek_est <- hajek_est[-10,]
nc_dat <- nc_dat[-10,]
iid_mean_smooth_est <-
  VSALM::iidMeanSmooth(Yhat = hajek_est$mean,
                       Vhat = hajek_est$var,
                       domain = hajek_est$domain)

iid_mean_smooth_logit_est <-
  VSALM::iidMeanSmoothLogit(Yhat = hajek_est$mean,
                            Vhat = hajek_est$var,
                            domain = hajek_est$domain)
iid_mean_smooth_unmatched_est <-
  VSALM::iidMeanSmoothUnmatched(Yhat = hajek_est$mean,
                                Vhat = hajek_est$var,
                                domain = hajek_est$domain)
spatial_mean_smooth_est <-
  VSALM::spatialMeanSmooth(Yhat = hajek_est$mean,
                           Vhat = hajek_est$var,
                           domain = hajek_est$domain,
                           adj_mat = admin1_mat)
spatial_mean_smooth_logit_est <-
  VSALM::spatialMeanSmoothLogit(Yhat = hajek_est$mean,
                                Vhat = hajek_est$var,
                                domain = hajek_est$domain,
                                adj_mat = admin1_mat)
spatial_mean_smooth_unmatched_res <-
  VSALM::spatialMeanSmoothUnmatched(Yhat = hajek_est$mean,
                                    Vhat = hajek_est$var,
                                    domain = hajek_est$domain,
                                    adj_mat = admin1_mat, 
                                    detailed_output = T)
spatial_mean_smooth_unmatched_est <- spatial_mean_smooth_unmatched_res$est

spatial_mean_smooth_unmatched_stan <- 
  summary(spatial_mean_smooth_unmatched_res$stan)$summary
spatial_mean_smooth_unmatched_param <- 
  spatial_mean_smooth_unmatched_stan[c("mu", "sigma_u", "phi"),]
ijs_rhat <- Inf
nrun <- 0
while (ijs_rhat >= 1.01 & nrun < 10) {
  iid_joint_smooth_res <-
    VSALM::iidJointSmooth(Yhat = hajek_est$mean,
                          Vhat = hajek_est$var,
                          domain = hajek_est$domain,
                          na = nc_dat$na,
                          df = nc_dat$df, detailed_output = T)
  temp <- summary(iid_joint_smooth_res$stan)$summary
  if (max(temp[, 10]) < 1.01) {
    ijs_rhat <- max(temp[, 10])
    iid_joint_smooth_est <- iid_joint_smooth_res$est
  }
  nrun <- nrun + 1
}
iid_joint_smooth_logit_est <-
  VSALM::iidJointSmoothLogit(Yhat = hajek_est$mean,
                             Vhat = hajek_est$var,
                             domain = hajek_est$domain,
                             na = nc_dat$na,
                             df = nc_dat$df)
iid_joint_smooth_unmatched_est <-
  VSALM::iidJointSmoothUnmatched(Yhat = hajek_est$mean,
                                 Vhat = hajek_est$var,
                                 domain = hajek_est$domain,
                                 na = nc_dat$na,
                                 df = nc_dat$df)
sjs_rhat <- Inf
nrun <- 0
while (sjs_rhat >= 1.01 & nrun < 10) {
  spatial_joint_smooth_res <-
    VSALM::spatialJointSmooth(Yhat = hajek_est$mean,
                              Vhat = hajek_est$var,
                              domain = hajek_est$domain,
                              na = nc_dat$na,
                              df = nc_dat$df,
                              adj_mat = admin1_mat,detailed_output = T)
  temp <- summary(spatial_joint_smooth_res$stan)$summary
  if (max(temp[, 10]) < 1.01) {
    sjs_rhat <- max(temp[, 10])
    spatial_joint_smooth_est <- spatial_joint_smooth_res$est
  }
  nrun <- nrun + 1
}
spatial_joint_smooth_logit_est <-
  VSALM::spatialJointSmoothLogit(Yhat = hajek_est$mean,
                                 Vhat = hajek_est$var,
                                 domain = hajek_est$domain,
                                 na = nc_dat$na,
                                 df = nc_dat$df,
                                 adj_mat = admin1_mat)
spatial_joint_smooth_unmatched_res <-
  VSALM::spatialJointSmoothUnmatched(Yhat = hajek_est$mean,
                                     Vhat = hajek_est$var,
                                     domain = hajek_est$domain,
                                     na = nc_dat$na,
                                     df = nc_dat$df,
                                     adj_mat = admin1_mat, 
                                     detailed_output = T)
spatial_joint_smooth_unmatched_est <- 
  spatial_joint_smooth_unmatched_res$est
spatial_joint_smooth_unmatched_stan <- 
  summary( spatial_joint_smooth_unmatched_res$stan)$summary
spatial_joint_smooth_unmatched_param <- 
  spatial_joint_smooth_unmatched_stan[c("mu", "sigma_u", "phi", "g0", "g1", "g2", "sigma_tau"),]
all_res <- bind_rows(hajek_est,
                     iid_mean_smooth_est,
                     iid_mean_smooth_logit_est,
                     iid_mean_smooth_unmatched_est,
                     iid_joint_smooth_est,
                     iid_joint_smooth_logit_est,
                     iid_joint_smooth_unmatched_est,
                     spatial_mean_smooth_est,
                     spatial_mean_smooth_logit_est,
                     spatial_mean_smooth_unmatched_est,
                     spatial_joint_smooth_est,
                     spatial_joint_smooth_logit_est,
                     spatial_joint_smooth_unmatched_est) %>%
  left_join(hajek_est %>%
              rename(direct = mean, dir_var = var) %>%
              mutate(dir_int_length = upper - lower) %>%
              dplyr::select(domain, direct, dir_var, dir_int_length))

dir.create(file.path("results/Malawi/"), showWarnings = FALSE)
param_list <- list(spatial_ms = spatial_mean_smooth_unmatched_param,
                   spatial_js = spatial_joint_smooth_unmatched_param)
saveRDS(all_res, file = "results/Malawi/Malawi_HIV_ests.rds")
saveRDS(param_list, file = "results/Malawi/Malawi_HIV_param_ests.rds")
all_res <- readRDS("../VSALM-paper/results/Malawi/Malawi_HIV_ests.rds")

