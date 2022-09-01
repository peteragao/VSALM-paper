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
library(raster)
library(patchwork)
library(rstan)
library(SUMMER)
library(VSALM)

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
  "spatialMeanSmoothUnmatched",
  "spatialJointSmoothUnmatched"
)
#### FILE MANAGEMENT ####
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/BALM-SAE/"))

res_dir <- "results/"
dir.create(file.path(res_dir), showWarnings = FALSE)
cluster_res_dir <- "results/cluster/"
dir.create(file.path(cluster_res_dir), showWarnings = FALSE)
sim_res_dir <- "results/cluster/sims/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)

figs_dir <- "results/figures/"
#### 0.1 COUNTRY INFO ##########################################################
country <- "Nigeria"
survey_year <- 2018
gadm_abbrev <- "NGA"
pop_abbrev <- 'nga'

# dhsStata, which contains survey data
dhs_file <- "dhsStata/NGKR7BDT/NGKR7BFL.DTA"
# survey GPS
dhsFlat_file <- "NGGE7BFL"

country_dir <- paste0("../", country, "/")
poly_path <- paste0(country_dir, "shapeFiles_gadm")
res_dir <- "./results"
proj_data_dir <-  paste0("./data/" , country, "/", sep="")
dir.create(file.path(proj_data_dir), showWarnings = FALSE)
country_res_dir <-  paste0(res_dir , country, "/", sep="")
dir.create(file.path(country_res_dir), showWarnings = FALSE)
# link to the main dropbox folder for population rasters 
pop_raw_dir <- paste0(country_dir, 'Population/')
cov_raw_dir <- paste0(country_dir, 'covariates/')
#### 1 GENERATE POPULATION #####################################################
#### 1.1 LOAD NIGERIA INFO #####################################################
country <- "Nigeria"
survey_year <- 2018
gadm_abbrev <- "NGA"
pop_abbrev <- 'nga'

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
poly_adm2$NAME_1[poly_adm2$NAME_1 == "Federal Capital Territory"] <- "Abuja"
poly_adm1$NAME_1[poly_adm1$NAME_1 == "Federal Capital Territory"] <- "Abuja"

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
if (TRUE) {
  svy_dat <- readRDS(paste0(proj_data_dir, "clean_DHS_data_DELETE-ME.rds"))
} else {
  #### 1.2.1 Load EA locations ####
  ea_locs_path <- paste0(country_dir, "dhsFlat/", dhsFlat_file)
  ea_locs <- readOGR(dsn = path.expand(ea_locs_path),
                     layer = as.character(dhsFlat_file))
  
  # Remove EAs with missing geo information
  missing_idx <- ea_locs$LATNUM == 0
  ea_locs <- ea_locs[!missing_idx, ]
  
  #### 1.2.2 Load Recode data ####
  in_dat = read.dta13(paste0(country_dir, dhs_file))
  
  svy_dat = data.frame(
    cluster = in_dat$v001, 
    hshold = in_dat$v002,
    stratum = in_dat$v023, 
    h9 = in_dat$h9,
    mcv_yes = 1 * (in_dat$h9 == "vaccination date on card" | 
                     in_dat$h9 == "reported by mother" | 
                     in_dat$h9 == "vaccination marked on card"),
    alive = in_dat$b5 == "yes",
    doi = in_dat$v008,
    dob = in_dat$b3,
    wt = in_dat$v005/1000000
  )
  svy_dat <- subset(svy_dat, doi-dob <= 23 & doi-dob >= 12)
  svy_dat <- subset(svy_dat, !is.na(mcv_yes))
  svy_dat = subset(svy_dat, alive)
  
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
          file = paste0(proj_data_dir, "clean_DHS_data_DELETE-ME.rds"))
}
svy_dat$DHSadm1 <- stringr::str_sub(svy_dat$stratum, end = -7)

#### 2 FIGURES #################################################################
#### 2.1 MAP OF NIGERIA ########################################################
set.seed(1204)
nga_map <- ggplot(data = st_as_sf(poly_adm1)) +
  geom_sf(lwd = .08, fill = NA) + 
  geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
  geom_point(data = svy_dat %>%
               mutate(urban = as.factor(ifelse(urban == "U", "Urban", "Rural"))),
             aes(x = lon, y = lat, color = urban),
             shape = 3, alpha = 1, size = .85) +
  #scale_fill_manual(values = sample(colorRampPalette(c("#008753", "#FFFFFF"))(37))) + 
  scale_color_manual(values = c("mediumblue", "gold"), name = NULL) + 
  guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
  theme_bw() + guides(fill="none") +
  theme(legend.position="bottom",
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/nga_map.pdf"), nga_map,
       width = 7, height = 7)
# ggsave(paste0("paper/figures/nga_map.tiff"),
#        nga_map, dpi = 200, width = 7, height = 7)
#### 2.2 MAP OF DIRECT ESTIMATES ###############################################
sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest=T, 
                        weights = ~wt, data=svy_dat)

dir_est <- svyby(~mcv_yes, ~admin1_char, sample_des, svymean)
hajek_est <- dir_est %>%
  rename(mean = mcv_yes) %>%
  mutate(median = mean, var = se ^ 2, domain = admin1_char) %>%
  dplyr::select(domain, mean, median, var) %>%
  mutate(lower = mean + qnorm((1-.9)/2) * sqrt(var),
         upper = mean + qnorm((1+.9)/2) * sqrt(var),
         method = "Hajek")  %>%
  arrange(match(domain, rownames(admin1_mat))) %>%
  left_join(admin1_names, by = c("domain" = "Internal"))
adm1_maps <- st_as_sf(poly_adm1) %>% 
  dplyr::select(NAME_1) %>%
  left_join(hajek_est, by = c("NAME_1" = "GADM"))
hajek_map <- ggplot(adm1_maps, aes(fill = median)) + geom_sf(lwd = 0) +
  scale_fill_viridis_c(direction = -1, name = "MCV")  +
  theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")

ggsave(paste0("../VSALM-paper/paper/figures/nga_mcv_direct_est_map.pdf"), hajek_map,
       width = 7, height = 7)

#### 2.3 MAP OF MODEL ESTIMATES ################################################

all_res <- readRDS("results/Nigeria/Nigeria_MCV_ests.rds") %>%
  left_join(admin1_names, by = c("domain" = "Internal")) %>%
  rename(State = GADM)
pub_order2 <- c("spatialMeanSmoothUnmatched", "spatialJointSmoothUnmatched")
sel_nga_res <- all_res %>%
  filter(method %in% pub_order2) %>%
  mutate(method = 
           pub_names$publication[match(method, pub_names$internal)]) %>%
  mutate(method = factor(method, 
                         levels = c("Spatial Unmatched MS",
                                    "Spatial Unmatched JS"))) %>%
  mutate(example = "Nigeria MCV prev.")



ggplot(all_res, aes(x = direct, y = median, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Nigeria-MCV-ests.pdf", width = 12, height = 6)
ggplot(all_res, aes(x = dir_var, y = var, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Nigeria-MCV-ests-var.pdf", width = 12, height = 6)
ggplot(all_res, aes(x = dir_int_length, y = upper - lower, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Nigeria-MCV-ests-int-length.pdf", width = 12, height = 6)



#### 2.4 TABLE OF MODEL ESTIMATES ##############################################
out_table <- all_res %>%
  filter(method %in% pub_order) %>%
  mutate(interval = paste0("(", round(lower, 2),", ", round(upper, 2), ")")) %>%
  dplyr::select(median, interval, method, State) %>%
  pivot_wider(values_from = c("median", "interval"),
              names_from = method,
              names_glue = "{method}_{.value}")
out_table <- out_table[, c(1, 2, 5,3, 6, 4, 7)]
out_table <- out_table[order(out_table[, 2], decreasing = T),]

out_table %>% 
  knitr::kable(digits = rep(2, 7), format = "latex", booktabs = T,
               linesep = "", col.names = c("State",
                                           rep(c("Point est.", "Int est."), 3))) %>%
  writeLines(paste0("../VSALM-paper/paper/figures/nga_mcv_est_table.tex"))

mcv_lims <- range(all_res$median)
length_lims <- range(all_res$upper - all_res$lower)
adm1_maps <- st_as_sf(poly_adm1) %>% 
  dplyr::select(NAME_1) %>%
  left_join(all_res, by = c("NAME_1" = "State"))
sel_maps <- adm1_maps %>%
  filter(method %in% pub_order) %>%
  mutate(method = 
           pub_names$publication[match(method, pub_names$internal)]) %>%
  mutate(method = factor(method, 
                         levels = c("Direct (Hájek)", "Spatial MS", 
                                    "Spatial Unmatched MS",
                                    "Spatial JS", 
                                    "Spatial Unmatched JS")))

sel_ests <- ggplot(sel_maps,
                   aes(fill = median)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, name = "MCV", limits = mcv_lims)  +
  facet_wrap(~method, nrow = 1) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 20)) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/nga_mcv_est_maps.pdf"), 
       sel_ests, width = 12, height = 3.5)
# ggsave(paste0("paper/figures/sel_adm1_mcv_est_maps.tiff"),
#        sel_ests, dpi = 200, width = 6, height = 7)

#### 2.5 MAPS OF ADMIN-1 CI LENGTHS ############################################
sel_lengths <- ggplot(sel_maps,
                      aes(fill = upper - lower)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, option = "magma",
                       name = "90% Int. length", limits = length_lims)  +
  facet_wrap(~method, nrow = 1) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 20)) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/nga_mcv_int_len_maps.pdf"), 
       sel_lengths, width = 12, height = 3.5)
# ggsave(paste0("paper/figures/sel_adm1_mcv_len_maps.tiff"),
#        sel_lengths, dpi = 200, width = 6, height = 7)

comb_maps <- (sel_ests +
                theme(plot.caption = element_text(hjust = .5 , size = 15))) /
  (sel_lengths +
     theme(plot.caption = element_text(hjust = .5 , size = 15)))
ggsave(paste0("../VSALM-paper/paper/figures/nga_mcv_combined_maps.pdf"), 
       comb_maps, width = 12, height = 7)
# ggsave(paste0("paper/figures/sel_adm1_mcv_maps.tiff"),
#        comb_maps, dpi = 200, width = 12, height = 7)

ht_comp <- adm1_est %>%
  filter(method != "Hájek") %>%
  left_join(adm1_est %>% 
              filter(method == "Hájek") %>% 
              dplyr::select(region, est, upper, lower) %>% 
              rename(hajek = est) %>%
              mutate(intlen = upper - lower) %>%
              dplyr::select(region, hajek, intlen),
            by = "region")

gg <- ggplot(ht_comp %>% 
               filter(method %in% selected_methods), 
             aes(x = hajek, y = est, color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) +
  
  scale_color_discrete(name = "") + 
  theme_classic() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("Hájek estimate") + ylab("Model estimate")
ggsave(plot = gg, height = 5, width = 8,
       filename = paste0("../SMA-SAE/paper/figures/all_adm1_mcv_est_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_est_scatter.tiff"))
gg <- ggplot(ht_comp %>% filter(method %in% selected_methods), 
             aes(x = intlen, y = upper-lower, color = method)) +
  geom_point(shape = 17) + geom_abline(slope = 1) + facet_wrap(~method) +
  scale_color_discrete(name = "") + 
  theme_classic() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("Hájek 90% interval length") + ylab("Model 90% interval length")
ggsave(plot = gg, height = 5, width = 8,
       filename = paste0("../SMA-SAE/paper/figures/all_adm1_mcv_se_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_se_scatter.tiff"))

#### 2.6 TABLE OF MODEL PARAMETERS #############################################
params <- 
  readRDS(file = "results/Nigeria/Nigeria_MCV_param_ests.rds")
data.frame(ms_med = c(as.character(round(params$spatial_ms[,6], 2)), rep("", 4)),
           ms_int = c(paste0("(", round(params$spatial_ms[,4], 2),", ", 
                             round(params$spatial_ms[,8], 2), ")"), rep("", 4)),
           js_med = as.character(round(params$spatial_js[,6], 2)),
           js_int = paste0("(", round(params$spatial_js[,4], 2),", ", 
                           round(params$spatial_js[,8], 2), ")")) %>%
  knitr::kable(format = "latex", booktabs = T,
               linesep = "") %>%
  writeLines(paste0("../VSALM-paper/paper/figures/nga_mcv_param_table.tex"))




#### ---- MALAWI EXAMPLE ---- ##################################################
home_dir <- '~/'
if (!("Dropbox" %in% list.files("~"))) {
  home_dir <- "~/../../mnt/beegfs/homes/petergao/"
}
setwd(paste0(home_dir, "Dropbox/BALM-SAE/"))

res_dir <- "results/"
dir.create(file.path(res_dir), showWarnings = FALSE)
cluster_res_dir <- "results/cluster/"
dir.create(file.path(cluster_res_dir), showWarnings = FALSE)
sim_res_dir <- "results/cluster/sims/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)

figs_dir <- "results/figures/"
country <- "Malawi"
gadm_abbrev <- "MWI"
pop_abbrev <- 'mwi'
country_dir <- paste0("data/", country, "/")
poly_path <- paste0(country_dir, "shapeFiles_gadm")

# dhsStata, which contains survey data

ir_file <- "data/Malawi/MW_DHS/MWIR7HDT/MWIR7HFL.DTA"
hh_file <- "data/Malawi/MW_DHS/MWHR7HDT/MWHR7HFL.DTA"
hiv_file <- "data/Malawi/MW_DHS/MWAR7ADT/MWAR7AFL.DTA"
# survey GPS
dhsFlat_file <- "data/Malawi/MW_DHS/MWGE7AFL"
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
  hh_dat <- read.dta13(paste0(hh_file))
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
#### 1.3.1 Hajek ####
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

#### 2 FIGURES #################################################################
#### 2.1 MAP OF MALAWI ########################################################
set.seed(1204)
mwi_map <- ggplot(data = st_as_sf(poly_adm1)) +
  geom_sf(lwd = .08, fill = NA) + 
  geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
  geom_point(data = svy_dat %>%
               mutate(urban = as.factor(ifelse(urban == "U", "Urban", "Rural"))),
             aes(x = lon, y = lat, color = urban),
             shape = 3, alpha = 1, size = .85) +
  #scale_fill_manual(values = sample(colorRampPalette(c("#008753", "#FFFFFF"))(37))) + 
  scale_color_manual(values = c("mediumblue", "gold"), name = NULL) + 
  guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
  theme_bw() + guides(fill="none") +
  theme(legend.position="none",
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/mwi_map.pdf"), mwi_map,
       width = 4, height = 7)
# ggsave(paste0("paper/figures/nga_map.tiff"),
#        nga_map, dpi = 200, width = 7, height = 7)
#### 2.2 MAP OF DIRECT ESTIMATES ###############################################
hajek_est <- hajek_est %>%
  left_join(admin1_names, by = c("domain" = "Internal"))
adm1_maps <- st_as_sf(poly_adm1) %>% 
  dplyr::select(NAME_1) %>%
  left_join(hajek_est, by = c("NAME_1" = "GADM"))
hajek_map <- ggplot(adm1_maps, aes(fill = median)) + geom_sf(lwd = 0) +
  scale_fill_viridis_c(direction = -1, name = "HIV prev.", option = "plasma")  +
  theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("") + ylab("")

ggsave(paste0("../VSALM-paper/paper/figures/mwi_hiv_direct_est_map.pdf"), hajek_map,
       width = 4.5, height = 7)


#### 2.3 MAP OF MODEL ESTIMATES ################################################

all_res <- readRDS("../BALM-SAE/results/Malawi/Malawi_HIV_ests.rds") %>%
  left_join(admin1_names, by = c("domain" = "Internal")) %>%
  rename(State = GADM)
pub_order2 <- c("spatialMeanSmoothUnmatched", "spatialJointSmoothUnmatched")
sel_mwi_res <- all_res %>%
  filter(method %in% pub_order2) %>%
  mutate(method = 
           pub_names$publication[match(method, pub_names$internal)]) %>%
  mutate(method = factor(method, 
                         levels = c("Spatial Unmatched MS",
                                    "Spatial Unmatched JS"))) %>%
  mutate(example = "Malawi HIV prev.")



ggplot(all_res, aes(x = direct, y = median, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Malawi-HIV-ests.pdf", width = 12, height = 6)
ggplot(all_res, aes(x = dir_var, y = var, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Malawi-HIV-ests-var.pdf", width = 12, height = 6)
ggplot(all_res, aes(x = dir_int_length, y = upper - lower, color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~method)
ggsave("../VSALM-paper/paper/figures/Malawi-HIV-ests-int-length.pdf", width = 12, height = 6)



#### 2.4 TABLE OF MODEL ESTIMATES ##############################################
out_table <- all_res %>%
  filter(method %in% pub_order) %>%
  mutate(interval = paste0("(", round(lower, 2),", ", round(upper, 2), ")")) %>%
  dplyr::select(median, interval, method, State) %>%
  pivot_wider(values_from = c("median", "interval"),
              names_from = method,
              names_glue = "{method}_{.value}")
out_table <- out_table[, c(1, 2, 5,3, 6, 4, 7)]
out_table <- out_table[order(out_table[, 2], decreasing = T),]

out_table %>% 
  knitr::kable(digits = rep(2, 7), format = "latex", booktabs = T,
               linesep = "", col.names = c("State",
                                           rep(c("Point est.", "Int est."), 3))) %>%
  writeLines(paste0("../VSALM-paper/paper/figures/mwi_hiv_est_table.tex"))

mcv_lims <- range(all_res$median)
length_lims <- range(all_res$upper - all_res$lower)
adm1_maps <- st_as_sf(poly_adm1) %>% 
  dplyr::select(NAME_1) %>%
  left_join(all_res, by = c("NAME_1" = "State"))
sel_maps <- adm1_maps %>%
  filter(method %in% pub_order) %>%
  mutate(method = 
           pub_names$publication[match(method, pub_names$internal)]) %>%
  mutate(method = factor(method, 
                         levels = c("Direct (Hájek)", "Spatial MS", 
                                    "Spatial Unmatched MS",
                                    "Spatial JS", 
                                    "Spatial Unmatched JS")))

sel_ests <- ggplot(sel_maps,
                   aes(fill = median)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, name = "HIV prev.", limits = mcv_lims)  +
  facet_wrap(~method, nrow = 1) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 9)) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/mwi_hiv_est_maps.pdf"), 
       sel_ests, width = 8, height = 5)
# ggsave(paste0("paper/figures/sel_adm1_mcv_est_maps.tiff"),
#        sel_ests, dpi = 200, width = 6, height = 7)

#### 2.5 MAPS OF ADMIN-1 CI LENGTHS ############################################
sel_lengths <- ggplot(sel_maps,
                      aes(fill = upper - lower)) + geom_sf(lwd = 0) + 
  scale_fill_viridis_c(direction = -1, option = "magma",
                       name = "90% Int. length", limits = length_lims)  +
  facet_wrap(~method, nrow = 1) + theme_bw() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 9)) + 
  xlab("") + ylab("")
ggsave(paste0("../VSALM-paper/paper/figures/mwi_hiv_int_len_maps.pdf"), 
       sel_lengths, width = 8, height = 5)
# ggsave(paste0("paper/figures/sel_adm1_mcv_len_maps.tiff"),
#        sel_lengths, dpi = 200, width = 6, height = 7)

comb_maps <- (sel_ests +
  theme(plot.caption = element_text(hjust = .5 , size = 11))) /
  (sel_lengths +
  theme(plot.caption = element_text(hjust = .5 , size = 11)))
ggsave(paste0("../VSALM-paper/paper/figures/mwi_hiv_combined_maps.pdf"), 
       comb_maps, width = 6, height = 8)
# ggsave(paste0("paper/figures/sel_adm1_mcv_maps.tiff"),
#        comb_maps, dpi = 200, width = 12, height = 7)

ht_comp <- adm1_est %>%
  filter(method != "Hájek") %>%
  left_join(adm1_est %>% 
              filter(method == "Hájek") %>% 
              dplyr::select(region, est, upper, lower) %>% 
              rename(hajek = est) %>%
              mutate(intlen = upper - lower) %>%
              dplyr::select(region, hajek, intlen),
            by = "region")

gg <- ggplot(ht_comp %>% 
               filter(method %in% selected_methods), 
             aes(x = hajek, y = est, color = method)) +
  geom_point() + geom_abline(slope = 1) + facet_wrap(~method) +
  
  scale_color_discrete(name = "") + 
  theme_classic() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("Hájek estimate") + ylab("Model estimate")
ggsave(plot = gg, height = 5, width = 8,
       filename = paste0("../SMA-SAE/paper/figures/all_adm1_mcv_est_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_est_scatter.tiff"))
gg <- ggplot(ht_comp %>% filter(method %in% selected_methods), 
             aes(x = intlen, y = upper-lower, color = method)) +
  geom_point(shape = 17) + geom_abline(slope = 1) + facet_wrap(~method) +
  scale_color_discrete(name = "") + 
  theme_classic() + 
  theme(strip.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank()) + 
  xlab("Hájek 90% interval length") + ylab("Model 90% interval length")
ggsave(plot = gg, height = 5, width = 8,
       filename = paste0("../SMA-SAE/paper/figures/all_adm1_mcv_se_scatter.pdf"))
# ggsave(plot = gg, dpi = 200, height = 6, width = 7,
#        filename = paste0("paper/figures/all_adm1_mcv_se_scatter.tiff"))

int_length_dat <- bind_rows(sel_nga_res, sel_mwi_res)
int_length_comp <- ggplot(int_length_dat, aes(x = dir_int_length, y = upper - lower,
                           color = method)) + geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~example, scales = 'free') + 
  xlab("Hájek 90% interval length") + ylab("Model-based 90% interval length") + 
  theme_bw() + guides(fill="none") +
  scale_color_discrete(name = "") + 
  theme(legend.position="bottom",
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggsave(paste0("../VSALM-paper/paper/figures/int_length_scatter.pdf"), 
       int_length_comp, width = 8, height = 4)
ggplot(int_length_dat, aes(x = direct, y = median, color = method)) + geom_point() +
  geom_abline(slope = 1)  +
  facet_wrap(~example, scales = 'free')



params <- 
  readRDS(file = "results/Malawi/Malawi_HIV_param_ests.rds")

data.frame(ms_med = c(as.character(round(params$spatial_ms[,6], 2)), rep("", 4)),
           ms_int = c(paste0("(", round(params$spatial_ms[,4], 2),", ", 
                             round(params$spatial_ms[,8], 2), ")"), rep("", 4)),
           js_med = as.character(round(params$spatial_js[,6], 2)),
           js_int = paste0("(", round(params$spatial_js[,4], 2),", ", 
                           round(params$spatial_js[,8], 2), ")")) %>%
  knitr::kable(format = "latex", booktabs = T,
               linesep = "") %>%
  writeLines(paste0("../VSALM-paper/paper/figures/mwi_hiv_param_table.tex"))

