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
library(rstan)
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
sim_res_dir <- "results/cluster/sims/area-level-models/"
dir.create(file.path(sim_res_dir), showWarnings = FALSE)
#### 1 GENERATE POPULATION #####################################################

#### 1.1 LOAD NIGERIA INFO #####################################################
country <- "Nigeria"
survey_year <- 2018
gadm_abbrev <- "NGA"
pop_abbrev <- 'nga'

country_dir <- paste0("../", country, "/")
poly_path <- paste0(country_dir, "shapeFiles_gadm")
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
}else{
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
}else{
  message("There is no Admin2 polygon file.")
}

#### 1.2 SAMPLE EA LOCATIONS ###################################################

gen_ea_locs <- function(N = 73000) {
  C <- 73
  pop_dat_list <- list()
  for (i in 1:37) {
    print(i)
    pix_dat_i <- 
      readRDS(paste0("../data-SMA-SAE/Nigeria/nga_pix_tbl_admin1_", i, ".rds")) %>%
      filter(!is.na(access) & !is.na(poverty))
    pop_dat_list[[i]] <- pix_dat_i %>%
      group_by(urban) %>%   
      nest() %>%
      mutate(ns_per_area = min(nrow(data), N / C)) %>% 
      ungroup() %>% 
      mutate(pop = 
               map2(data, ns_per_area,
                    function(x, y) 
                      slice_sample(x, n = y,
                                   weight_by = cen_pop))) %>%
      dplyr::select(-data, -ns_per_area) %>%
      unnest(pop)
  }
  return(do.call(rbind, pop_dat_list))
}

#' Simulate from Besag model
#'
#' @param Amat radjacency matrix
#' @param mean mean of effect
#' @param marginal_var variance parameter (following Riebler et al)
#'
#' @return data frame of effects for each row in Amat
#' @export
#'
#' @examples
simulate_besag <-
  function(Amat, mean = 0, marginal_var = 1) {
    n <- nrow(Amat)
    
    # build (singular) precision matrix
    Q <- -1 * (Amat != 0)
    diag(Q) <- -rowSums(Q)
    Q <- inla.scale.model(Q, constr = list(A = matrix(1, nrow = 1, ncol = n), e = 0))
    Q <- Q / sqrt(marginal_var)
    Q <- as.matrix(Q)
    
    eigen_Q <- eigen(Q)
    which_to_sim <- which(eigen_Q$values > .0001)
    y <- rnorm(length(which_to_sim), sd = sqrt(1 / eigen_Q$values[which_to_sim]))
    out_sample <- eigen_Q$vectors[, which_to_sim] %*% y
    out <- data.frame(region = 1:n,
                      effect = out_sample + mean,
                      stringsAsFactors = F)
    out
  }

#' Simulate from BYM2 model
#'
#' @param Amat radjacency matrix
#' @param mean mean of effect
#' @param marginal_var variance parameter (following Riebler et al)
#' @param phi spatial proportion parameter (following Riebler et al)
#'
#' @return data frame of effects for each row in Amat
#' @export
#'
#' @examples
simulate_bym2 <-
  function(Amat, mean = 0, marginal_var = 1, phi = 1) {
    n <- nrow(Amat)
    out <- simulate_besag(Amat, mean = mean, 
                          marginal_var = marginal_var) %>%
      mutate(effect = sqrt(phi) * effect + 
               sqrt(1 - phi) * rnorm(n) + mean)
    out
  }

#' Title
#'
#' @param locs 
#' @param mesh 
#' @param eff_range 
#' @param marginal_var 
#' @param max_edge 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
simulate_SPDE = function(locs, mesh = NULL,
                         eff_range = (max(locs[,1]) - min(locs[,1])) / 3,
                         marginal_var = 1, max_edge = eff_range / 5, ...) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    obs_bnd = inla.nonconvex.hull(locs, ...)
    mesh = inla.mesh.2d(boundary = obs_bnd, max.edge = max_edge, ...)
  }
  # calculate SPDE model parameters based on Lindgren Rue (2015)
  # "Bayesian Spatial Modelling with R-INLA"
  # meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 
  # then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  # from page 5 of the paper listed above:
  log_kappa = 0.5 * log(8)
  log_tau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - log_kappa
  theta = c(log(sqrt(marginal_var)), log(eff_range))
  spde <- inla.spde2.matern(mesh, 
                            B.tau = cbind(log_tau, -1, +1),
                            B.kappa = cbind(log_kappa, 0, -1),
                            theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, locs)
  # generate simulations
  sim_gmrf = inla.qsample(1, Q)
  sim_obs = as.matrix(A %*% sim_gmrf)
  as.vector(sim_obs)
}

gen_pop_X <- function(ea_locs) {
  ea_locs_sf <- st_as_sf(ea_locs, coords = c("lon", "lat"))
  pop_dat <- ea_locs
  N <- nrow(pop_dat)
  A <- length(unique(pop_dat$admin1_name))
  
  #### non-spatial with areal variation ####
  pop_dat$x1_ns <-
    rbinom(N, size = 1, prob = .5)
  pop_dat$x1_ns_sig <-
    rbinom(N, size = 1, prob = 0.3 + 0.5 *
             match(pop_dat$admin1_name, admin1_names$GADM) / A)
  
  #### BYM2 covariates ####
  x1_BYM2 <- simulate_bym2(admin1_mat) 
  pop_dat$x1_BYM2 <- 
    x1_BYM2$effect[match(pop_dat$admin1_name, admin1_names$GADM)]
  x2_BYM2 <- simulate_bym2(admin2_mat) 
  pop_dat$x2_BYM2 <- 
    x2_BYM2$effect[match(pop_dat$admin2_name, admin2_names$GADM)]
  #### SPDE covariates ####
  # make mesh
  nc_hull <-
    inla.nonconvex.hull(st_coordinates(ea_locs_sf))
  mesh <- inla.mesh.2d(boundary = nc_hull, max.edge = c(0.75, 2))
  pop_dat$x1_SPDE <- simulate_SPDE(st_coordinates(ea_locs_sf), mesh = mesh) 
  pop_dat$x2_SPDE <- simulate_SPDE(st_coordinates(ea_locs_sf), mesh = mesh) 
  
  pop_dat$admin1_char <- 
    admin1_names$Internal[match(pop_dat$admin1_name, admin1_names$GADM)]
  pop_dat$admin2_char <- 
    admin2_names$Internal[match(pop_dat$admin2_name, admin2_names$GADM)]
  pop_dat$stratum <- 
    paste0(pop_dat$admin1_name, "-", ifelse(pop_dat$urban, "urban", "rural"))
  pop_dat$id_stratum <- match(pop_dat$stratum, unique(pop_dat$stratum))
  pop_dat$id_cluster <- 1:N
  pop_dat$cl_size <- rpois(N, 15)
  return(pop_dat)
}

#### 1.2 GENERATE X ############################################################

gen_pop_Y <- function(pop_dat, eta_a_sd = .25,
                      e_ac_sd = .5,
                      beta = c(SUMMER::logit(.5), .25, -.25,  .5, .25, .25, 0, 0, 0)) {
  N <- nrow(pop_dat)
  pop_dat$l1a = log(1 + pop_dat$access)
  pop_X <- model.matrix(~ x1_ns + x1_ns_sig + x1_BYM2 + x2_BYM2 + 
                          x1_SPDE + x2_SPDE + l1a + poverty, pop_dat)
  pop_dat <- pop_dat %>%
    group_by(admin1_name) %>%
    mutate(eta_a = rep(rnorm(1, sd = eta_a_sd), each = n())) %>%
    ungroup() %>%
    mutate(e_ac = rnorm(N, sd = e_ac_sd),
           p = SUMMER::expit(pop_X %*% beta + eta_a + e_ac))
  pop_dat <- pop_dat[rep(seq_len(nrow(pop_dat)), times = pop_dat$cl_size), ]
  pop_dat$y <- rbinom(nrow(pop_dat), 1, prob = pop_dat$p)
  return(pop_dat)
}

gen_sample_ind <- 
  function(
    pop_dat, strata,
    n_cl_vec = 
      rep(10, length(unique(as.vector(model.matrix(strata, pop_dat)[, 2]))))
  ) {
    strata_vec <- as.vector(model.matrix(strata, pop_dat)[, 2])
    N_S <- length(unique(strata_vec))
    pop_dat$n_cl <- n_cl_vec[match(pop_dat$id_stratum, 1:N_S)]
    pop_dat <- pop_dat %>%
      group_by(id_stratum) %>%
      mutate(pop_s = sum(cl_size)) %>%
      ungroup() %>%
      group_by(id_stratum, n_cl, cl_size) %>%
      mutate(pi = n_cl * cl_size / pop_s, 
             wt = 1 / pi) %>%
      ungroup() 
    sample_dat <- pop_dat %>%
      group_by(id_stratum, n_cl) %>%
      nest() %>% 
      mutate(samp = 
               map(data,
                   function(x) 
                     slice_sample(x, n = n_cl,
                                  weight_by = pi, replace = F)))  %>%
      dplyr::select(-data) %>%
      unnest(samp) %>%
      ungroup()

    return(list(ind = sample_dat$id_cluster, pop_dat = pop_dat))
  }
gen_sample_ind_inf <- 
  function(
    pop_dat, strata,
    prop_sample_per_strata = .02
  ) {
    strata_vec <- as.vector(model.matrix(strata, pop_dat)[, 2])
    n_strata <- length(unique(strata_vec))
    pop_dat$strata <- strata_vec
    pop_dat$idx <- 1:nrow(pop_dat)
    ns_per_area_table <- pop_dat %>%
      group_by(strata) %>%
      summarize(n = round(n() * prop_sample_per_strata))
    ns_per_area <- ns_per_area_table$n
    pop_dat <- pop_dat %>%
      mutate(a = 1 + 1 * exp(x2_SPDE) + 1 * exp(x1_BYM2)) %>%
      left_join(ns_per_area_table, by = "strata") %>%
      group_by(strata) %>%  
      mutate(wt = sum(a) / a / n) %>%
      ungroup()
    sample_dat <- pop_dat %>%
      group_by(strata) %>%  
      nest() %>%             
      ungroup() %>% 
      mutate(samp = 
               map2(data, ns_per_area,
                    function(x, y) 
                      slice_sample(x, n = y,
                                   weight_by = 1/x$wt))) %>%
      dplyr::select(-data) %>%
      unnest(samp)
    return(list(ind = sample_dat$idx, pop_dat = pop_dat))
  }

#### RUN SIMULATIONS ####
run_one_sim <- function(sample_ind, pop_dat_X, pop_gen_fn, i, ...) {
  pop_dat_Y <- pop_gen_fn(pop_dat_X, ...)
  #pop_dat_Y <- gen_pop_Y(pop_dat)
  pop_dat_Y <- pop_dat_Y %>%
    mutate(domain = as.character(admin1_char)) 
  sample_dat <- pop_dat_Y[pop_dat_Y$id_cluster %in% sample_ind, ]
  pop_means <- pop_dat_Y %>%
    group_by(domain) %>%
    summarize(pop_mean = mean(y)) %>%
    as.data.frame() %>%
    arrange(match(domain, rownames(admin1_mat)))
  sample_dat$id <- 1:nrow(sample_dat)
  sample_des <- svydesign(ids = ~id_cluster, strata = ~id_stratum,
                          weights = ~wt, data = sample_dat)
  
  dir_est <- svyby(~y, ~domain, sample_des, svymean)
  hajek_est <- dir_est %>%
    rename(mean = y) %>%
    mutate(median = mean, var = se ^ 2) %>%
    dplyr::select(domain, mean, median, var) %>%
    mutate(lower = mean + qnorm((1-.9)/2) * sqrt(var),
           upper = mean + qnorm((1+.9)/2) * sqrt(var),
           method = "Hajek")  %>%
    arrange(match(domain, rownames(admin1_mat)))
  nc_dat <- sample_dat %>%
    group_by(domain) %>% 
    summarize(na = n(),
              df = length(unique(id_cluster)) - 1,
              c_size = n() / length(unique(id_cluster))) %>%
    arrange(match(domain, rownames(admin1_mat)))
  print(hajek_est)
  temp <- data.frame(
    domain = hajek_est$domain,
    Yhat = hajek_est$mean,
    Vhat = hajek_est$var
  )
  smoth <- temp$Vhat > 1e-5
  print(smoth)
  iid_mean_smooth_est <-
    VSALM::iidMeanSmooth(Yhat = hajek_est$mean,
                         Vhat = hajek_est$var,
                         domain = hajek_est$domain,
                         verbose = T)
  print(iid_mean_smooth_est)
  iid_mean_smooth_logit_est <-
    VSALM::iidMeanSmoothLogit(Yhat = 
                                hajek_est$mean,
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
  spatial_mean_smooth_unmatched_est <-
    VSALM::spatialMeanSmoothUnmatched(Yhat = hajek_est$mean,
                                      Vhat = hajek_est$var,
                                      domain = hajek_est$domain,
                                      adj_mat = admin1_mat)
  
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
  spatial_joint_smooth_unmatched_est <-
    VSALM::spatialJointSmoothUnmatched(Yhat = hajek_est$mean,
                                       Vhat = hajek_est$var,
                                       domain = hajek_est$domain,
                                       na = nc_dat$na,
                                       df = nc_dat$df,
                                       adj_mat = admin1_mat)
  res <- bind_rows(hajek_est,
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
                   spatial_joint_smooth_unmatched_est) 
  res$id_sim = i
  res <- left_join(res, pop_means, by = "domain")
  return(res)
}

if (!file.exists(paste0(sim_res_dir, "spatial_pop_dat.rds"))) {
  set.seed(0412)
  ea_locs <- gen_ea_locs()
  pop_dat <- gen_pop_X(ea_locs)
  saveRDS(pop_dat, paste0(sim_res_dir, "spatial_pop_dat.rds"))
} else {
  pop_dat <- readRDS(paste0(sim_res_dir, "spatial_pop_dat.rds"))
}
set.seed(0412)
sample_ind_res <- gen_sample_ind(pop_dat, ~id_stratum, rep(4, 73))
sample_ind <- sample_ind_res$ind
pop_dat <- sample_ind_res$pop_dat

args = commandArgs(TRUE)
# supplied at the command line
k = as.numeric(args[1])
set.seed(k)
res <- do.call(
  rbind,
  lapply(1:5, function(i) run_one_sim(sample_ind,  pop_dat, gen_pop_Y, i + (k - 1) * 10))
  )
saveRDS(res, paste0(sim_res_dir, "area-level-model_res_", k, ".rds"))
