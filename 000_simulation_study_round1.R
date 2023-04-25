gc()
options(scipen=999)
#libraries
libs <- c("data.table", "dplyr")
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib)
    library(lib, character.only = TRUE)
  } else {
    library(lib, character.only = TRUE)
  }
}

# Load utility functions
source("utils/utils.R")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # DEFINE TESTING GROUND # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# round 1
n.tests = 100
n_rand = round(runif(n = n.tests,min = 100 ,max = 2000))
pi.hat.naive_rand = runif(n = n.tests,min = 0.01,max = 0.99)
pi_rand =  runif(n = n.tests,min = 1/1000000,max = 1/10)
spatial_structure_draw = sample(c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia"),size = n.tests,replace = T)
Moran.I.corr_rand = runif(n = n.tests,min = 0.01,max = 0.99)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# round 2
#n.tests = 200
#n_rand = round(runif(n = n.tests,min = 100 ,max = 2000))
#pi.hat.naive_rand = runif(n = n.tests,min = 0.01,max = 0.99)
#pi_rand =  runif(n = n.tests,min = 1/10,max = 5/10)
#spatial_structure_draw = sample(c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia"),size = n.tests,replace = T)
#Moran.I.corr_rand = runif(n = n.tests,min = 0.01,max = 0.99)

new_sims = T
begin_iter = 0
optimal.iter = 1000
    # Complex Model
if(new_sims == T){
  mu.mean_complex = list()
  beta.mean_complex = list()
  gamma.mean_complex = list()
  y_gen.mean_complex = list()
  rhat_complex = list()
  ess_bulk_complex = list()
  data_list_complex = list()
}else{
  load(file = 'generated_data/Sims/simulation_results_list_complex.RData')
  begin_iter = length(results_list_complex$mu.mean_complex)
  mu.mean_complex = results_list_complex$mu.mean_complex
  beta.mean_complex = results_list_complex$beta.mean_complex
  gamma.mean_complex = results_list_complex$gamma.mean_complex
  y_gen.mean_complex = results_list_complex$y_gen.mean_complex
  rhat_complex = results_list_complex$rhat_complex
  ess_bulk_complex = results_list_complex$ess_bulk_complex
  data_list_complex = results_list_complex$data_list_complex
}
for(iter in (begin_iter+1):n.tests){

  sims = dgf(n.sims = 1,
                  n = n_rand[iter],
                  pi.hat.naive = pi.hat.naive_rand[iter],
                  p = 1,
                  X_corr = 0,
                  pi = pi_rand[iter],
                  Moran.I.corr = Moran.I.corr_rand[iter],
                  spatial_structure = spatial_structure_draw[iter] )
  data_sim = sims[-which(names(sims)=="sp.object")]
  data_sim$y = as.numeric(data_sim$y)
  data_sim$spatial_map = spatial_structure_draw[iter]

    # fit a complex model with `computationally cheap` settings
fit_object = fit(data = data_sim,
               show_code = T,
               contamination = T,
               offset = T,
               beta_prior = "cauchy",
               small_area_prior = "BYM2",
               intercept_scale = 10,
               large_area_prior = NA,
               # stan functionalities
               iter = optimal.iter*4,
               warmup = optimal.iter*4*(2/3),
               thin = 4,
               cores = 4,
               chains = 4,
               control = list(max_treedepth = 10, adapt_delta = 0.8),
               verbose = TRUE)

    # extract mu and check stability over iters
    summary_fit_object = summary(fit_object)$summary
    mu.mean_complex = append(mu.mean_complex,list(summary_fit_object[grep("mu",rownames(summary_fit_object )),][,"mean"]))
    beta.mean_complex = append(beta.mean_complex,list(summary_fit_object[grep("beta",rownames(summary_fit_object )),][,"mean"]))
    gamma.mean_complex = append(gamma.mean_complex,list(summary_fit_object[grep("gamma\\[",rownames(summary_fit_object )),][,"mean"]))
    y_gen.mean_complex = append(y_gen.mean_complex,list(summary_fit_object[grep("y_gen",rownames(summary_fit_object )),][,"mean"]))
    rhat_complex = append(rhat_complex,list(summary_fit_object[,"Rhat"]))
    ess_bulk_complex = append(ess_bulk_complex,list(summary_fit_object[,"n_eff"]))
    data_list_complex = append(data_list_complex,list(data_sim))

    results_list_complex = list( mu.mean_complex = mu.mean_complex,
                                 beta.mean_complex = beta.mean_complex,
                                 gamma.mean_complex = gamma.mean_complex,
                                 y_gen.mean_complex = y_gen.mean_complex,
                                 rhat_complex = rhat_complex,
                                 ess_bulk_complex = ess_bulk_complex,
                                 data_list_complex = data_list_complex)
    save(results_list_complex ,file = 'generated_data/Sims/simulation_results_list_complex.RData')
}

new_sims = T
begin_iter = 0
optimal.iter = 1000
    # Simple Model
if(new_sims == T){
    mu.mean_simple = list()
    beta.mean_simple = list()
    gamma.mean_simple = list()
    y_gen.mean_simple = list()
    rhat_simple = list()
    ess_bulk_simple = list()
    data_list_simple = list()
  }else{
    load(file = 'generated_data/Sims/simulation_results_list_simple.RData') ;
    begin_iter = length(results_list_simple$mu.mean_simple)
    mu.mean_simple = results_list_simple$mu.mean_simple
    beta.mean_simple = results_list_simple$beta.mean_simple
    gamma.mean_simple = results_list_simple$gamma.mean_simple
    y_gen.mean_simple = results_list_simple$y_gen.mean_simple
    rhat_simple = results_list_simple$rhat_simple
    ess_bulk_simple = results_list_simple$rhat_simple
    data_list_simple = results_list_simple$data_list_simple
  }
for(iter in (begin_iter+1):n.tests){

      sims = dgf(n.sims = 1,
                      n = n_rand[iter],
                      pi.hat.naive = pi.hat.naive_rand[iter],
                      p = 1,
                      X_corr = 0,
                      pi = pi_rand[iter],
                      Moran.I.corr = Moran.I.corr_rand[iter],
                      spatial_structure = spatial_structure_draw[iter] )
      data_sim = sims[-which(names(sims)=="sp.object")]
      data_sim$y = as.numeric(data_sim$y)
      data_sim$spatial_map = spatial_structure_draw[iter]

    # fit a complex model with `computationally cheap` settings
fit_object = fit(data = data_sim,
               show_code = T,
               contamination = F,
               offset = F,
               beta_prior = "normal",
               small_area_prior = "fixed",
               intercept_scale = 10,
               large_area_prior = NA,
               # stan functionalities
               iter = optimal.iter*4,
               warmup = optimal.iter*4*(2/3),
               thin = 4,
               cores = 4,
               chains = 4,
               control = list(max_treedepth = 10, adapt_delta = 0.8),
               verbose = TRUE)
    # extract mu and check stability over iters
    summary_fit_object = summary(fit_object)$summary
    mu.mean_simple = append(mu.mean_simple,list(summary_fit_object[grep("mu",rownames(summary_fit_object )),][,"mean"]))
    beta.mean_simple = append(beta.mean_simple,list(summary_fit_object[grep("beta",rownames(summary_fit_object )),][,"mean"]))
    gamma.mean_simple = append(gamma.mean_simple,list(summary_fit_object[grep("gamma\\[",rownames(summary_fit_object )),][,"mean"]))
    y_gen.mean_simple = append(y_gen.mean_simple,list(summary_fit_object[grep("y_gen",rownames(summary_fit_object )),][,"mean"]))
    rhat_simple = append(rhat_simple,list(summary_fit_object[,"Rhat"]))
    ess_bulk_simple = append(ess_bulk_simple,list(summary_fit_object[,"n_eff"]))
    data_list_simple = append(data_list_simple,list(data_sim))

    results_list_simple = list( mu.mean_simple = mu.mean_simple,
                                 beta.mean_simple = beta.mean_simple,
                                 gamma.mean_simple = gamma.mean_simple,
                                 y_gen.mean_simple = y_gen.mean_simple,
                                 rhat_simple = rhat_simple,
                                 ess_bulk_simple = ess_bulk_simple,
                                 data_list_simple = data_list_simple)
    save(results_list_simple ,file = 'generated_data/Sims/simulation_results_list_simple.RData')
}

new_sims = T
begin_iter = 0
optimal.iter = 1000
    # Offset-only Model
if(new_sims == T){
    mu.mean_offset = list()
    beta.mean_offset = list()
    gamma.mean_offset = list()
    y_gen.mean_offset = list()
    rhat_offset = list()
    ess_bulk_offset = list()
    data_list_offset = list()
    }else{
      load(file = 'generated_data/Sims/simulation_results_list_offset.RData') ;
      begin_iter = length(results_list_offset$mu.mean_offset)
      mu.mean_offset = results_list_offset$mu.mean_offset
      beta.mean_offset = results_list_offset$beta.mean_offset
      gamma.mean_offset = results_list_offset$gamma.mean_offset
      y_gen.mean_offset = results_list_offset$y_gen.mean_offset
      rhat_offset = results_list_offset$rhat_offset
      ess_bulk_offset = results_list_offset$ess_bulk_offset
      data_list_offset = results_list_offset$data_list_offset
      }
for(iter in (begin_iter + 1):n.tests){

      sims = dgf(n.sims = 1,
                      n = n_rand[iter],
                      pi.hat.naive = pi.hat.naive_rand[iter],
                      p = 1,
                      X_corr = 0,
                      pi = pi_rand[iter],
                      Moran.I.corr = Moran.I.corr_rand[iter],
                      spatial_structure = spatial_structure_draw[iter] )
      data_sim = sims[-which(names(sims)=="sp.object")]
      data_sim$y = as.numeric(data_sim$y)
      data_sim$spatial_map = spatial_structure_draw[iter]

      # fit a complex model with `computationally cheap` settings
fit_object = fit(data = data_sim,
                 show_code = T,
                 contamination = F,
                 offset = T,
                 beta_prior = "normal",
                 small_area_prior = "fixed",
                 intercept_scale = 10,
                 large_area_prior = NA,
                 # stan functionalities
                 iter = optimal.iter*4,
                 warmup = optimal.iter*4*(2/3),
                 thin = 4,
                 cores = 4,
                 chains = 4,
                 control = list(max_treedepth = 10, adapt_delta = 0.8),
                 verbose = TRUE)
      # extract mu and check stability over iters
      summary_fit_object = summary(fit_object)$summary
      mu.mean_offset = append(mu.mean_offset,list(summary_fit_object[grep("mu",rownames(summary_fit_object )),][,"mean"]))
      beta.mean_offset = append(beta.mean_offset,list(summary_fit_object[grep("beta",rownames(summary_fit_object )),][,"mean"]))
      gamma.mean_offset = append(gamma.mean_offset,list(summary_fit_object[grep("gamma\\[",rownames(summary_fit_object )),][,"mean"]))
      y_gen.mean_offset = append(y_gen.mean_offset,list(summary_fit_object[grep("y_gen",rownames(summary_fit_object )),][,"mean"]))
      rhat_offset = append(rhat_offset,list(summary_fit_object[,"Rhat"]))
      ess_bulk_offset = append(ess_bulk_offset,list(summary_fit_object[,"n_eff"]))
      data_list_offset = append(data_list_offset,list(data_sim))

      results_list_offset = list( mu.mean_offset = mu.mean_offset,
                                  beta.mean_offset = beta.mean_offset,
                                  gamma.mean_offset = gamma.mean_offset,
                                  y_gen.mean_offset = y_gen.mean_offset,
                                  rhat_offset = rhat_offset,
                                  ess_bulk_offset = ess_bulk_offset,
                                  data_list_offset = data_list_offset)
      save(results_list_offset ,file = 'generated_data/Sims/simulation_results_list_offset.RData')
}
