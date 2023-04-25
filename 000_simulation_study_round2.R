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
# #n.tests = 100
# #n_rand = round(runif(n = n.tests,min = 100 ,max = 2000))
# #pi.hat.naive_rand = runif(n = n.tests,min = 0.01,max = 0.99)
# #pi_rand =  runif(n = n.tests,min = 1/1000000,max = 1/10)
# #spatial_structure_draw = sample(c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia"),size = n.tests,replace = T)
# #Moran.I.corr_rand = runif(n = n.tests,min = 0.01,max = 0.99)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# round 2
n.tests = 200
n_rand = round(runif(n = n.tests,min = 100 ,max = 2000))
pi.hat.naive_rand = runif(n = n.tests,min = 0.01,max = 0.99)
pi_rand =  runif(n = n.tests,min = 1/10,max = 5/10)
spatial_structure_draw = sample(c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia"),size = n.tests,replace = T)
Moran.I.corr_rand = runif(n = n.tests,min = 0.01,max = 0.99)

new_sims = F
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
               control = list(max_treedepth = 20, adapt_delta = 0.8),
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

new_sims = F
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
new_sims =F
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

    load(file = 'generated_data/Sims/simulation_results_list_complex.RData')
    load(file = 'generated_data/Sims/simulation_results_list_simple.RData')
    load(file = 'generated_data/Sims/simulation_results_list_offset.RData')

# # # As noted is the manuscript (Simulation study, seciton C, last sentence of page 6),
# very rarely the chains will diverge for m.3 under high-levels of contamination.
# When this happens, we drop these simulations from the analysis and re-run the model.

drop <- which(
  sapply(1:length(results_list_complex$mu.mean_complex),
    function(x){
      any(abs(results_list_complex$beta.mean_complex[[x]])>1000)
      } ) )

 results_list_complex <-
    lapply(
      results_list_complex,
        function(x){
          x[-drop]
        }
      )

   results_list_simple <-
    lapply(
      results_list_simple,
        function(x){
          x[-drop]
        }
      )

   results_list_offset <-
    lapply(
      results_list_offset,
        function(x){
          x[-drop]
        }
      )
    # we have to get out-of-sample mu for each observation - this is the true value produced by the system,
    # the log-offset makes it adhere to the sampling design
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # residual  # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    residual.oob.mu_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      log_offset = log(sum(results_list_complex$data_list_complex[[x]]$y)/(results_list_complex$data_list_complex[[x]]$pi*(results_list_complex$data_list_complex[[x]]$n-sum(results_list_complex$data_list_complex[[x]]$y))) +1)
      (results_list_complex$data_list_complex[[x]]$mu-log_offset)-(results_list_complex$mu.mean_complex[[x]]-log_offset)
    })
    residual.oob.mu_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      log_offset = log(sum(results_list_simple$data_list_simple[[x]]$y)/(results_list_simple$data_list_simple[[x]]$pi*(results_list_simple$data_list_simple[[x]]$n-sum(results_list_simple$data_list_simple[[x]]$y))) +1)
      (results_list_simple$data_list_simple[[x]]$mu-log_offset)-results_list_simple$mu.mean_simple[[x]]
    })
    residual.oob.mu_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      log_offset = log(sum(results_list_offset$data_list_offset[[x]]$y)/(results_list_offset$data_list_offset[[x]]$pi*(results_list_offset$data_list_offset[[x]]$n-sum(results_list_offset$data_list_offset[[x]]$y))) +1)
      log_offset_king = log(((1-results_list_offset$data_list_offset[[x]]$pi)/results_list_offset$data_list_offset[[x]]$pi)*(mean(results_list_offset$data_list_offset[[x]]$y)/(1-mean(results_list_offset$data_list_offset[[x]]$y))))

      (results_list_offset$data_list_offset[[x]]$mu-log_offset)-(results_list_offset$mu.mean_offset[[x]]-log_offset_king)
    })

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # residual I  # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    residual.moran.I.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      W_nb =  results_list_complex$data_list_complex[[x]]$nb[results_list_complex$data_list_complex[[x]]$small_area_id,
                                                             results_list_complex$data_list_complex[[x]]$small_area_id]
      c(ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_complex[[x]])),weight = W_nb)$observed,
        ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_complex[[x]])),weight = W_nb)$p.value)
      } )

    residual.moran.I.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      W_nb =  results_list_simple$data_list_simple[[x]]$nb[results_list_simple$data_list_simple[[x]]$small_area_id,
                                                           results_list_simple$data_list_simple[[x]]$small_area_id]
      c(ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_simple[[x]])),weight = W_nb)$observed,
        ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_simple[[x]])),weight = W_nb)$p.value)
    })
    residual.moran.I.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      W_nb =  results_list_offset$data_list_offset[[x]]$nb[results_list_offset$data_list_offset[[x]]$small_area_id,
                                                           results_list_offset$data_list_offset[[x]]$small_area_id]
      c(ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_offset[[x]])),weight = W_nb)$observed,
        ape::Moran.I(x = as.numeric(unlist(residual.oob.mu_offset[[x]])),weight = W_nb)$p.value)
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # mu  # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # the idea here is that we have:
    # mu.obs = contaminated.offset + system
    # mu.pred.complex = contaminated.offset + system
    # mu.pred.simple =  system
    # mu.pred.offset = offset + system


    # bias(mu.pred.complex) = mean((mu.pred.complex-contaminated.offset)-(mu.obs - contaminated.offset))
    bias.oob.mu_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      log_offset = log(sum(results_list_complex$data_list_complex[[x]]$y)/(results_list_complex$data_list_complex[[x]]$pi*(results_list_complex$data_list_complex[[x]]$n-sum(results_list_complex$data_list_complex[[x]]$y))) +1)
      mean((results_list_complex$mu.mean_complex[[x]]-log_offset)-(results_list_complex$data_list_complex[[x]]$mu-log_offset))
      })
    # bias(mu.pred.simple) = mean((mu.pred.simple)-(mu.obs - contaminated.offset))
    bias.oob.mu_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      log_offset = log(sum(results_list_simple$data_list_simple[[x]]$y)/(results_list_simple$data_list_simple[[x]]$pi*(results_list_simple$data_list_simple[[x]]$n-sum(results_list_simple$data_list_simple[[x]]$y))) +1)
      mean((results_list_simple$mu.mean_simple[[x]])-(results_list_simple$data_list_simple[[x]]$mu-log_offset))
    })
    # bias(mu.pred.offset) = mean((mu.pred.offset - offset)-(mu.obs - contaminated.offset))
    bias.oob.mu_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      log_offset = log(sum(results_list_offset$data_list_offset[[x]]$y)/(results_list_offset$data_list_offset[[x]]$pi*(results_list_offset$data_list_offset[[x]]$n-sum(results_list_offset$data_list_offset[[x]]$y))) +1)
      log_offset_king = log(((1-results_list_offset$data_list_offset[[x]]$pi)/results_list_offset$data_list_offset[[x]]$pi)*(mean(results_list_offset$data_list_offset[[x]]$y)/(1-mean(results_list_offset$data_list_offset[[x]]$y))))
      mean((results_list_offset$mu.mean_offset[[x]]-log_offset_king)-(results_list_offset$data_list_offset[[x]]$mu-log_offset))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    rmse.oob.mu_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      log_offset = log(sum(results_list_complex$data_list_complex[[x]]$y)/(results_list_complex$data_list_complex[[x]]$pi*(results_list_complex$data_list_complex[[x]]$n-sum(results_list_complex$data_list_complex[[x]]$y))) +1)
      sqrt(mean(((results_list_complex$mu.mean_complex[[x]]-log_offset)-(results_list_complex$data_list_complex[[x]]$mu-log_offset))^2))
    })
    rmse.oob.mu_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      log_offset = log(sum(results_list_simple$data_list_simple[[x]]$y)/(results_list_simple$data_list_simple[[x]]$pi*(results_list_simple$data_list_simple[[x]]$n-sum(results_list_simple$data_list_simple[[x]]$y))) +1)
      sqrt(mean(((results_list_simple$mu.mean_simple[[x]])-(results_list_simple$data_list_simple[[x]]$mu-log_offset))^2))
    })
    rmse.oob.mu_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      log_offset = log(sum(results_list_offset$data_list_offset[[x]]$y)/(results_list_offset$data_list_offset[[x]]$pi*(results_list_offset$data_list_offset[[x]]$n-sum(results_list_offset$data_list_offset[[x]]$y))) +1)
      log_offset_king = log(((1-results_list_offset$data_list_offset[[x]]$pi)/results_list_offset$data_list_offset[[x]]$pi)*(mean(results_list_offset$data_list_offset[[x]]$y)/(1-mean(results_list_offset$data_list_offset[[x]]$y))))
      sqrt(mean(((results_list_offset$mu.mean_offset[[x]]-log_offset_king)-(results_list_offset$data_list_offset[[x]]$mu-log_offset))^2))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    cor.oob.mu_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      log_offset = log(sum(results_list_complex$data_list_complex[[x]]$y)/(results_list_complex$data_list_complex[[x]]$pi*(results_list_complex$data_list_complex[[x]]$n-sum(results_list_complex$data_list_complex[[x]]$y))) +1)
      cor(x = (results_list_complex$mu.mean_complex[[x]]-log_offset), y = (results_list_complex$data_list_complex[[x]]$mu-log_offset))
    })
    cor.oob.mu_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      log_offset = log(sum(results_list_simple$data_list_simple[[x]]$y)/(results_list_simple$data_list_simple[[x]]$pi*(results_list_simple$data_list_simple[[x]]$n-sum(results_list_simple$data_list_simple[[x]]$y))) +1)
      cor(x = results_list_simple$mu.mean_simple[[x]],y = (results_list_simple$data_list_simple[[x]]$mu-log_offset))
    })
    cor.oob.mu_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      log_offset = log(sum(results_list_offset$data_list_offset[[x]]$y)/(results_list_offset$data_list_offset[[x]]$pi*(results_list_offset$data_list_offset[[x]]$n-sum(results_list_offset$data_list_offset[[x]]$y))) +1)
      log_offset_king = log(((1-results_list_offset$data_list_offset[[x]]$pi)/results_list_offset$data_list_offset[[x]]$pi)*(mean(results_list_offset$data_list_offset[[x]]$y)/(1-mean(results_list_offset$data_list_offset[[x]]$y))))
      cor(x = (results_list_offset$mu.mean_offset[[x]]-log_offset_king),y = (results_list_offset$data_list_offset[[x]]$mu-log_offset))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # gamma # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    bias.oob.gamma_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      mean(results_list_complex$gamma.mean_complex[[x]]-(results_list_complex$data_list_complex[[x]]$gamma))
    })
    bias.oob.gamma_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      mean(results_list_simple$gamma.mean_simple[[x]]-(results_list_simple$data_list_simple[[x]]$gamma))
    })
    bias.oob.gamma_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      mean(results_list_offset$gamma.mean_offset[[x]]-(results_list_offset$data_list_offset[[x]]$gamma))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    rmse.oob.gamma_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      sqrt(mean((results_list_complex$gamma.mean_complex[[x]]-(results_list_complex$data_list_complex[[x]]$gamma))^2))
    })
    rmse.oob.gamma_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      sqrt(mean((results_list_simple$gamma.mean_simple[[x]]-(results_list_simple$data_list_simple[[x]]$gamma))^2))
    })
    rmse.oob.gamma_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      sqrt(mean((results_list_offset$gamma.mean_offset[[x]]-(results_list_offset$data_list_offset[[x]]$gamma))^2))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    cor.oob.gamma_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      cor(x = results_list_complex$gamma.mean_complex[[x]], y = (results_list_complex$data_list_complex[[x]]$gamma))
    })
    cor.oob.gamma_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      cor(x = results_list_simple$gamma.mean_simple[[x]],y = (results_list_simple$data_list_simple[[x]]$gamma))
    })
    cor.oob.gamma_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      cor(x = results_list_offset$gamma.mean_offset[[x]],y = (results_list_offset$data_list_offset[[x]]$gamma))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # beta  # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    bias.oob.beta0_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      mean(results_list_complex$beta.mean_complex[[x]][1]-(results_list_complex$data_list_complex[[x]]$beta[1]))
    })
    bias.oob.beta0_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      mean(results_list_simple$beta.mean_simple[[x]][1]-(results_list_simple$data_list_simple[[x]]$beta[1]))
    })
    bias.oob.beta0_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      mean(results_list_offset$beta.mean_offset[[x]][1]-(results_list_offset$data_list_offset[[x]]$beta[1]))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    rmse.oob.beta0_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      sqrt(mean((results_list_complex$beta.mean_complex[[x]][1]-(results_list_complex$data_list_complex[[x]]$beta[1]))^2))
    })
    rmse.oob.beta0_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      sqrt(mean((results_list_simple$beta.mean_simple[[x]][1]-(results_list_simple$data_list_simple[[x]]$beta[1]))^2))
    })
    rmse.oob.beta0_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      sqrt(mean((results_list_offset$beta.mean_offset[[x]][1]-(results_list_offset$data_list_offset[[x]]$beta[1]))^2))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    bias.oob.beta1_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      mean(results_list_complex$beta.mean_complex[[x]][2]-(results_list_complex$data_list_complex[[x]]$beta[2]))
    })
    bias.oob.beta1_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      mean(results_list_simple$beta.mean_simple[[x]][2]-(results_list_simple$data_list_simple[[x]]$beta[2]))
    })
    bias.oob.beta1_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      mean(results_list_offset$beta.mean_offset[[x]][2]-(results_list_offset$data_list_offset[[x]]$beta[2]))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    rmse.oob.beta1_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      sqrt(mean((results_list_complex$beta.mean_complex[[x]][2]-(results_list_complex$data_list_complex[[x]]$beta[2]))^2))
    })
    rmse.oob.beta1_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      sqrt(mean((results_list_simple$beta.mean_simple[[x]][2]-(results_list_simple$data_list_simple[[x]]$beta[2]))^2))
    })
    rmse.oob.beta1_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      sqrt(mean((results_list_offset$beta.mean_offset[[x]][2]-(results_list_offset$data_list_offset[[x]]$beta[2]))^2))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    bias.oob.beta2_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      mean(results_list_complex$beta.mean_complex[[x]][3]-(results_list_complex$data_list_complex[[x]]$beta[3]))
    })
    bias.oob.beta2_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      mean(results_list_simple$beta.mean_simple[[x]][3]-(results_list_simple$data_list_simple[[x]]$beta[3]))
    })
    bias.oob.beta2_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      mean(results_list_offset$beta.mean_offset[[x]][3]-(results_list_offset$data_list_offset[[x]]$beta[3]))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    rmse.oob.beta2_complex = sapply(1:length(results_list_complex$data_list_complex),function(x){
      sqrt(mean((results_list_complex$beta.mean_complex[[x]][3]-(results_list_complex$data_list_complex[[x]]$beta[3]))^2))
    })
    rmse.oob.beta2_simple = sapply(1:length(results_list_simple$data_list_simple),function(x){
      sqrt(mean((results_list_simple$beta.mean_simple[[x]][3]-(results_list_simple$data_list_simple[[x]]$beta[3]))^2))
    })
    rmse.oob.beta2_offset = sapply(1:length(results_list_offset$data_list_offset),function(x){
      sqrt(mean((results_list_offset$beta.mean_offset[[x]][3]-(results_list_offset$data_list_offset[[x]]$beta[3]))^2))
    })
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    n.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){length(results_list_complex$data_list_complex[[x]]$y)})
    n.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){length(results_list_simple$data_list_simple[[x]]$y)})
    n.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){length(results_list_offset$data_list_offset[[x]]$y)})

    pi.hat.naive.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){mean(results_list_complex$data_list_complex[[x]]$y)})
    pi.hat.naive.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){mean(results_list_simple$data_list_simple[[x]]$y)})
    pi.hat.naive.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){mean(results_list_offset$data_list_offset[[x]]$y)})

    pi.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){results_list_complex$data_list_complex[[x]]$pi})
    pi.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){results_list_simple$data_list_simple[[x]]$pi})
    pi.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){results_list_offset$data_list_offset[[x]]$pi})

    map.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){results_list_complex$data_list_complex[[x]]$spatial_map})
    map.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){results_list_simple$data_list_simple[[x]]$spatial_map})
    map.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){results_list_offset$data_list_offset[[x]]$spatial_map})

    moran.I.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){results_list_complex$data_list_complex[[x]]$moranI.mu$observed})
    moran.I.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){results_list_simple$data_list_simple[[x]]$moranI.mu$observed})
    moran.I.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){results_list_offset$data_list_offset[[x]]$moranI.mu$observed})

    moran.I.p.complex = sapply(1:length(results_list_complex$data_list_complex),function(x){results_list_complex$data_list_complex[[x]]$moranI.mu$p.value})
    moran.I.p.simple = sapply(1:length(results_list_simple$data_list_simple),function(x){results_list_simple$data_list_simple[[x]]$moranI.mu$p.value})
    moran.I.p.offset = sapply(1:length(results_list_offset$data_list_offset),function(x){results_list_offset$data_list_offset[[x]]$moranI.mu$p.value})

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

   x_list = list(
   x_n = list(complex = n.complex,simple = n.simple,offset = n.offset),
   x_pi.minus.pihat = list(complex = pi.complex-pi.hat.naive.complex,simple = pi.simple-pi.hat.naive.simple,offset = pi.offset-pi.hat.naive.offset),
   x_pi = list(complex = pi.complex,simple = pi.simple,offset = pi.offset),
   x_moranI = list(complex = moran.I.complex,simple = moran.I.simple,offset = moran.I.offset)
   )
   names(x_list) = c("n","pi-widehat(pi)","pi","Moran I")


   y_list = list(
     bias = list(
       y_bias_mu_complex = list(complex = bias.oob.mu_complex,simple = bias.oob.mu_simple,offset = bias.oob.mu_offset),
       y_bias_gamma_complex = list(complex = bias.oob.gamma_complex,simple = bias.oob.gamma_simple,offset = bias.oob.gamma_offset),
       y_bias_beta0_complex = list(complex = bias.oob.beta0_complex,simple = bias.oob.beta0_simple,offset = bias.oob.beta0_offset),
       y_bias_beta1_complex = list(complex = bias.oob.beta1_complex,simple = bias.oob.beta1_simple,offset = bias.oob.beta1_offset),
       y_bias_beta2_complex = list(complex = bias.oob.beta2_complex,simple = bias.oob.beta2_simple,offset = bias.oob.beta2_offset)
     ),
     rmse = list(
       y_rmse_mu_complex = list(complex = rmse.oob.mu_complex,simple = rmse.oob.mu_simple,offset = rmse.oob.mu_offset),
       y_rmse_gamma_complex = list(complex = rmse.oob.gamma_complex,simple = rmse.oob.gamma_simple,offset = rmse.oob.gamma_offset),
       y_rmse_beta0_complex = list(complex = rmse.oob.beta0_complex,simple = rmse.oob.beta0_simple,offset = rmse.oob.beta0_offset),
       y_rmse_beta1_complex = list(complex = rmse.oob.beta1_complex,simple = rmse.oob.beta1_simple,offset = rmse.oob.beta1_offset),
       y_rmse_beta2_complex = list(complex = rmse.oob.beta2_complex,simple = rmse.oob.beta2_simple,offset = rmse.oob.beta2_offset)
     ),
     cor = list(
       y_cor_mu_complex = list(complex = cor.oob.mu_complex,simple = cor.oob.mu_simple,offset = cor.oob.mu_offset),
       y_cor_gamma_complex = list(complex = cor.oob.gamma_complex,simple = cor.oob.gamma_simple,offset = cor.oob.gamma_offset)
     )
    )

   zlab_list = c("mu","gamma","beta[0]","beta[1]")

   plot_simulation_diagnostics = function(y,x,ylab,xlab,zlab){

     if(ylab=="cor"){
       y.range = c(min(unlist(y)),max(unlist(y)))
     }else{
     if(all(unlist(y)>0)){
       y.range = c(0,max(unlist(y))) }else{
         y.range = c(-max(abs(unlist(y))),max(abs(unlist(y)))) }
     }

     if(all(unlist(x)>0)){
       x.range = c(0,max(unlist(x))) }else{
         x.range = c(min(unlist(x)),max(unlist(x)))
       }

     col = list(complex = "purple",
                simple = "black",
                offset = "orange")
     pch = list(complex = as.integer(as.factor(map.complex)),
                simple = as.integer(as.factor(map.simple)),
                offset = as.integer(as.factor(map.offset)))
     plot(y = unlist(y),
          x = unlist(x),
          ylim = y.range,
          xlim = x.range,
          pch = NA,cex.main = 1.5,cex.lab=1.5,
          xlab = parse(text = gsub(" ","~",xlab)),
          ylab = ylab,
          main = parse(text = paste(ylab,"(widehat(",zlab,"))~by~",gsub(" ","~",xlab),sep="")),
          bty="n")
     for(i in names(y)){
       points(y = y[[i]],
              x = x[[i]],
              pch = pch[[i]],
              col = adjustcolor(col[[i]],0.25)
       )
       fit <- predict(loess( y[[i]] ~ x[[i]]),se = T)
       j <- order( x[[i]] )
       lines(x = (x[[i]])[j],y = fit$fit[j],col = col[[i]],lwd = 2)
       lines(x = (x[[i]])[j],y = fit$fit[j]+2*fit$se.fit[j],col = col[[i]],lwd = 2,lty = 2)
       lines(x = (x[[i]])[j],y = fit$fit[j]-2*fit$se.fit[j],col = col[[i]],lwd = 2,lty = 2)
     }
     if(ylab =="bias"){abline(h = 0,lty = 1,lwd = 2,col = 'darkgreen')}
   }

   pdf(file = paste('plots/Sims/mu_sims_diagnostics.pdf',sep=""),height = 10,width = 12.5)
   par(mfrow = c(length(y_list),length(x_list)))
   for(i in names(y_list)){
     for(j in names(x_list)){
   plot_simulation_diagnostics(y = y_list[[i]][[which(zlab_list =="mu")]],
                               x = x_list[[which(names(x_list)==j)]],
                               ylab = i,
                               xlab = j,
                               zlab = "mu")
     } }
   dev.off()

   pdf(file = paste('plots/Sims/intercept_sims_diagnostics.pdf',sep=""),height = 5,width = 12.5)
   par(mfrow = c(2,length(x_list)))
   for(i in names(y_list)[-which(names(y_list)=="cor")]){
     for(j in names(x_list)){
       plot_simulation_diagnostics(y = y_list[[i]][[which(zlab_list =="beta[0]")]],
                                   x = x_list[[which(names(x_list)==j)]],
                                   ylab = i,
                                   xlab = j,
                                   zlab = "beta[1]")
     } }
   dev.off()

   pdf(file = paste('plots/Sims/covariate_sims_diagnostics.pdf',sep=""),height = 5,width = 12.5)
   par(mfrow = c(2,length(x_list)))
   for(i in names(y_list)[-which(names(y_list)=="cor")]){
     for(j in names(x_list)){
       plot_simulation_diagnostics(y = y_list[[i]][[which(zlab_list =="beta[1]")]],
                                   x = x_list[[which(names(x_list)==j)]],
                                   ylab = i,
                                   xlab = j,
                                   zlab = "beta[2]")
     } }
   dev.off()

   pdf(file = paste('plots/Sims/gamma_sims_diagnostics.pdf',sep=""),height = 10,width = 12.5)
   par(mfrow = c(length(y_list),length(x_list)))
   for(i in names(y_list)){
     for(j in names(x_list)){
       plot_simulation_diagnostics(y = y_list[[i]][[which(zlab_list =="gamma")]],
                                   x = x_list[[which(names(x_list)==j)]],
                                   ylab = i,
                                   xlab = j,
                                   zlab = "gamma")
     } }
   dev.off()



   pdf(file = paste('plots/Sims/mu_sims_diagnostics_rec3_v_KingZeng.pdf',sep=""),height = 10,width = 12.5)
   par(mfrow = c(length(y_list),length(x_list)))
   for(i in names(y_list)){
     for(j in names(x_list)){
       plot_simulation_diagnostics(y = y_list[[i]][[1]][-which(names(y_list[[i]][[1]])=="simple")],
                                   x = x_list[[which(names(x_list)==j)]][-which(names(y_list[[i]][[1]])=="simple")],
                                   ylab = i,
                                   xlab = j,
                                   zlab = "mu")
     } }
   dev.off()

   pdf(file = paste('plots/Sims/intercept_sims_diagnostics_rec3_v_KingZeng.pdf',sep=""),height = 5,width = 12.5)
   par(mfrow = c(2,length(x_list)))
   for(i in names(y_list)[-which(names(y_list)=="cor")]){
     for(j in names(x_list)){
       plot_simulation_diagnostics(y = y_list[[i]][[which(zlab_list =="beta[0]")]][-which(names(y_list[[i]][[which(zlab_list =="beta[0]")]])=="simple")],
                                   x = x_list[[which(names(x_list)==j)]][-which(names(y_list[[i]][[which(zlab_list =="beta[0]")]])=="simple")],
                                   ylab = i,
                                   xlab = j,
                                   zlab = "beta[1]")
     } }
   dev.off()