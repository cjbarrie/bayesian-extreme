gc()
options(scipen=999)

#libraries
libs <- c("data.table", "dplyr", "xtable")
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

# load individual-level training data 
individual.data_path <- 'data/matchdata_isistunsc.csv'
survey = fread(file = individual.data_path)

# load shape file 
shapefile_path <- "shapefiles/tunisia_shapefiles/TUN_adm2.shp"
shape = sf::st_read(shapefile_path)

# get denominator for prevalence adjustment 
denom = fread(file = 'data/mena_pops.csv')

# some changes to names to ensure clean matching
survey$gov_name_f = as.factor(survey$gov_name_f) 
levels(survey$gov_name_f)[levels(survey$gov_name_f)=="B√©ja"] = "Beja"
levels(survey$gov_name_f)[levels(survey$gov_name_f)=="Gab√®s"] = "Gabes"
levels(survey$gov_name_f)[levels(survey$gov_name_f)=="M√©denine"] = "Mednine"
levels(survey$gov_name_f)[levels(survey$gov_name_f)=="Zaghouane"] = "Zaghouan"
levels(survey$gov_name_f)[levels(survey$gov_name_f)=="Manouba"] = "Manubah"

shape$gov_name_f = as.factor(shape$gov_name_f)
levels(shape$gov_name_f) [which(levels(shape$gov_name_f)=="Béja")] = "Beja"
levels(shape$gov_name_f) [which(levels(shape$gov_name_f)=="Gabès")] = "Gabes"
levels(shape$gov_name_f) [which(levels(shape$gov_name_f)=="Médenine")] = "Mednine"
levels(shape$gov_name_f) [which(levels(shape$gov_name_f)=="Kassérine")] = "Kasserine"
levels(shape$gov_name_f) [which(levels(shape$gov_name_f)=="Sidi Bou Zid")] = "Sidi Bouzid"

# prep data 
data.list = 
  data.prep(shape = shape, survey = survey[,!c("deleg_na_1","gov_name_f")],
            shape_large.area_id_name = "gov_name_f", shape_large.area_id_num = NA,
            shape_small.area_id_name = "deleg_na_1",shape_small.area_id_num = "adm_id",
            survey_small.area_id_name = NA, survey_small.area_id_num = "adm_id",
            drop.incomplete.records = T,
            colnames_X = c("coledu","age","married","student","lowstat",
                           "pop_10plus_2014","population_density","pct_agri_2014",
                           "pct_higher_edu_2014","unemp_rate_2014","dip_unemp_rate_2014",
                           "sqrt_distance_to_libya","sqrt_post_rev_protest_events","pct_nahdha_2014",
                           "pct_nahdha_2011"),
            interactions_list = list(age2 = "age*age",coledu_lowstat = "coledu*lowstat"),
            scale_X = "1sd",
            colname_y = "case",
            pi = 7000/denom$male_18_sunni[which(denom$country=="Tunisia")],
            # source: What Explains the Flow of Foreign Fighters to ISIS? Efraim Benmelech & Esteban F. Klor, 2020
            large_area_shape = T
  )

# check that the nb object is fully connected 
testconnected(data.list$nb_object)

# plot edges and nodes over map
pdf(file = 'Plots/Tunisia/connected.pdf',width = 5, height = 5)
plot_network(shape = shape, nb_object = data.list$nb_object)
dev.off()

# fit the model 
begin.time = Sys.time()
fit_object = fit(data = data.list,
                 show_code = T,
                 contamination = T,offset = T,beta_prior = "cauchy",small_area_prior = "BYM2",intercept_scale = 10,large_area_prior = "random",
                 iter = 25000,warmup = 22500,thin = 4,cores = 4,chains = 4,control = list(max_treedepth = 25, adapt_delta = 0.99),verbose = T)
time.taken = Sys.time() - begin.time

# save the fit object
save(fit_object,file = 'generated_data/Tunisia/fit_object.RData',compress = TRUE)
load(file = 'generated_data/Tunisia/fit_object.RData')
# load sources
source('utils/utils_plot.R')

# # # # # # # # # # # # # # # 
# # # CONVERGENCE PLOTS # # #
# # # # # # # # # # # # # # # 

# # # get convergence monitoring quantities for each parameter
mon <- monitor(fit_object)

pdf(file = 'plots/Tunisia/convergence_min.ESS_bulk.and.tail.pdf',width =10,height =7.5)
plotConvergence_minESS(mon = mon)
dev.off()

pdf(file = 'plots/Tunisia/convergence_mixing_min.ESS_bulk.and.tail.pdf',width =15,height =10)
plotConvergence_mixing_rank.hist_minESS(mon = mon,breaks = 20)
dev.off()

pdf(file = 'plots/Tunisia/convergence_rhat.pdf',width =10,height =10)
plotConvergence_Rhat(fit_object = fit_object)
dev.off()

pdf(file = 'plots/Tunisia/convergence_ESS.pdf',width =10,height =12.5)
plotConvergence_ESS(fit_object = fit_object)
dev.off()

pdf(file = 'plots/Tunisia/ppc.pdf',width =7.5,height =5)
plot_ppc(fit_object = fit_object,data.list = data.list)
dev.off()

# traceplots
pdf(file = paste('plots/Tunisia/convergence_traceplot_beta.pdf',sep=""),height = 5,width = 15)
traceplot(object = fit_object,pars = 'beta')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_phi.pdf',sep=""),height = 20,width = 30)
traceplot(object = fit_object,pars = 'phi')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_psi.pdf',sep=""),height = 20,width = 30)
traceplot(object = fit_object,pars = 'psi')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_gamma.pdf',sep=""),height = 20,width = 30)
traceplot(object = fit_object,pars = 'gamma')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_eta.pdf',sep=""),height = 10,width = 15)
traceplot(object = fit_object,pars = 'eta')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_lambda.pdf',sep=""),height = 5,width = 15)
traceplot(object = fit_object,pars = 'lambda')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_sigma_gamma.pdf',sep=""),height = 5,width = 15)
traceplot(object = fit_object,pars = 'sigma_gamma')
dev.off()
pdf(file = paste('plots/Tunisia/convergence_traceplot_sigma_eta.pdf',sep=""),height = 5,width = 15)
traceplot(object = fit_object,pars = 'sigma_eta')
dev.off()

# # # # # # # # # # # # # # # 
# # # PARAMETER PLOTS # # # # 
# # # # # # # # # # # # # # # 

# extract MCMC samples for each parameter if interest
extracted = rstan::extract(fit_object , 
                    pars = c('beta','phi','psi','gamma','sigma_gamma','lambda','eta','sigma_eta','mu','y_gen'), 
                    permuted = TRUE, 
                    inc_warmup = FALSE,
                    include = TRUE)

# define names of covariates according to the design matrix
pred.names = colnames(data.list$X)
# identify individual-level covariates amongst these 
plot.names.ind = c("intercept",
                   "married",
                   "student",
                   "coledu",
                   "lowstat",
                   "coledu_lowstat",
                   "age",
                   "age2")
# define new plot-friendly names for each variable
plot.names = c(
  "Intercept",
  "College edu.",
  "Age",
  "Married",
  "Student",
  "Low status",
  "Population",
  "Population density",
  "% Agriculture",
  "% College edu.",
  "Unemp. rate",
  "Graduate Unemp. rate",
  "Sqrt. Distance to Libya",
  'Sqrt. N. Post-Revoltion Protests',
  "% Votes Nahdha 2014",
  "% Votes Nahdha 2011",
  "Age^2",
  "College edu.*Low status"
  )

# identify individual-level variables amongst these
plot.names.new.ind = c("Intercept",
                       "Married",
                       "Student",
                       "College edu.",
                       "Low status",
                       "College edu.*Low status",
                       "Age",
                       "Age^2")

# # # calculate relative deprivation marginal effects 
# logit scale
pdf(file = paste('plots/Tunisia/marginal.effect_relative.deprivation_logit.pdf',sep=""),height = 5,width = 7.5)
plotPredicted_2way(scale = "logit",
                   n = denom$male_18_sunni[which(denom$country=="Tunisia")],
                   interaction.vars.table = data.table(expand.grid(coledu = c(0, 1),lowstat = c(0, 1))),
                   data.list = data.list,
                   pred.names = pred.names,
                   extracted_beta = extracted$beta,
                   interaction.name = "coledu_lowstat",
                   outcome.lab = "recruitment",
                   interaction.states.cols = c("skyblue", "lightgreen", "lightcoral", "grey"),
                   interaction.states.names =  c("no college education & high status",
                                                 "college educated & high status",
                                                 "no college education & low status",
                                                 "college educated & low status"),
                   y.max = 1.35)
dev.off()
# n
pdf(file = paste('plots/Tunisia/marginal.effect_relative.deprivation_n.pdf',sep=""),height = 7.5,width = 15)
plotPredicted_2way(scale = "n",
                   n = denom$male_18_sunni[which(denom$country=="Tunisia")],
                   interaction.vars.table = data.table(expand.grid(coledu = c(0, 1),lowstat = c(0, 1))),
                   data.list = data.list,
                   pred.names = pred.names,
                   extracted_beta = extracted$beta,
                   interaction.name = "coledu_lowstat",
                   outcome.lab = "recruits",
                   interaction.states.cols = c("skyblue", "lightgreen", "lightcoral", "grey"),
                   interaction.states.cols.dark = c("blue", "darkgreen", "red", "black"),
                   interaction.states.names =  c("no college education & high status",
                                                 "college educated & high status",
                                                 "no college education & low status",
                                                 "college educated & low status"),
                   separate.plots = T)
dev.off()
# relative probability
pdf(file = paste('plots/Tunisia/marginal.effect_relative.deprivation_relprob.pdf',sep=""),height = 7.5,width = 12.5)
plotPredicted_2way(scale = "rel",
                   n = denom$male_18_sunni[which(denom$country=="Tunisia")],
                   interaction.vars.table = data.table(expand.grid(coledu = c(0, 1),lowstat = c(0, 1))),
                   data.list = data.list,
                   pred.names = pred.names,
                   extracted_beta = extracted$beta,
                   interaction.name = "coledu_lowstat",
                   outcome.lab = "recruits",
                   interaction.states.cols = c("skyblue", "lightgreen", "lightcoral", "grey"),
                   interaction.states.cols.dark = c("blue", "darkgreen", "red", "black"),
                   interaction.states.names =  c("no college education & high status",
                                                 "college educated & high status",
                                                 "no college education & low status",
                                                 "college educated & low status"),
                   separate.plots = T)
dev.off()

# # # regression coefficients
# individual level
pdf(file = 'plots/Tunisia/regression.coefficients_ind_standardized.pdf',height = 6,width = 7.5)
plotbetas(extractedbetas = extracted$beta,
          mfrow = c(3,3),
          beta.names = pred.names,
          beta.names.plot = plot.names.ind,
          beta.names.plot.new = plot.names.new.ind ,
          include_legend = T)
dev.off()
pdf(file = 'plots/Tunisia/regression.coefficients_ind_original.scale.pdf',height = 6,width = 7.5)
plotbetas(extractedbetas = extracted$beta,
          mfrow = c(3,3),
          beta.names = pred.names,
          beta.names.plot = plot.names.ind,
          beta.names.plot.new = plot.names.new.ind ,
          include_legend = T,
          descale = data.list$scale_X,
          mu.X = data.list$mu_X,
          sd.X = data.list$sd_X)
dev.off()
# contextual level
pdf(file = 'plots/Tunisia/regression.coefficients_context_standardized.pdf',height = 5,width = 15)
plotbetas(extractedbetas = extracted$beta,
          mfrow = c(2,5),
          beta.names = pred.names,
          beta.names.plot = pred.names[-which(pred.names %in% plot.names.ind)],
          beta.names.plot.new = plot.names[-which(plot.names %in% plot.names.new.ind)],
          include_legend = T)
dev.off()
pdf(file = 'plots/Tunisia/regression.coefficients_context_original.scale.pdf',height = 5,width = 15)
plotbetas(extractedbetas = extracted$beta,
          mfrow = c(2,5),
          beta.names = pred.names,
          beta.names.plot = pred.names[-which(pred.names %in% plot.names.ind)],
          beta.names.plot.new = plot.names[-which(plot.names %in% plot.names.new.ind)],
          include_legend = T,
          descale = data.list$scale_X,
          mu.X = data.list$mu_X,
          sd.X = data.list$sd_X)
dev.off()
# lambda density
pdf(file = 'plots/Tunisia/lambda.pdf',height = 5,width = 5)
plotlambdas(extractedlambdas = extracted$lambda)
dev.off()

# plot Moran I 
pdf(file = 'plots/Tunisia/moranI.pdf',height = 5,width = 5)
plotMoran(data.list = data.list,extractedy = extracted$y)
dev.off()

# plot BYM2 small-area effects
pdf(file = 'plots/Tunisia/totalBYM2_smallarea.pdf',height = 7.5,width = 20)
plotTotalBYM2(data.list = data.list,extracted = extracted, 
              oma = c(5, 0, 3, 0), mar = c(5, 5, 7.5, 2))
dev.off()

# plot large-area effects
pdf(file = 'plots/Tunisia/largearea.pdf',height = 7.5,width = 10)
plotLargeArea(data.list = data.list,
              extracted = extracted, 
              oma = c(5, 0, 3, 0), mar = c(5, 5, 7.5, 2))
dev.off()

# plot outcome on a map
pdf(file = 'plots/Tunisia/map_outcome_smallarea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "y",shape_area_id = "adm_id",area = "small")
dev.off()

# plot residual on a map
pdf(file = 'plots/Tunisia/map_residual_smallarea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "residual",shape_area_id = "adm_id",area = "small")
dev.off()

# plot BYM2 small-area effect on a map
pdf(file = 'plots/Tunisia/map_BYM2_smallarea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "gamma",shape_area_id = "adm_id",area = "small")
dev.off()

# plot spatial small-area effect on a map
pdf(file = 'plots/Tunisia/map_ICAR_smallarea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "psi",shape_area_id =  "adm_id",area = "small")
dev.off()

# plot unstructured small-area effect on a map
pdf(file = 'plots/Tunisia/map_unstructured_smallarea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "phi",shape_area_id = "adm_id",area = "small")
dev.off()

# plot unstructured large-area effect on a map
pdf(file = 'plots/Tunisia/map_unstructured_largearea.pdf',height = 5,width = 7.5)
map.param(data.list = data.list,extracted = extracted,param = "eta",shape_area_id = "gov_name_f",area = "large")
dev.off()

# generate distribution of recruitment by all possible profiles - holding location constant 
X.pred = 
  expand.grid(married = c(0,1),
              student = c(0,1),
              coledu = c(0,1),
              lowstat = c(0,1),
              age = round(seq(18,max(data.list$survey_complete$age),length.out =10)))
X.pred.original = X.pred
X.pred $coledu_lowstat = X.pred$coledu*X.pred $lowstat
X.pred$age2 = X.pred$age^2
# prepare labels
X.pred.labels = X.pred
X.pred.labels$married = ifelse(X.pred.labels$married ==1,"married","not married")
X.pred.labels$student = ifelse(X.pred.labels$student==1,"studying","not studying")
X.pred.labels$coledu = ifelse(X.pred.labels$coledu==1,"college educated","not college educated")
X.pred.labels$lowstat = ifelse(X.pred.labels$lowstat==1,"low status","high status")
X.pred.labels$age = paste(X.pred.labels$age,"years of age")
X.pred.labels = apply(X.pred.labels[,c("coledu","lowstat","age","married","student")] ,1,paste0,collapse = "\n")
# standardize appropriately
if(data.list$scale_X=="1sd"){
  X.pred = (X.pred-t(array(data.list$mu_X[names(X.pred)],c(dim(X.pred)[2],dim(X.pred)[1]))))/t(array(data.list$sd_X[names(X.pred)],c(dim(X.pred)[2],dim(X.pred)[1])))
}
if(data.list$scale_X=="2sd"){
  X.pred = (X.pred-t(array(data.list$mu_X[names(X.pred)],c(dim(X.pred)[2],dim(X.pred)[1]))))/t(array(2*data.list$sd_X[names(X.pred)],c(dim(X.pred)[2],dim(X.pred)[1])))
}
# add intercept
X.pred$intercept = 1
# make logit-level predictions
logit.pred.effect =  sapply(1:dim(extracted$beta)[1],function(s){
  as.matrix(X.pred[,plot.names.ind]) %*% as.numeric(extracted$beta[s,match(plot.names.ind,pred.names)]) 
})
logit.pred.effect.quantiles  = apply(logit.pred.effect ,1,function(x){quantile(x,c(0.05,0.5,0.95))})
# make probability-level predictions
prob.pred.effect = inv_logit(logit.pred.effect )
prob.pred.effect.quantiles = inv_logit(logit.pred.effect.quantiles )
prob.pred.effect.baseline = t(array(inv_logit(extracted$beta[,1]),c(dim(prob.pred.effect)[2],dim(prob.pred.effect)[1])))
# make odds-level predictions
odds = prob.pred.effect /prob.pred.effect.baseline
odds.quantiles = apply(odds ,1,function(x){quantile(x,c(0.05,0.5,0.95))})

# plot ranked logit profiles
pdf(file = 'plots/Tunisia/profile_ranks_logit.pdf',height = 6,width = 10)
profile_ranks_logit (sims = logit.pred.effect,
                     sims.quantiles = logit.pred.effect.quantiles,
                     xlim = c(-5,165),
                     ylim = c(-60,15),
                     outcome  = 'recruitment',
                     labels = X.pred.labels)
dev.off() 

to_frac = function(x, den) {paste0(round(x * den), "/", den)}

temp = data.table(apply(X.pred.original,2,as.integer),
                  median_prob = prob.pred.effect.quantiles[2,],
                  median_prob_frac = to_frac(prob.pred.effect.quantiles[2,],10000),
                  median_odds = odds.quantiles[2,],
                  median_log_odds = logit.pred.effect.quantiles[2,]
)
temp = temp[rev(order(temp$median_prob))[1:10],]
tab.tex <- xtable(temp , include.rownames=FALSE,digits = c(0,rep(6,dim(temp)[2]-2),3,3))
print(tab.tex, file = "tables/profiles_tun.tex", compress = FALSE)

# plot ranked odds profiles
pdf(file = 'plots/Tunisia/profile_ranks_odds_v_average.pdf',height = 10,width = 15)
plot.ranks.quantiles (sims = odds,
                      sims.quantiles = odds.quantiles,
                      labels = X.pred.labels,
                      outcome = 'recruitment',
                      unit = 'odds',
                      ylim.point = 500,
                      ylim.min = 175,
                      ylim.max = 2000,
                      xlim.min = -40,
                      xlim.max = 165)
dev.off()

# plot ranked probability profiles
pdf(file = 'plots/Tunisia/profile_ranks_absolute_probability.pdf',height = 10,width = 17.5)
plot.ranks.quantiles (sims = prob.pred.effect,
                      sims.quantiles = prob.pred.effect.quantiles,
                      labels = X.pred.labels,
                      outcome = 'recruitment',
                      unit = 'probability',
                      ylim.point = 10/100,
                      ylim.min = 30/1000,
                      ylim.max = 3/10,
                      xlim.min = -40,
                      xlim.max = 164)
dev.off()





