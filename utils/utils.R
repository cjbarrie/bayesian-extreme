# NOTE: You need to install Stan, JAGS and INLA to run these functions properly

library(INLA) # necessary for calculation of scaling factor 

library(MASS) # necessary to simulate from multivariate normal 

# necessary for handling shapefiles and various spatial data tools
library(sf) 
library(spdep)

library(geodist) # necessary for calculation of distance between coordiates 

library(ape) # contains Moran I calculator
# Calculate Moran's  I with bootleg function, over the simulated residuals
# (row-normalized weight matrix) John L. Gittleman, Mark Kot,
# Adaptation: Statistics and a Null Model for Estimating Phylogenetic Effects
Moran.I.bootleg = function(x,W){ 
  N <- dim(W)[1];
  x.bar <- mean(x);
  dx <- x - x.bar;
  g <- expand.grid(dx, dx);
  xixj <- g[,1] * g[,2];
  num <- N*sum(W*matrix(xixj, ncol=N));
  denom <- sum(W)*sum(dx^2);
  I = num/denom;
  return(I)
}

library(SpatialEpi) # to get maps 
library(R2jags) # requires JAGS to simulate conditionally autoregressive data 

library(parallel)
library(rstan) # need stan to fit model 
# Rstan Options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# utils
inv_logit = function(x){exp(x)/(1+exp(x))}

# sample function is useful but buggy - if you specify a single integer it returns a sequence up to that integer
sample = function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}

# fully connect the graph
addnbs <- function(sp.sample,ID,D){
  queen_nb <- poly2nb(sp.sample, row.names=ID, queen=TRUE)
  count = card(queen_nb)
  if(!any(count==0)){
    return(queen_nb)
  }
  ## get nearest neighbour:
  for(i in which(count==0)){
    queen_nb[[i]] =  order(D[i,])[2]
    queen_nb[[order(D[i,])[2]]] = sort(c(i,queen_nb[[order(D[i,])[2]]]))
  }
  return(queen_nb)
}

# add specific neighbours - useful to merge shapefiles and get correct neighbourhood
add_specific_nbs <- function(nb,name1,name2,IDs){
  nb[[which(IDs==name1)]] = sort(c(nb[[which(IDs==name1)]],which(IDs==name2)))
  nb[[which(IDs==name2)]] = sort(c(nb[[which(IDs==name2)]],which(IDs==name1)))
  return(nb)
}

# Collect and connect shapefiles
collect_and_connect_shapes <- function(shape_paths_list,
                                       adm_colnames_list,
                                       plotnames_list,
                                       plot_path){
  shape_total = st_sf(ADM1_EN = character(0),geometry = st_sfc(),ADM0_EN=character(0))
  for(i in 1:length(paths_list)){
    shape_temp = sf::st_read(dsn = shape_paths_list[[i]])
    shape_temp$ADM1_EN <- as.factor(shape_temp[[adm_colnames_list[[i]]]])
    shape_temp_all = shape_temp[,c("ADM1_EN","geometry")]
    shape_temp_all = 
      shape_temp_all %>% 
      group_by(ADM1_EN) %>%
      summarise(geometry = sf::st_union(geometry)) %>%
      ungroup()
    shape_temp_all$ADM0_EN = str_to_title(gsub("_shape","",plotnames_list[[i]]))
    pdf(file = paste(plot_path,plotnames_list[[i]],'.pdf',sep=''),width =10, height = 10)
    plot(shape_temp_all)
    dev.off()
    st_crs(shape_total) = st_crs(shape_temp_all) # added this line to avoid "arguments have different crs" error message
    shape_total = bind_rows(shape_total,shape_temp_all)
  }
  shape_total = st_make_valid(shape_total)
  return(shape_total)
}

# The following functions come from teh file 'nb_data_utils.R', created by Mitzi Morris. 
# Source: https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R

# check that graph is fully connected 
isDisconnected = function(x) {
  return(n.comp.nb(x)[[1]] > 1);
}

# test if nb object is fully connected
testconnected <- function(nb_object) {
  if (!isDisconnected(nb_object)) {
    print("Success! The Graph is Fully Connected")
  } else{
    warning("Failure... Some parts of the graph are still disconnected...")
  }
}

# validateNb
# check that nbObject is symmetric, has connected components
validateNb = function(x) {
  if (is.symmetric.nb(x) && n.comp.nb(x)[[1]] < length(x)) return(TRUE);
  return(FALSE);
}

# indexByComponent
#
# input: vector of component ids
# returns: vector of per-component consecutive node ids
#
indexByComponent = function(x) {
  y = x;
  comps = as.matrix(table(x));
  num_comps = nrow(comps);
  for (i in 1:nrow(comps)) {
    idx = 1;
    rel_idx = 1;
    while (idx <= length(x)) {
      if (x[idx] == i) {
        y[idx] = rel_idx;
        rel_idx = rel_idx + 1;
      }
      idx = idx + 1;
    }
  }
  return(y);
}

# nb2graph
#
# input: nb_object
# returns: dataframe containing num nodes, num edges,
#          and a list of graph edges from node1 to node2.
#
nb2graph = function(x) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (x[[i]][1] != 0) {
      n_links = n_links + length(x[[i]]);
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (x[[i]][1] > 0) {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i;
          node2[idx] = n2;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}


# nb2subgraph
# for a given subcomponent, return graph as lists of node1, node2 pairs
#
# inputs:
# x: nb object
# c_id: subcomponent id
# comp_ids: vector of subcomponent ids
# offsets: vector of subcomponent node numberings
# returns: list of node1, node2 ids
#
nb2subgraph = function(x, c_id, comp_ids, offsets) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        n_links = n_links + length(x[i]);
      }
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        for (j in 1:length(x[[i]])) {
          n2 = unlist(x[[i]][j]);    
          if (i < n2) {
            idx = idx + 1;
            node1[idx] = offsets[i];
            node2[idx] = offsets[n2];
          }
        }
      }
    }
  }
  return (list("node1"=node1,"node2"=node2));
}


# orderByComp
# given nbObject, reorder nb object so that all components are contiguous
# singletons moved to end of list
# returns list containing:
#   - new nbObject
#   - vector which maps old index to new index
#
orderByComponent = function(x) {
  if (!validateNb(x)) return(list());
  N = length(x);
  if (!isDisconnected(x)) return(list(x,seq(1:N)));
  
  rMap = rep(integer(0),N);
  comp_ids = n.comp.nb(x)[[2]];
  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  comp_sizes = as.vector(comps[,1]);
  idx = 1;
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] > 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] == 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  new_ids = vector("character", length=N);
  for (i in 1:N) {
    idx_old = rMap[i];
    new_ids[i] = attributes(x)$region.id[idx_old];
  }
  
  # generate new nb list
  new_nb = structure(vector("list", length=N),class="nb");
  attr(new_nb, "region.id") = new_ids;  
  attr(new_nb, "type") = attributes(x)$type;
  attr(new_nb, "sym") = attributes(x)$sym;
  attr(new_nb, "region.id") = new_ids;  
  for (i in 1:N) {
    idx_old = rMap[i];
    old_nbs = x[[idx_old]];
    num_nbs = length(old_nbs);
    new_nb[[i]] = vector("integer", length=num_nbs);
    for (j in 1:num_nbs) {
      old_id = old_nbs[j];
      if (old_id == 0) {
        new_nb[[i]][j] = as.integer(0);
      } else {
        new_id = which(rMap == old_id);
        new_nb[[i]][j] = as.integer(new_id);
      }
    }
  }
  return(list(new_nb,rMap));
}  


# reorderVector
#
# input: data vector, offsets vector
# returns: vector of same length as input data vector
#          reordered according to offsets
#
reorderVector = function(x, rMap) {
  if (!is.vector(x)) return(NULL);
  N = length(x);              
  result = vector("numeric", length=N);
  for (i in 1:N) {
    result[i]= x[rMap[i]];
  }
  return(result);
}


# reorderMatrix
#
# input: data matrix, offsets vector
# returns: matrix of same shape as input data matrix,
#          rows reordered according to offsets
#
reorderMatrix = function(x, rMap) {
  if (!is.matrix(x)) return(NULL);
  N = nrow(x);
  result = matrix(nrow=N, ncol=ncol(x));
  for (i in 1:N) {
    result[i,]= x[rMap[i],];
  }
  return(result);
}


# scale_nb_components
#
# input: nb_object
# returns: vector of per-component scaling factor (for BYM2 model)
# scaling factor for singletons is 0
#
scale_nb_components = function(x) {
  N = length(x);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = indexByComponent(comp_ids);
  
  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  scales = vector("numeric", length=num_comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    scales[i] = 0.0;
    if (N_subregions > 1) {
      # get adj matrix for this component
      drops = comp_ids != i;
      nb_tmp = droplinks(x, drops);      
      nb_graph = nb2subgraph(nb_tmp, i, comp_ids, offsets);
      adj.matrix = sparseMatrix( i=nb_graph$node1, j=nb_graph$node2, x=1, dims=c(N_subregions,N_subregions), symmetric=TRUE);
      # compute ICAR precision matrix
      Q =  Diagonal(N_subregions, rowSums(adj.matrix)) - adj.matrix;
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N_subregions) * max(diag(Q)) * sqrt(.Machine$double.eps)
      # Compute the diagonal elements of the covariance matrix subject to the 
      # constraint that the entries of the ICAR sum to zero.
      Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N_subregions),e=0))
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor = exp(mean(log(diag(Q_inv))))
      scales[i] = scaling_factor;
    }
  }
  return(scales);
}

# explicit function for calculating scale-factor 
# (alternative to INLA's scaling calculator, not recommended but useful to see how the scaling factor is calculated)
# get_scaling_factor = function(nbs) {
#  #Build the adjacency matrix
#  adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
#  #The ICAR precision matrix (note! This is singular)
#  Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix
#  
#  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
#  Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
#  
#  # Compute the diagonal elements of the covariance matrix subject to the
#  # constraint that the entries of the ICAR sum to zero.
#  #See the function help for further details.
#  Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))
#  
#  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
#  return((mean(log(diag(Q_inv)))))
#}

# # # extremeR main functions 

# Function to prep the data;
# currently only supports single pi - for muultiple pi's you need to munge the data separately (see Bird's eye example)
data.prep = function(shape, # a shapefile 
                     survey, # data which can be connected to the shape file by some small-area id;
                     shape_large.area_id_name= NA,shape_large.area_id_num = NA, # large area identifiers in the shapefile - this helps augment the survey data via merging;
                     shape_small.area_id_name = NA,shape_small.area_id_num = NA, # small area identifiers in the shapefile - these are used to match survey and shape to get proper factor levels;
                     survey_small.area_id_num = NA,survey_small.area_id_name = NA, # small area identifiers in the survey;
                     drop.incomplete.records = NA, # should the function return complete data ? uses simple complete.cases function;
                     colnames_X = NA, #  a character vector of covariates definining the design matrix X. Note these should be numeric;
                     interactions_list = NA, # a list object where each element is a string of the form "a*b" where a and be are the names of two variables in colnames_X ; 
                     # the names of each element correspond to the name of the interaction variable in X; if no names are specified, the names of a and b will be used to define the variable name ;
                     scale_X = NA,# takes values "1sd" or "2sd" corresponding namely to (x-xbar)/sd(x) and (x-xbar)/(2sd(x)). No scaling if NA ;
                     colname_y = NA,# variable name for the outcome variable - note this should be numeric ;
                     contamination = T,# the function spits out an offset based on log(P1/P0) (see Rota et al. 2013) - should this offset account for contamination ? ;
                     pi = NA, # a scalar defining the prevalence of the outcome in the population of interest ;
                     # pi_large_area = NA,# a vector of length = n. large areas in shapefile, indicating prevalence at the large-area level;
                     large_area_shape = F# should the function return a large-area shapefile ? useful for plots ;
                     ){ 
  
  if(all(is.na(c(shape_large.area_id_name,shape_large.area_id_num)))){warning("no large area ID (nominal or numeric) to search for in shapefile - preparing data exclusively at the small-area level")}
  
  if(!is.na(shape_large.area_id_name)){# factor the large area english names in shape file 
    shape[[shape_large.area_id_name]] <- as.factor(shape[[shape_large.area_id_name]])
  }
  if(!is.na(shape_large.area_id_num)){# factor the large area english names in shape file 
    shape[[shape_large.area_id_num]] <- as.factor(shape[[shape_large.area_id_num]])
  }
  
  
  
  if(all(is.na(c(shape_small.area_id_name,shape_small.area_id_num)))){stop("no small area ID (nominal or numeric) to search for in shapefile")}
  if(all(is.na(c(survey_small.area_id_name,survey_small.area_id_num)))){stop("no small area ID (nominal or numeric) to search for in survey data")}
  
  if(!is.na(survey_small.area_id_name) ){
    if( !survey_small.area_id_name %in% colnames(survey)){stop("the small area nominal ID provided does not exist in survey data")}
    if( !shape_small.area_id_name %in% colnames(shape)){stop("the small area nominal ID provided does not exist in shapefile")}
    # cleaning name id - if empty, give NA value, then factor 
    survey[[survey_small.area_id_name]] = as.factor(ifelse(survey[[survey_small.area_id_name]]=="",NA,survey[[survey_small.area_id_name]]))
  }
  if(!is.na(survey_small.area_id_num)){
    if( !survey_small.area_id_num %in% colnames(survey)){stop("the small area numerical ID provided does not exist in survey data")}
    if( !shape_small.area_id_num %in% colnames(shape)){stop("the small area numerical ID provided does not exist in shapefile")}
    # cleaning num id- if empty, give NA value, then factor 
    survey[[survey_small.area_id_num]] = as.factor(ifelse(survey[[survey_small.area_id_num]]=="",NA,survey[[survey_small.area_id_num]]))
  }
  
  
  # merge to add legitimate large-area to survey and/or ensure small-area IDs match shapefile
  index <- c(shape_small.area_id_name = shape_small.area_id_name,
             shape_small.area_id_num = shape_small.area_id_num,
             shape_large.area_id_name = shape_large.area_id_name,
             shape_large.area_id_num = shape_large.area_id_num)
  
  index = index[which(!is.na(index))]
  
  by.x = c(num = survey_small.area_id_num,name = survey_small.area_id_name)
  by.y = c(num = shape_small.area_id_num,name = shape_small.area_id_name)
  
  by.y = by.y[which(!is.na(by.x))]
  by.x = by.x[which(!is.na(by.x))]
  
  if(length(by.y)==0 | length(by.x)==0){stop("we cannot match numbers with names - if you provide a numerical/nominal id for the survey, you must also provide an equilvalent for the shapefile, and vice-versa") }
  
  
  survey = merge(survey,
                 as.data.table(shape[,index])[,1:length(index)],
                 by.x = by.x,
                 by.y =  by.y,
                 all.x = TRUE)

  # keep only relevant columns 
  survey = survey[,c(by.x,index[-which(index %in% by.y)],colnames_X,colname_y),with=F]
  # and columns in X must all be numeric
  if(!all(survey[,lapply(.SD,class),.SDcols = c(colnames_X,colname_y)] == "numeric" | survey[,lapply(.SD,class),.SDcols = c(colnames_X,colname_y)] == "integer" )){stop("one of your variables in colnames_X or colname_y is not numeric - please change this in the survey data and re-run the data.prep function")}
  
  # if there are any interactions
  names(interactions_list) = ifelse(names(interactions_list)=="",gsub("\\*","\\_",unlist(interactions_list)),names(interactions_list))
  if(length(interactions_list)>0){
    # make sure all variables needed are in X 
    if(!all(unlist(strsplit(unlist(interactions_list),split = "\\*")) %in% colnames_X)){stop("some variables needed for interactions are missing from colnames_X - please include all interaction variables in colnames_X")}
    # calculate the interactions
    int.temp = data.table()
    for(i in 1:length(interactions_list)){
      int.temp = cbind(int.temp,
                       apply(survey[,unlist(strsplit(interactions_list[[i]],split = "\\*")),with = F],1,function(x){x[1]*x[2]}))
    }
    names(int.temp) = names(interactions_list)
    # and add them to the survey
    survey = as.data.table(cbind(survey,int.temp))
  }
  
  # remove incomplete records
  if(drop.incomplete.records){
    if(sum(complete.cases(survey))>0){
      survey_complete = survey[complete.cases(survey),]
    }else{stop("cannot subset complete cases - a missing value exists for every record\nmaybe something went wrong with the matching - check that levels of shape and survey areas match (spelling could be a problem)")}
  }else{
    survey_complete = survey
  }
  
  # re-factor now that we have removed data 
  if(!is.na(survey_small.area_id_name)){
    survey_complete[[ survey_small.area_id_name]] = as.factor(as.character(unlist(survey_complete[[ survey_small.area_id_name]])))
  }
  if(!is.na(survey_small.area_id_num)){
    survey_complete[[ survey_small.area_id_num]] = as.factor(as.character(unlist(survey_complete[[ survey_small.area_id_num]])))
  }
  
  # and extract area ids
  small_area_id = match(survey_complete[[by.x[1]]],shape[[by.y[1]]])
  small_area_names = shape[[by.y]]
  
  if(!all(is.na(c(shape_large.area_id_name,shape_large.area_id_num)))){
    by.y.large = c(num = shape_large.area_id_num,name = shape_large.area_id_name)
    by.y.large = by.y.large[which(!is.na(by.y.large))]
    #  large area is typically added on from shape from small-area
    large_area_id =  match(survey_complete[[by.y.large[1]]],levels(shape[[by.y.large[1]]]))
    large_area_names = levels(shape[[by.y.large[1]]])
    N_large_area = nlevels(shape[[by.y.large[1]]])
  }
  
  
  # get covariates in a matrix object 
  X = survey_complete[,c(colnames_X,names(interactions_list)),with=F]
  mu_X = apply(X,2,mean,na.rm=T)
  sd_X = apply(X,2,sd,na.rm=T)
  # choose and apply standardizaton type
  if(scale_X=="1sd"){
    X = (X - t(array(mu_X,c(dim(X)[2],dim(X)[1]))))/t(array(sd_X,c(dim(X)[2],dim(X)[1])))
  }
  if(scale_X=="2sd"){
    X_dich = which(sapply(X,function(x){all(sort(unique(x))==c(0,1))}))
    X.tmp = (X - t(array(mu_X,c(dim(X)[2],dim(X)[1]))))/t(array(2*sd_X,c(dim(X)[2],dim(X)[1])))
    for(i in 1:length(X_dich)){
    X.tmp[[X_dich[i]]] = X[[X_dich[i]]]
    }
    X = X.tmp
  }
  X = cbind(intercept = 1,X)
  p = dim(X)[2]
  
  # get outcome
  if(!is.na(colname_y)){
    if(!all(sort(unique(survey_complete[[colname_y]])) == c(0,1)) | !is.numeric(survey_complete[[colname_y]])){stop("outcome must be a numeric, dichotomous variable")}
    y = survey_complete[[colname_y]]
  }else{stop("please enter a string for colname_y so we can identify the outcome in the surey data")}
  
  # number of observations
  n = length(y)
  # number of cases
  n1 = sum(y)
  # number of controls
  n0 = n-n1
  
  if(n1==n |n0==n){stop("no variance in outcome variable")}
  
  # Now we calculate the offset :
  
  # offset = log(P1/P0) =
  #        = log( n1/(pi*N) / n0/((1-pi)*N) = 
  #        = log( n1/pi / n0/(1-pi) )
  
  # under contamination: 
  # offset = log( (n1 + pi*n0)/(pi*N)) / ((n0-pi*n0)/((1-pi)*N)) ) = 
  #          = log( (n1 + pi*n0)(1-pi)/(pi*(n0-pi*n0)) ) = 
  #          = log( (n1 + pi*n0)(1-pi)/pi*n0(1-pi) ) = 
  #          = log( (n1 + pi*n0)/pi*n0 ) = 
  #          = log( n1/pi*n0 + 1 )
  
  if(contamination){
    theta_0_large_area = NA; theta_1_large_area = NA ; n1_large_area = NA; n0_large_area = NA; offset_large_area = NA
    theta_0 = NA; theta_1 = NA ;
    # calculate offset and mixture probabilities under contaminationa: 
    # at the area-level
    # if(!is.na(pi_large_area)){
    #   n1_large_area = as.numeric(table(survey_complete[[by.y.large[1]]])[,"1"])
    #    n0_large_area = as.numeric(table(survey_complete[[by.y.large[1]]])[,"0"])
    #    offset_large_area = c(); for(i in 1:length(n1_large_area)){ offset_large_area = c(offset_large_area, log( n1_large_area[i]/(pi_large_area[i]*n0_large_area[i]) + 1 ) ) }
    #    theta_0_large_area = 0
    #    theta_1_large_area = n1_large_area/(n1_large_area+pi_large_area*n0_large_area)
    #  }else{warning("no large-area level prevalence was defined - will not calculate large-area level contaminated offest")} 
    # at the population level 
    if(!is.na(pi)){
      offset = log( n1/(pi*n0) + 1 ) 
      theta_0 = 0
      theta_1 = n1/(n1+pi*n0)
    }else{warning("no population-level prevalence was defined - will not calculate contaminated offest")}
    if(is.na(theta_0)){stop("could not calculate contaminated offset - likely missing prevalence pi")}
    
    # if(is.na(theta_0) & is.na(theta_0_large_area)){stop("could not calculate contaminated offset - likely missing prevalence (pi or pi_large_area)")}
  }else{
    # calculate offset under simple case-control (King and Zeng, 2001, 'prior correction'): 
    # at the area-level 
    #if(!is.na(pi_large_area)){
    #   offset_large_area = c(); for(i in 1:length(n1_large_area)){ offset_large_area = c(offset_large_area, log( (n1_large_area[i]/pi_large_area[i])/(n0_large_area[i]/(1-pi_large_area[i]))  ) ) }
    # }else{warning("no large-area level prevalence was defined - will not calculate large-area level offest")} 
    # at the population level
    if(!is.na(pi)){
      offset = log( (n1/pi)/(n0/(1-pi))  )
    }else{warning("no population-level prevalence was defined - will not calculate offest")}
  }
  
  
  # now extract neighbourhood object 
  # ensure you can take lon-lat centroids over weird spherical objects  
  sf::sf_use_s2(FALSE)
  # get distance matrix 
  lat_lon = st_coordinates(st_centroid(shape))
  colnames(lat_lon) = c("lat","lon")
  D = geodist(lat_lon,measure = "geodesic") #1m
  # fully connect the graph based on shortest distance 
  nb = addnbs(sp.sample = shape,ID = shape[[by.y]],D=D);
  # extract graph information
  nb_graph = nb2graph(nb)

  if(length(grep("large",names(index)))>0){
    # get a shapefile with areas collapsed into large areas - useful to plot large-area effects 
    shape_large = shape
    shape_large$large_area_id = shape_large[[index [grep("large",names(index))[1]]]]
    shape_large = 
      shape_large %>% 
      group_by(large_area_id) %>%
      summarise(geometry = sf::st_union(geometry)) %>%
      ungroup()
    names(shape_large)[1] = index [grep("large",names(index))[1]]
  }else{shape_large = NA}
  
  output = list(survey_complete = survey_complete,
                
                n = n,
                y = y,
                
                pi = pi,
                
                n1 = n1, 
                n0 = n0,
                offset = offset,
                theta_1 = theta_1,
                theta_0 = theta_0,
                #offset_large_area = offset_large_area,
                #theta_1_large_area = theta_1_large_area,
                #theta_0_large_area = theta_0_large_area,
                
                X = X,
                p = p,
                mu_X = mu_X,
                sd_X = sd_X,
                scale_X = scale_X,
                
                small_area_id = small_area_id,
                small_area_names = small_area_names,
                large_area_id =  large_area_id,
                large_area_names = large_area_names,
                
                N_large_area = N_large_area,
                nb_object = nb,
                N_small_area = nb_graph$N,
                node1_small_area = nb_graph$node1,
                node2_small_area = nb_graph$node2,
                N_small_area_edges = nb_graph$N_edges,
                scaling_factor = scale_nb_components(poly2nb(shape,queen = TRUE))[1],
                nb.matrix  = nb2mat(neighbours = nb ,style = 'B',zero.policy = TRUE),
                
                shape_small = shape,
                shape_large = shape_large
  )
  
  return(output)
}

# Function to plot nodes and edges (if you want, over a map from a shapefile)
plot_network <- function(shape, nb_object,
                         plot.over.map = T,
                         map.border.col = adjustcolor('grey',0.5), map.plygon.col = adjustcolor('lightcoral',0.10),map.border.lwd = 0.1,
                         node.cex=0.01,edge.lwd=0.25,node.pch=0,edge.lty=1,edge.col = adjustcolor('blue',0.25)){
  if(plot.over.map){connectedplot = plot(st_geometry(shape),border = map.border.col,lwd = map.border.lwd,col = map.plygon.col)}
  connectedplot = plot(nb_object,st_coordinates(st_centroid(shape)),
                       add=plot.over.map,col = edge.col,cex = node.cex,lwd = edge.lwd,pch = node.pch,lty = edge.lty)
}

# Function to simulate data from a relevant data generating process (DGP) ;
# currently this function only supports the creation of DGPs with 1 layer of area effects (i.e. small-area effects); 
# we plan to include generation capabilities for hierarchical areal (i.e. small-area & large-area) processes in the near future. 
dgf = function(
  n.sims = 1, # how many samples of simulated data would you like ?
  n = 100, # how large (sample size) should each sample be ? 
  pi.hat.naive = 0.5, # what should be the fraction of cases in the sample ? 
  p = 1, # how many normally-distributed covariates should the DGP have ?
  X_corr = 0, # what should be the average correlation among these covariates ? 
  pi = 0.05, # what is the population-level probability of sampling a case ? 
  Moran.I.corr = 0.8, # what degree of global spatial autocorrelation (Moran I) should the underlying DGP have ?
  spatial_structure="scotland_lipcancer" 
  # on which map should the data be simulated ? - choose between scotland_lipcancer, pennsylvania_lungcancer and newyork_leukemia)
  ){
  
  # some errors 
  if(p<1){stop('number of regression coefficients excluding the intercept, `p`, must be at least 1')}
  if(any(X_corr < 0 | X_corr > 1)) stop('correlation between covariates, `X_corr`, must be between 0 and 1')
  if(!spatial_structure %in% c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia")){stop('we do not have the map you have specified\npelase choose between scotland_lipcancer, pennsylvania_lungcancer and newyork_leukemia')}
  
  
  # define some observable quantities to calculate for the process to be generated
  n1 = n*pi.hat.naive
  n0 = n-n1
  P1 = (n1+pi*n0)/pi
  P0 = n0
  theta_0 = 0
  theta_1 = n1/(n1+pi*n0)

  # simulate from a CAR pocess on an existing map
  
  sf::sf_use_s2(FALSE)
  if(spatial_structure=="scotland_lipcancer"){
    temp.sp= scotland$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(scotland$spatial.polygon)
  }
  if(spatial_structure=="pennsylvania_lungcancer"){
    temp.sp= pennLC$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(pennLC$spatial.polygon)
    temp.sf$y = as.data.table(pennLC$data)[,lapply(.SD,sum),by = "county",.SDcols = c("cases")]$cases
  }
  if(spatial_structure=="newyork_leukemia"){
    temp.sp = NYleukemia$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(NYleukemia$spatial.polygon)
    temp.sf$y = NYleukemia$data$cases
  }
  # get distance matrix 
  lat_lon = st_coordinates(st_centroid(temp.sf))
  colnames(lat_lon) = c("lat","lon")
  D = geodist(lat_lon,measure = "geodesic") 
  # ensure graph is fully connexted 
  nb = addnbs(sp.sample = temp.sf,ID = temp.sf$small_area,D=D);
  # Have we achieved a fully connected graph ? 
  if(isDisconnected(nb)==FALSE){print("Success! The Graph is Fully Connected")}else{print("Failure... Some parts of the graph are still disconnected...")}
  # get neighborhood objects - these will be useful later for the fitting of the model
  nbs = nb2graph(nb);
  
  nb_objects = list(
    N_small_area = nbs$N,
    node1_small_area = nbs$node1,
    node2_small_area = nbs$node2,
    N_small_area_edges = nbs$N_edges,
    scaling_factor = scale_nb_components(poly2nb(temp.sf,queen = TRUE))[1],
    nb  = nb2mat(neighbours = nb ,style = 'B',zero.policy = TRUE),
    sp.object = temp.sp
  )
  # # # simulate conditionally autoregressive effects 
  jags.model = "model{
  phi ~ dmnorm(zero, tau * (D - alpha*W));
  tau ~ dgamma(1,1); 
  }"
  data = list(W = nb_objects$nb,
              D = diag(rowSums(nb_objects$nb)),
              zero = rep(0,nb_objects$N_small_area),
              alpha = 0.99999)
  CAR.sample = 
    jags(data=data,
         parameters.to.save=c("phi",'tau'),
         model.file=textConnection(jags.model),
         n.iter=100,
         DIC = F)
  # # # pick simulated data closest to the wanted Moran I 
  sample.Moran.I = apply(CAR.sample$BUGSoutput$sims.list$phi,1,function(x){ape::Moran.I(x,nb_objects$nb)}$observed)
  temp.sf$y = as.numeric(scale(CAR.sample$BUGSoutput$sims.list$phi[which(abs(sample.Moran.I-Moran.I.corr)==min(abs(sample.Moran.I-Moran.I.corr))),]))

  # scale the effects for cross-simulation consistency
  gamma = as.numeric(scale(temp.sf$y))
  small_area_id = as.integer(as.factor(round(runif(n = n,min = 1,max = nb_objects$N_small_area))))
  
  # draw covariates and covariate coefficients
  Sigma = array(1,c(p,p))*X_corr
  diag(Sigma) = 1
  X = cbind(1,mvrnorm(n = n,mu = rep(0,p) ,Sigma = Sigma))
  beta = sapply(X = 1:p,function(x){rnorm(n = 1,mean = 0,sd = 1)})
  # calculate the total number of regression coefficients, including the intercept
  p = dim(X)[2]
  
  # calibrate intercept of linear DGP to match pi.hat.naive
  intercept.optimizer = 
    optimize(f = function(beta0_true){
      
      beta_temp = c(beta0_true,beta)
      mu = log(P1/P0) + X %*% beta_temp +  gamma[small_area_id]
      rho = inv_logit(mu)
      
      ybar_gen = (1-rho)*theta_0 + rho*theta_1
      temp = pi.hat.naive - mean(ybar_gen)
      
      return(abs(temp))
    },
    interval = c(-30,30),
    tol =  1/10000000000)
  
  # after calibration, generate sample 
  beta_temp = c(intercept.optimizer $minimum,beta)
  mu = log(P1/P0) + X %*% beta_temp +  gamma[small_area_id]
  rho = inv_logit(mu)
  
  y_gen = sapply(1:n.sims,function(x){
    r_gen = rbinom(n = n, size = 1, prob = rho);
    y_gen = ifelse(r_gen==0,
                   rbinom(n = n,size = 1,prob = theta_0),
                   rbinom(n = n,size = 1,prob = theta_1));
    return(y_gen)
  } )
  
  # generate output object with samples and ground-truth for parameters 
  output = list(y = y_gen, 
                pi = pi,
                n = n, 
                p = p, # remember that this now includes intercept
                X = X,
                X_corr = X_corr,
                small_area_id = small_area_id,
                beta = beta_temp,
                gamma = gamma,
                mu = mu,
                moranI.mu = ape::Moran.I(as.numeric(mu),nb_objects$nb[small_area_id,small_area_id]),
                moranI.gamma = ape::Moran.I(as.numeric(gamma),nb_objects$nb)
  )
  
  output = append(output,nb_objects)
  
  return(output)
}
# # # output: a list object with the following variables: 
    # y = see input ; 
    # pi = see input ; 
    # n =  see input ; 
    # p =  see input ; 
    # X =  matrix of simulated covariates, including a vector of 1s for the intercept ; 
    # X_corr = see input ; 
    # small_area_id = see input for 'fit' function ; 
    # beta = true regression coefficient values ;
    # gamma = true small-area effects ;
    # mu = true latent propensity of recruitment ; 
    # moranI.mu = empirical Moran I at the individual-level for true latent propensity of recruitment ;
    # moranI.gamma = empirical Moran I for true small-area effects 

# # # Function to analyse generated data. 
# # # Note this function uses STAN, and the working environment needs to be optimized for STAN fitting purposes, in order to efficiently run the model.  
fit = function(data = data_sim, 
               # data is a list object containing ad-minimum the following elements: 
               # data = 
               # list(n = integer denoting sample size (recommended n>=100) ;
               #      y = vector of size n of Bernoulli distributed values y ∈ {0,1} ;
               #      pi = integer denoting the % of cases in the population pi ∈ [0,1] ;
               #      p = number of fixed effects regression coefficients (including the intercept) ;
               #      X = (n x p) matrix of covariates (first column is a column of 1s for intercept) - ensure columns of X are scaled for optimal performance ;
               # The following are relevant if there are supposed to be any small-area-level effects :     
               #      N_small_area = integer denoting the number of small-areas in the map;
               #      small_area_id = numerical vector (length = n ) of integers denoting the numerical ID of the area to which each individual belongs to - small_area_id ∈ {0,...,N_small_area};
               # The following are relevant if there are supposed to be BYM2-style small-area-level effects : 
               # to get these, simply use the nb2graph function, followed by the scale_nb_components to derive the scaling factor
               #      node1_small_area = ...;
               #      node2_small_area = ...;
               #      N_small_area_edges = ...;
               #      scaling_factor = ...;
               # The following are relevant if there are supposed to be any large-area-level effects :     
               #      N_large_area = integer denoting the number of large-areas in the map;
               #      large_area_id = numerical vector (length = n ) of integers denoting the numerical ID of the area to which each individual belongs to - large_area_id ∈ {0,...,N_large_area};
               # )
               show_code = T, # should the stan code print on the screen at the end of the run ? (useful for checks)
               contamination = F, # should the model include the rota et al. (2013) style contamination layer ?
               offset = F, # should the model include a King and Zeng (2001) style offset ? if contamination is also specified, this will be a contaminated-offset ;
               beta_prior = "normal", # what prior should the regression coefficients have ? choice between normal and cauchy ;
               small_area_prior = "fixed", # what should be the small-area effects type ? choice between area fixed effects, area random effects, area ICAR effects or area BYM2 effects or NA = no area effects ;
               large_area_prior = NA, # what should be the small-area effects type ? choice between area fixed effects, area random effects, or NA = no area effects ;
               intercept_scale = NA,
               # stan functionalities - arguments get passed to stan; see the stan function for an explenation ;
               iter = 100, 
               warmup = 50, 
               thin = 4, 
               cores = 4, 
               chains = 4, 
               control = list(max_treedepth = 7, adapt_delta = 0.8), 
               verbose = TRUE ){

  model_skeleton = "data{
int<lower = 1> n; // total number of observations
int<lower = 1> p; // number of covariates in design matrix
int<lower = 0> y[n]; // vector of labels
matrix[n, p] X; // design matrix
}//end_data

parameters{//begin_param

}//end_param

transformed parameters{//begin_transform
vector[n] mu;
mu = X * beta; 
}//end_transform

model{//begin_model

// likelihood
y  ~  bernoulli_logit(mu);
}//end_model

generated quantities{//begin_gen
vector[n] log_lik;
vector[n] y_gen;
  
for (i in 1:n) {
log_lik[i] = bernoulli_logit_lpmf(y[i] | mu[i]);
y_gen[i] = bernoulli_logit_rng(mu[i]);
} }//end_gen
";
  
  # # # parameters to monitor
  pars = c("log_lik","y_gen","mu")
  
  # # # define beta prior
  
  if(beta_prior=="normal"){
    
    pars = c(pars,"beta")
    
    model_skeleton = 
      gsub("\\}//end_param",
           "vector[p] beta;\n}//end_param",
           model_skeleton)
    
    if(!is.na(intercept_scale)){
    model_skeleton = 
      gsub("// likelihood",
           paste("beta[1] ~ normal(0,",intercept_scale,");\n//beta ~ normal(0,1);\n// likelihood",sep=""),
           model_skeleton)
    }else{
      model_skeleton = 
        gsub("// likelihood",
             "beta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    }
    
  };
  # cauchy by gamma: https://betanalpha.github.io/assets/case_studies/fitting_the_cauchy.html
  if(beta_prior=="cauchy"){
    
    pars = c(pars,"beta")
    
    model_skeleton = 
      gsub("begin_param",
           "begin_param\nvector[p] aux_a;\nvector<lower = 0>[p] aux_b;\n",
           model_skeleton)
    
    model_skeleton = 
      gsub("begin_transform",
           "begin_transform\nvector[p] beta = aux_a ./ sqrt(aux_b);\n",
           model_skeleton)
    
    if(!is.na(intercept_scale)){
      model_skeleton = 
        gsub("// likelihood",
             paste("aux_a ~ normal(0,1);\naux_b[1] ~ gamma(0.5,",intercept_scale^2,"*0.5);\naux_b ~ gamma(0.5,0.5);\n// likelihood",sep=""),
             model_skeleton)
    }else{
      model_skeleton = 
        gsub("// likelihood",
             paste("aux_a ~ normal(0,1);\naux_b ~ gamma(0.5,0.5);\n// likelihood",sep=""),
             model_skeleton)
    }
  };
  
  # # # define small_area effects prior
  if(!is.na(small_area_prior)){
    if(small_area_prior=="fixed"){
      
      pars = c(pars,"gamma")
      
      data$Z = model.matrix(~as.factor(data$small_area_id)-1)
      colnames(data$Z) = gsub("as\\.factor\\(data\\$small_area_id\\)","",colnames(data$Z))
      for(i in 1:data$N_small_area){
        if(! as.character(i) %in% colnames(data$Z)){
          temp = data.table(x =  rep(0,data$n))
          names(temp) = as.character(i)
          data$Z = cbind(data$Z,temp)
        }
      }
      data$Z = as.matrix(data$Z)[,order(as.numeric(colnames(data$Z)))]
      data$pZ = dim(data$Z)[2]
      
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> pZ;\nmatrix[n, pZ] Z;\nint<lower = 1> N_small_area;}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[pZ] gamma;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "mu = Z * gamma +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "gamma ~ normal(0,1);\nsum(gamma) ~ normal(0, 0.01 * N_small_area);// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="random"){
      
      pars = c(pars,"gamma","sigma_gamma")
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];\nint<lower = 1> N_small_area; \n}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[N_small_area] gamma;\n real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "mu = gamma[small_area_id]*sigma_gamma +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "gamma ~ normal(0,1);\nsigma_gamma ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="BYM2"){
      
      pars = c(pars,"gamma","sigma_gamma","lambda","phi","psi")
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];
int<lower = 1> N_small_area;
int<lower = 1> N_small_area_edges;
int<lower=1, upper=N_small_area> node1_small_area[N_small_area_edges];
int<lower=1, upper=N_small_area> node2_small_area[N_small_area_edges];
real scaling_factor;
}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[N_small_area] phi;
vector[N_small_area] psi;
real<lower = 0,upper = 1> lambda;
real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "vector[N_small_area] gamma = (sqrt(1-lambda) * phi + sqrt(lambda / scaling_factor) * psi)*sigma_gamma;\n
mu = gamma[small_area_id] +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "target += -0.5 * dot_self(psi[node1_small_area] - psi[node2_small_area]);
phi ~ normal(0,1); 
sum(psi) ~ normal(0, 0.01 * N_small_area);
lambda ~ beta(0.5,0.5);
sigma_gamma ~ normal(0,1);
// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="ICAR"){
      
      pars = c(pars,"gamma","sigma_gamma")
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];
int<lower = 1> N_small_area;
int<lower = 1> N_small_area_edges;
int<lower=1, upper=N_small_area> node1_small_area[N_small_area_edges];
int<lower=1, upper=N_small_area> node2_small_area[N_small_area_edges];
real scaling_factor;
}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[N_small_area] psi;
real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "vector[N_small_area] gamma = psi*sigma_gamma;\n
mu = gamma[small_area_id] +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "target += -0.5 * dot_self(psi[node1_small_area] - psi[node2_small_area]);
sum(psi) ~ normal(0, 0.01 * N_small_area);
sigma_gamma ~ normal(0,1);
// likelihood",
             model_skeleton)
    };
  }
  # # # define large_area effects prior
  if(!is.na(large_area_prior)){
    if(large_area_prior=="fixed"){
      
      pars = c(pars,"eta")
      
      data$Q = model.matrix(~as.factor(data$large_area_id)-1)
      colnames(data$Q) = gsub("as\\.factor\\(data\\$large_area_id\\)","",colnames(data$Q))
      data$pQ = dim(data$Q)[2]
      
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> pQ;\nmatrix[n, pQ] Q; \n}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[pQ] eta;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "mu = Q * eta +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "eta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
    if(large_area_prior=="random"){
      
      pars = c(pars,"eta","sigma_eta")
      
      model_skeleton = 
        gsub("\\}//end_data",
             "int<lower = 1> large_area_id[n];\nint<lower = 1> N_large_area; \n}//end_data",
             model_skeleton)
      
      model_skeleton = 
        gsub("\\}//end_param",
             "vector[N_large_area] eta;\n real<lower = 0> sigma_eta;\n}//end_param",
             model_skeleton)
      
      model_skeleton = 
        gsub("mu =",
             "mu = eta[large_area_id]*sigma_eta +",
             model_skeleton)
      
      model_skeleton = 
        gsub("// likelihood",
             "eta ~ normal(0,1);\nsigma_eta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
  }
  # # # define mu with offset
  if(!is.null(data$pi_large_area) & any(is.na(data$pi_large_area))){warning("some missing values in pi_large_areas - cannot be used to make offset")}
  if(offset & all(is.null(data$pi_large_area))){
    
    if(contamination){ data$log_offset = log(sum(data$y)/(data$pi*(data$n-sum(data$y))) +1) }else{
      data$log_offset = log(((1-data$pi)/data$pi)*(mean(data$y)/(1-mean(data$y)))) }
    
    model_skeleton = 
      gsub("\\}//end_data",
           "real log_offset;\n}//end_data",
           model_skeleton)
    
    model_skeleton = 
      gsub("mu =",
           "mu = log_offset +",
           model_skeleton)
  };
  if(offset & all(!is.null(data$pi_large_area))){
    
    if(contamination){ 
      if(any(names(data$offset_large_area)!=data$large_area_names)){stop("Order of offsets in data doesn't match order of levels in large-area.")}
        data$log_offset = data$offset_large_area
      }else{
      if(any(names(data$offset_large_area_no.contamination)!=data$large_area_names)){stop("Order of offsets in data doesn't match order of levels in large-area.")}
        data$log_offset = data$offset_large_area_no.contamination
      }
    
    model_skeleton = 
      gsub("\\}//end_data",
           "vector[N_large_area] log_offset;\n}//end_data",
           model_skeleton)
    
    model_skeleton = 
      gsub("mu =",
           "mu = log_offset[large_area_id] +",
           model_skeleton)
  };
  # # # redefine the model with contamination
  if(contamination & all(is.na(data$pi_large_area))){
    data$theta = c()
    data$theta[1] = 0 
    data$theta[2] = sum(data$y)/(sum(data$y) + data$pi*(data$n-sum(data$y)))
    
    model_skeleton = 
      gsub("\\}//end_data",
           "vector[2] theta;\n}//end_data",
           model_skeleton)
    
    model_skeleton = 
      gsub("y  ~  bernoulli_logit\\(mu\\);",
           "for (i in 1:n) { 
  target += log_mix(1-inv_logit(mu[i]),
                    bernoulli_lpmf(y[i] | theta[1]),
                    bernoulli_lpmf(y[i] | theta[2]));
                    }",
           model_skeleton)
    
    model_skeleton = 
      gsub("vector\\[n\\] y\\_gen;",
           "vector[n] y_gen;\nint r_gen[n];",
           model_skeleton)
    
    model_skeleton = 
      gsub("log_lik\\[i\\] = bernoulli_logit_lpmf\\(y\\[i\\] \\| mu\\[i\\]\\);",
           "log_lik[i] = log_mix(1-inv_logit(mu[i]),bernoulli_lpmf(y[i] | theta[1] ),bernoulli_lpmf(y[i] | theta[2] ) );",
           model_skeleton)
    
    model_skeleton = 
      gsub("y_gen\\[i\\] = bernoulli_logit_rng\\(mu\\[i\\]\\);",
           "r_gen[i] = bernoulli_rng(inv_logit(mu[i]));\ny_gen[i] = bernoulli_rng(theta[r_gen[i]+1]);",
           model_skeleton)
    
  };
  if(contamination & all(!is.null(data$pi_large_area))){
    data$theta = as.matrix(t(data.table(theta_0 = data$theta_0_large_area, 
                                        theta_1 = data$theta_1_large_area)))
    
    model_skeleton = 
      gsub("\\}//end_data",
           "matrix[2,N_large_area] theta;\n}//end_data",
           model_skeleton)
    
    model_skeleton = 
      gsub("y  ~  bernoulli_logit\\(mu\\);",
           "for (i in 1:n) { 
  target += log_mix(1-inv_logit(mu[i]),
                    bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),
                    bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]));
                    }",
           model_skeleton)
    
    model_skeleton = 
      gsub("vector\\[n\\] y\\_gen;",
           "vector[n] y_gen;\nint r_gen[n];",
           model_skeleton)
    
    model_skeleton = 
      gsub("log_lik\\[i\\] = bernoulli_logit_lpmf\\(y\\[i\\] \\| mu\\[i\\]\\);",
           "log_lik[i] = log_mix(1-inv_logit(mu[i]),bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]) );",
           model_skeleton)
    
    model_skeleton = 
      gsub("y_gen\\[i\\] = bernoulli_logit_rng\\(mu\\[i\\]\\);",
           "r_gen[i] = bernoulli_rng(inv_logit(mu[i]));\ny_gen[i] = bernoulli_rng(theta[r_gen[i]+1,large_area_id[i]]);",
           model_skeleton)
    
  };
  
  
  # # # avoid stan recompiling of code if possible - do not amend code string
  if(exists("model_skeleton_store")){
    if(model_skeleton == model_skeleton_store){
      
    }else{model_skeleton_store <<- model_skeleton}}else{model_skeleton_store <<- model_skeleton}
  
  fit <- stan(model_code = model_skeleton_store, 
              data = data, 
              iter = iter,
              warmup = warmup,
              thin = thin,
              cores = cores,
              pars = pars, 
              chains = chains,
              control = control,
              verbose = verbose);
  
  
  cat(model_skeleton_store)
  # # # should we show the code at the end ? 
  if(show_code == T){ cat(model_skeleton_store) }
  
  return(fit)
}
# # # output: a STAN model fit object , which can be analysed with traditional STAN functions. 
