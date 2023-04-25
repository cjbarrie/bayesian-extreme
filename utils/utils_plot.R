# load some useful functions for analysis of convergence and fit 
source('utils/monitorplot.R') # https://github.com/avehtari/rhat_ess/blob/master/code/monitorplot.R
source('utils/monitornew.R') # https://github.com/avehtari/rhat_ess/blob/master/code/monitornew.R
libs <- c("ggplot2", "grid", "gridExtra", "mgcv", "caret", "RColorBrewer", "modeest", "stringr")
# Check if each library is already installed
# If not, install it
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib)
    library(lib, character.only = TRUE)
  } else {
    library(lib, character.only = TRUE)
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # CONVERGENCE PLOTS # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # look at ESS properties of worst-converging parameter (via Bulk ESS and Tail ESS)
plotConvergence_minESS <- function(mon){
which_min_bulk_ess <- sample(which(mon[, 'Bulk_ESS']==min(mon[, 'Bulk_ESS'],na.rm=T)),size = 1) # if multiple with same ESS, sample 
local_ess.min.bulk = plot_local_ess(fit = fit_object, par = which_min_bulk_ess, nalpha = 50)
quantile_ess.min.bulk = plot_quantile_ess(fit = fit_object, par = which_min_bulk_ess, nalpha = 40)
change_ess.min.bulk = plot_change_ess(fit = fit_object, par = which_min_bulk_ess)
which_min_tail_ess <- sample(which(mon[, 'Tail_ESS']==min(mon[, 'Tail_ESS'],na.rm=T)),size = 1) # if multiple with same ESS, sample 
local_ess.min.tail = plot_local_ess(fit = fit_object, par = which_min_tail_ess, nalpha = 50)
quantile_ess.min.tail = plot_quantile_ess(fit = fit_object, par = which_min_tail_ess, nalpha = 40)
change_ess.min.tail = plot_change_ess(fit = fit_object, par = which_min_tail_ess)

# # # plot convergence diagnostics for worst bulk ESS parameter (they are self-explanatory)
grid.arrange(local_ess.min.bulk,quantile_ess.min.bulk,change_ess.min.bulk,
             local_ess.min.tail,quantile_ess.min.tail,change_ess.min.tail,
             nrow = 2,ncol = 3,
             top = textGrob(paste("convergence diagnostics for lowest 'bulk ESS' param. (",names(which_min_bulk_ess),", top) and 'tail ESS' param. (",names(which_min_tail_ess),", bottom)"),
                            gp=gpar(fontsize=10,font=2)))
}
# # # check mixing properties via uniformity of ranked draws histograms
plotConvergence_mixing_rank.hist_minESS <- function(mon,breaks = 20){
  
which_min_bulk_ess <- sample(which(mon[, 'Bulk_ESS']==min(mon[, 'Bulk_ESS'],na.rm=T)),size = 1) # if multiple with same ESS, sample 
local_ess.min.bulk = plot_local_ess(fit = fit_object, par = which_min_bulk_ess, nalpha = 50)
quantile_ess.min.bulk = plot_quantile_ess(fit = fit_object, par = which_min_bulk_ess, nalpha = 40)
change_ess.min.bulk = plot_change_ess(fit = fit_object, par = which_min_bulk_ess)
which_min_tail_ess <- sample(which(mon[, 'Tail_ESS']==min(mon[, 'Tail_ESS'],na.rm=T)),size = 1) # if multiple with same ESS, sample 
local_ess.min.tail = plot_local_ess(fit = fit_object, par = which_min_tail_ess, nalpha = 50)
quantile_ess.min.tail = plot_quantile_ess(fit = fit_object, par = which_min_tail_ess, nalpha = 40)
change_ess.min.tail = plot_change_ess(fit = fit_object, par = which_min_tail_ess)
  
pars = c(names( which_min_bulk_ess),names( which_min_tail_ess))
par(mfrow = c(2,dim(as.array(fit_object))[2]),oma = c(0,0,3,0),xpd = NA)
for(p in 1:length(pars)){
  ranked_sims = r_scale(as.array(fit_object)[,, pars[p]])
  tmp = hist(as.numeric(unlist(ranked_sims)),plot =F,breaks = breaks)
  ylim = c(0,max(tmp$counts))
  for(c in 1:dim(ranked_sims)[2]){
    hist(ranked_sims[,c],breaks= 20,main = str_to_title(paste('chain',c)) ,col = p+1,ylim = ylim)
  }}
mtext(text = paste("mixing diagnostics for chains of lowest 'bulk ESS' param. (",names(which_min_bulk_ess),", top) and 'tail ESS' param. (",names(which_min_tail_ess),", bottom)"),
      cex = 0.8,outer = T,font = 2,side = 3)
}
# # # Plot overall Rhat diagnostics
plotConvergence_Rhat <- function(fit_object){
samp <- as.array(fit_object) 
res <- monitor_extra(samp)

par(mfrow = c(2,2))
plot(res$Rhat,ylim = c(.99 ,1.1),ylab = "Rhat",main = str_to_title("Rhat"),xlab = 'estimated parameters')
abline(h = 1.01,col= 'red',lty = 2)
plot(res$sRhat,ylim = c(.99 ,1.1),ylab = "sRhat",main = str_to_title("classic split-Rhat"),xlab = 'estimated parameters')
abline(h = 1.01,col= 'red',lty = 2)
plot(res$zsRhat,ylim = c(.99 ,1.1),ylab = "zsRhat",main = str_to_title("rank normalised split-Rhat"),xlab = 'estimated parameters')
abline(h = 1.01,col= 'red',lty = 2)
plot(res$zfsRhat,ylim = c(.99 ,1.1),ylab = "zfsRhat",main = str_to_title("folded rank normalised split-Rhat"),xlab = 'estimated parameters')
abline(h = 1.01,col= 'red',lty = 2)
}
# # # Plot overall ESS diagnostics
plotConvergence_ESS <- function(fit_object){
  samp <- as.array(fit_object) 
  res <- monitor_extra(samp)
par(mfrow = c(3,2))
plot(res$seff,ylim = c(0,max(res$seff,na.rm=T)),ylab = "seff",main = str_to_title("classic ESS"), xlab = 'estimated parameters')
abline(h = 400,col= 'red',lty = 2)
plot(res$zseff,ylim = c(0,max(res$zseff,na.rm=T)),ylab = "zseff",main = str_to_title("bulk ESS"),xlab = 'estimated parameters')
abline(h = 400,col= 'red',lty = 2)
plot(res$tailseff,ylim = c(0,max(res$tailseff,na.rm=T)),ylab = "tailseff",main = str_to_title("tail ESS"),xlab = 'estimated parameters')
abline(h = 400,col= 'red',lty = 2)
plot(res$medsseff,ylim = c(0,max(res$medsseff,na.rm=T)),ylab = "medsseff",main = str_to_title("median ESS"),xlab = 'estimated parameters')
abline(h = 400,col= 'red',lty = 2)
plot(res$madsseff,ylim = c(0,max(res$madsseff,na.rm=T)),ylab = "madsseff",main = str_to_title("mad ESS"),xlab = 'estimated parameters')
abline(h = 400,col= 'red',lty = 2)
}
# # # Plot Posterior Predictive Check
plot_ppc <- function(fit_object,data.list ){
y_gen = extract(fit_object , pars = c("y_gen"))$y_gen
eval <- apply(y_gen,1,function(x){
  eval = confusionMatrix(data=as.factor(x), reference=as.factor(data.list$y),positive = "1")
  return(c(eval$overall,eval$byClass))
  })
eval = apply(eval,1,function(x){quantile(x = x,prob = c(0.05,0.95))})


plot(x = 1:dim(y_gen)[2],
     y = colMeans(y_gen)[order(colMeans(y_gen))],
     pch = NA,
     ylim = c(-0.1,1.1),
     xlab = "observation order by predicted probability of label",
     ylab = 'label likelihood',
     main = str_to_title('posterior predictive check'))
points(x = 1:dim(y_gen)[2],
       y = colMeans(y_gen)[order(colMeans(y_gen))],
       col = adjustcolor('blue',0.25))
points(x = 1:dim(y_gen)[2],
       y = jitter(data.list$y[order(colMeans(y_gen))],factor = 0.05),
       col = adjustcolor(col = 'black',0.25),
       cex = 0.75)
lines(x = 1:dim(y_gen)[2],y = colMeans(y_gen)[order(colMeans(y_gen))],col =  adjustcolor(col = 'blue',0.95),lwd=1.5)
legend("topleft",
       legend = c("observed label","average predicted label","observed smooth fit (mgcv)",
                  paste(
                  c("accuracy:", # a measure that agreement between prediction and observation
                    "precision:",# how many actual recruits in the people the model thinks are recruits ? 
                    "recall / sensitivity:",# how many people the model thinks are recruits, over the total number of recruits 
                    "specificity:", # how many people the model thinks are not recruits, over the total number of non-recruits
                    "F1-score:"),# a statistic of accuracy that accounts for class imbalance - harmonic mean privileges importance of false negatives and false positives over truee ones. 
                  paste("(",apply(round(eval[,c("Accuracy","Precision","Recall","Specificity","F1")],2),
                        2,function(x){paste0(x,collapse = ',')}),")",sep = ""),sep=" ")
                  ),
       pch = c(1,1,NA,rep(NA,5)),lty = c(NA,1,1,rep(NA,5)),col = c("black","blue","black",rep(NA,5)),cex = 0.8)

fit <- gam(formula = data.list$y[order(colMeans(y_gen))] ~ s(c(1:dim(y_gen)[2])))
lines(x = 1:dim(y_gen)[2],y = predict(fit)[1:dim(y_gen)[2]],col =  adjustcolor(col = 'black',0.65),lwd = 1.5)
abline(h = 0.5,lty = 2,col = adjustcolor('blue',0.5))
abline(v = mean(which(abs(colMeans(y_gen)[order(colMeans(y_gen))]-0.5)==min(abs(colMeans(y_gen)[order(colMeans(y_gen))]-0.5)))),
       lty = 2,col = adjustcolor('blue',0.5))
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # PARAMETER PLOTS   # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # fixed effects plot
plotbetas <- function(extractedbetas, 
                      mfrow, 
                      beta.names, 
                      beta.names.plot = NA, 
                      beta.names.plot.new, 
                      descale = NA,# # currently only supports de-scaling from 1sd standardisation 
                      mu.X, sd.X,
                      dich_vars_spot,
                      include_legend = TRUE) {
  
  if(!is.na(descale)){
  if(descale=="1sd"){ 
    extractedbetas [,1] =  extractedbetas [,1]- rowSums(extractedbetas [,-1] *  t(array(mu.X/sd.X,c(dim(extractedbetas)[2]-1,dim(extractedbetas)[1]))))
    extractedbetas[,-1] = extractedbetas[,-1]*t(array(1/sd.X,c(dim(extractedbetas)[2]-1,dim(extractedbetas)[1])))
  }else{warning("currently this function only supports de-standardisation for he standard 1sd scaling;\nyour coefficients have not been de-scaled;\n2sd de-scaling will be added in the future.")}}
    
  if(!any(is.na(beta.names.plot))){
    extractedbetas = extractedbetas[,match(beta.names.plot,beta.names)]
  }
  
  par(mfrow = mfrow)
  
  for (i in 1:dim(extractedbetas)[2]) {
    
    dense = density(extractedbetas[, i])
    where = c("topleft","topright")[which(abs(mean(extractedbetas[, i]) - range(dense$x))==max(abs(mean(extractedbetas[, i]) - range(dense$x))))]

    plot(
      dense,
      zero.line = FALSE,
      main = str_to_title(beta.names.plot.new[i]),
      col = NA,
      xlab = 'log-odds of recruitment'
    )
    d <- density(extractedbetas[, i])
    polygon(
      d,
      main = '',
      xaxt = "n",
      yaxt = "n",
      ylab = "",
      xlab = "",
      border = NA,
      col = adjustcolor("orange", 0.75)
    )
    abline(
      v = colMeans(extractedbetas)[i],
      lty = 2,
      col = adjustcolor("purple", 0.75),
      lwd = 2
    )
    abline(v = 0, lty = 1, lwd = 1)
    
    if (i == 1) {
      legend.var = c(
        paste("E[beta_", i, "] = ", round(mean(
          extractedbetas[, i]
        ), 3), sep = ""),
        paste("sd(beta_", i, ") = ", round(sd(
          extractedbetas[, i]
        ), 3), sep = ""),
        ifelse(
          mean(extractedbetas[, i]) > 0,
          paste("Pr(beta_", i, ">0) = ", round(
            sum(extractedbetas[, i] > 0) / length(extractedbetas[, i]),
            3
          ), sep = ""),
          paste("Pr(beta_", i, "<0) = ", round(
            sum(extractedbetas[, i] < 0) / length(extractedbetas[, i]),
            3
          ), sep = "")
        )
      )
    } else{
      legend.var = c(
        paste("E[beta_", i, "] = ", round(mean(
          extractedbetas[, i]
        ), 3), sep = ""),
        paste("sd(beta_", i, ") = ", round(sd(
          extractedbetas[, i]
        ), 3), sep = ""),
        ifelse(
          mean(extractedbetas[, i]) > 0,
          paste("Pr(beta_", i, ">0) = ", round(
            sum(extractedbetas[, i] > 0) / length(extractedbetas[, i]),
            3
          ), sep = ""),
          paste("Pr(beta_", i, "<0) = ", round(
            sum(extractedbetas[, i] < 0) / length(extractedbetas[, i]),
            3
          ), sep = "")
        )
      )
    }
    if(isTRUE(include_legend)){
      legend(
        where,
        legend = legend.var,
        cex = .75,
        bg = "transparent",
        box.lty=0
      ) }
  } }
# # # lambda plot
plotlambdas <- function(extractedlambdas) {
  plot(
    density(extractedlambdas),
    xlim = c(0, 1),
    zero.line = FALSE,
    main = "",
    col = NA,
    xlab = '% variance explained by ICAR model'
  )
  polygon(
    density(extractedlambdas),
    main = '',
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    border = NA,
    col = adjustcolor("orange", 0.75)
  )
  abline(
    v = mean(extractedlambdas),
    lty = 2,
    col = adjustcolor("purple", 0.75),
    lwd = 2
  )
  legend(
    ifelse(mean(extractedlambdas)>0.5,"topleft","topright") ,
    legend = paste("MCMC mean:", round(mean(
      extractedlambdas
    ), 3)),
    lty = 2,
    lwd = 2,
    col = 'purple'
  )
  
}
# # # Plot Total Effect 
plotPredicted_2way <- function(scale = "logit",
                               n = NA,
                               interaction.vars.table,
                               data.list,
                               interaction.name,
                               interaction.states.names, 
                               outcome.lab,
                               interaction.states.cols,
                               interaction.states.cols.dark,
                               extracted_beta,
                               pred.names,
                               y.max = 0.6,
                               separate.plots = F){
  
  interaction = interaction.vars.table[,1] * interaction.vars.table[,2]
  names(interaction) = interaction.name
  interaction.vars.table = cbind(interaction.vars.table,interaction)
  interaction.vars.table.center = data.list$mu_X[names(interaction.vars.table)]
  interaction.vars.table.scale = data.list$sd_X[names(interaction.vars.table)]
  
  if(is.na(data.list$scale_X)){
    interaction.vars.table.std = interaction.vars.table
  }
  if(data.list$scale_X=="1sd"){
    interaction.vars.table.std = (interaction.vars.table - t(array(interaction.vars.table.center, rev(dim(interaction.vars.table)))))/t(array(interaction.vars.table.scale, rev(dim(interaction.vars.table))))
  }
  if(data.list$scale_X=="2sd"){
    interaction.vars.table_dich = which(sapply(interaction.vars.table,function(x){all(sort(unique(x))==c(0,1))}))
    interaction.vars.table.tmp = (interaction.vars.table - t(array(interaction.vars.table.center,c(dim(interaction.vars.table)[2],dim(interaction.vars.table)[1]))))/t(array(2*interaction.vars.table.scale,c(dim(interaction.vars.table)[2],dim(interaction.vars.table)[1])))
    for(i in 1:length(interaction.vars.table_dich)){
      interaction.vars.table.tmp[[interaction.vars.table_dich[i]]] = interaction.vars.table[[interaction.vars.table_dich[i]]]
    }
    interaction.vars.table.std = interaction.vars.table.tmp
  }
  
  logit.total.effect = sapply(1:dim(extracted_beta)[1], function(i) {as.matrix(interaction.vars.table.std) %*% as.numeric(extracted_beta[i, match(colnames(interaction.vars.table), pred.names)])})
  # add oob intercept 
  intercept = t(array(extracted_beta[,1],c(dim(extracted_beta)[1],dim(logit.total.effect)[1])))
  logit.total.effect = logit.total.effect + intercept
  
  row.names(logit.total.effect) = interaction.states.names
  
  if(scale == "logit"){
  pred.effect = logit.total.effect
  xlab = "log-odds of"
  }
  if(scale == "prob"){
  pred.effect = inv_logit(logit.total.effect)
  xlab = "probability of"
  }
  if(scale == "n"){
  pred.effect = n*inv_logit(logit.total.effect)
  xlab = "expected number of"
  }
  
  if(scale == "rel"){
    pred.effect = inv_logit(logit.total.effect)/inv_logit(intercept)
    xlab = "probability relative to average profile"
  }
  
  if(separate.plots){
    rows = ceiling(dim( pred.effect)[1]/2)
    par(mfrow = c(rows,dim( pred.effect)[1]/2))
  for(i in 1:dim(interaction.vars.table)[1]){
      d <- hist(pred.effect[i, ],plot = F,breaks = 100)

        plot(
        d,
        col = adjustcolor(interaction.states.cols[i], 0.75),
        main = str_to_title(rownames(pred.effect)[i]),
        xlab = paste(xlab,outcome.lab),
        border = NA)
    
      abline(
        v = mean(pred.effect[i, ]),
        lty = 1,
        col = interaction.states.cols.dark[i],
        lwd = 1.5
      )
      
      abline(
        v = median(pred.effect[i, ]),
        lty = 2,
        col = interaction.states.cols.dark[i],
        lwd = 1.5
      )
      
      abline(
        v = mlv(pred.effect[i, ], method = "meanshift"),
        lty = 3,
        col = interaction.states.cols.dark[i],
        lwd = 1.5
      )

 legend("topright",
        legend = c(paste("mean:",round( mean(pred.effect[i, ]),2)),
                   paste("median:",round( median(pred.effect[i, ]),2)),
                   paste('mode:',round(mlv(pred.effect[i, ], method = "meanshift"),2))
                   ),
        lty = c(1,2,3)
      )
    }
  }else{

plot(
  density(pred.effect[1, ]),
  xlim = c(
    min(pred.effect),
    max(pred.effect)),
  ylim = c(0,y.max),
  main = "",
  col = NA,
  xlab = paste(xlab,outcome.lab),
  zero.line = FALSE
)
for (i in 1:dim(pred.effect)[1]) {
  d <- density(pred.effect[i, ])
  polygon(
    d,
    main = '',
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    border = NA,
    col = adjustcolor(c(
      interaction.states.cols
    )[i], 0.75)
  )
}
abline(
  v = rowMeans(pred.effect),
  lty = 2,
  col = adjustcolor(interaction.states.cols, 0.5),
  lwd = 2
)
text(
  x = rowMeans(pred.effect),
  y = 0.75*y.max,
  labels = rownames(pred.effect),
  srt = 50,
  cex = 0.85,
  pos = 1,
  col = interaction.states.cols
) }
}

# Moran's I
plotMoran <- function(data.list, extractedy,sample.pct=NA) {

W_nb = data.list$nb.matrix[data.list$small_area_id,data.list$small_area_id]

# Calculate Moran's I of residual
I.res = sapply(1:dim(extractedy)[1], function(j) {
  
  if(!is.na(sample.pct)){
  sample_id <- sample(x = 1:dim(extractedy)[2],size = sample.pct*dim(extractedy)[2],replace = F)
  temp.y = data.list$y[sample_id]
  temp.yhat = t(extractedy)[sample_id,j]
  temp.W_nb = W_nb[sample_id,sample_id]
  }else{
  temp.y = data.list$y
  temp.yhat = t(extractedy)[,j]
  temp.W_nb =  W_nb
  }
  
  Moran.I.bootleg(x = (temp.y -  temp.yhat),
                  W = apply(temp.W_nb, 1, function(x) { if (sum(x) > 0) { return(x / sum(x)) } else{ return(x) } } ) ) 
  } )

# Calculate Moran's I of outcome 
I.Y = Moran.I.bootleg(x = data.list$y , 
                      W =  apply(W_nb, 1, function(x) { if (sum(x) > 0) { return(x / sum(x)) } else{ return(x) } } ) )

d <- density(I.res)

# plot 
plot(d,
  main = str_to_title('global auto-correlation'),
  xlab = 'Moran\'s I',zero.line = FALSE,col = NA,
  xlim = c(min(c(d$x,1.1*I.Y)),max(c(d$x,1.1*I.Y)))) 


polygon(
  d,
  main = '',
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  border = NA,
  col = adjustcolor("orange", 0.75)
)

abline(
  v = c(-1 / (data.list$n - 1), I.Y, mean(I.res)),
  lty = c(1, 2, 3),
  col = 'purple',lwd = 1.5
)

legend(
  'topleft',
  legend = c(
    paste('E[I(Y)]:', round(-1 / (data.list$n - 1), 3)),
    paste('mean I(Y):', round(I.Y, 3)),
    paste('mean I(res.):', round(mean(I.res), 3))
  ),
  lty = c(1, 2, 3),
  lwd = 1,
  cex = 0.75,
  bty = "n"
)

}
# plot total small-area effects
plotTotalBYM2 <- function(data.list,  extracted, oma, mar,cex.axis = 0.5) {

extractedgamma = extracted$gamma
  
par(oma = oma, mar = mar)

test_under = apply( extractedgamma, 2, function(x) { sum(x < 0) / dim( extractedgamma)[1] }) 
test_over = apply( extractedgamma, 2, function(x) { sum(x > 0) / dim( extractedgamma)[1] })


plot(
  1:dim( extractedgamma)[2],
  colMeans( extractedgamma)[order(test_over)],
  bty = "n",
  xlab = '',
  xaxt = "n",
  ylab = 'effect size (log-odds of recruitment)',
  main = "",
  cex = 1.2,
  lwd = 2,
  pch = NA,
  ylim = c(min(
    colMeans( extractedgamma)[order(test_over)] - 1.96 * apply( extractedgamma, 2, sd)[order(test_over)]
  ),
  max(
    colMeans( extractedgamma)[order(test_over)] + 1.96 * apply( extractedgamma, 2, sd)[order(test_over)]
  ))
)
abline(h = 0, col = 'black', lty = 2)



segments(
  x0 = 1:dim( extractedgamma)[2],
  y0 = colMeans( extractedgamma)[order(test_over)] - apply( extractedgamma, 2, sd)[order(test_over)],
  x1 = 1:dim( extractedgamma)[2],
  y1 = colMeans( extractedgamma)[order(test_over)] + apply( extractedgamma, 2, sd)[order(test_over)],
  lty = 1,
  lwd = 2,
  col = ifelse(
    test_under > 0.8,
    'skyblue',
    ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
  )[order(test_over)]
)
segments(
  x0 = 1:dim( extractedgamma)[2],
  y0 = colMeans( extractedgamma)[order(test_over)] - 1.96 * apply( extractedgamma, 2, sd)[order(test_over)],
  x1 = 1:dim( extractedgamma)[2],
  y1 = colMeans( extractedgamma)[order(test_over)] + 1.96 * apply( extractedgamma, 2, sd)[order(test_over)],
  lty = 1,
  lwd = 1,
  col = ifelse(
    test_under > 0.8,
    'skyblue',
    ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
  )[order(test_over)]
)


points(
  1:dim( extractedgamma)[2],
  colMeans( extractedgamma)[order(test_over)],
  pch = 0,
  lwd = 2,
  col = ifelse(
    test_under > 0.8,
    'skyblue',
    ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
  )[order(test_over)]
)

axis(
  side = 1,
  at = seq(1, dim( extractedgamma)[2], by = 2),
  labels = (data.list$small_area_names)[order(test_over)][seq(1, dim( extractedgamma)[2], by = 2)],
  las = 3 ,
  cex.axis = cex.axis
)
axis(
  side = 3,
  at = seq(2, dim( extractedgamma)[2], by = 2),
  labels = (paste(
    data.list$large_area_names, data.list$small_area_names, sep = ' - '
  )[order(test_over)])[seq(2, dim( extractedgamma)[2], by = 2)],
  las = 3 ,
  cex.axis = cex.axis
)
mtext(
  side = 3 ,
  outer = TRUE,
  text = 'small-area effect (unstructured + spatial)',
  xpd = NA,
  cex = 1.25,
  font = 2
)

legend(
  'topleft',
  legend = c("Pr(gamma<0) > 0.8", "Pr(gamma>0) > 0.8", 'zero'),
  pch = c(0, 0, NA),
  lty = c(1, 1, 2),
  col = c("skyblue", 'lightcoral', 'black'),
  lwd = 2,
  bg = "transparent",
  box.lty = 0
)


}
# plot total small-area effects
plotLargeArea <- function(data.list,  extracted, oma, mar,cex.axis = 1) {
  
  extractedeta <- extracted$eta*array(extracted$sigma_eta,c(dim( extracted$eta)[1],dim( extracted$eta)[2]))
  
  par(oma = oma, mar = mar)
  
  test_under = apply( extractedeta, 2, function(x) { sum(x < 0) / dim( extractedeta)[1] }) 
  test_over = apply( extractedeta, 2, function(x) { sum(x > 0) / dim( extractedeta)[1] })
  
  
  plot(
    1:dim( extractedeta)[2],
    colMeans( extractedeta)[order(test_over)],
    bty = "n",
    xlab = '',
    xaxt = "n",
    ylab = 'effect size (log-odds of recruitment)',
    main = "",
    cex = 1.2,
    lwd = 2,
    pch = NA,
    ylim = c(min(
      colMeans( extractedeta)[order(test_over)] - 1.96 * apply( extractedeta, 2, sd)[order(test_over)]
    ),
    max(
      colMeans( extractedeta)[order(test_over)] + 1.96 * apply( extractedeta, 2, sd)[order(test_over)]
    ))
  )
  abline(h = 0, col = 'black', lty = 2)
  
  
  
  segments(
    x0 = 1:dim( extractedeta)[2],
    y0 = colMeans( extractedeta)[order(test_over)] - apply( extractedeta, 2, sd)[order(test_over)],
    x1 = 1:dim( extractedeta)[2],
    y1 = colMeans( extractedeta)[order(test_over)] + apply( extractedeta, 2, sd)[order(test_over)],
    lty = 1,
    lwd = 2,
    col = ifelse(
      test_under > 0.8,
      'skyblue',
      ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
    )[order(test_over)]
  )
  segments(
    x0 = 1:dim( extractedeta)[2],
    y0 = colMeans( extractedeta)[order(test_over)] - 1.96 * apply( extractedeta, 2, sd)[order(test_over)],
    x1 = 1:dim( extractedeta)[2],
    y1 = colMeans( extractedeta)[order(test_over)] + 1.96 * apply( extractedeta, 2, sd)[order(test_over)],
    lty = 1,
    lwd = 1,
    col = ifelse(
      test_under > 0.8,
      'skyblue',
      ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
    )[order(test_over)]
  )
  
  
  points(
    1:dim( extractedeta)[2],
    colMeans( extractedeta)[order(test_over)],
    pch = 0,
    lwd = 2,
    col = ifelse(
      test_under > 0.8,
      'skyblue',
      ifelse(test_over > 0.8, 'lightcoral', 'darkgrey')
    )[order(test_over)]
  )
  
  axis(
    side = 1,
    at = seq(1, dim( extractedeta)[2], by = 2),
    labels = (paste(
      data.list$large_area_names, sep = ''
    )[order(test_over)])[seq(1, dim( extractedeta)[2], by = 2)],
    las = 3 ,
    cex.axis = cex.axis
  )
  axis(
    side = 3,
    at = seq(2, dim( extractedeta)[2], by = 2),
    labels = (paste(
      data.list$large_area_names, sep = ''
    )[order(test_over)])[seq(2, dim( extractedeta)[2], by = 2)],
    las = 3 ,
    cex.axis = cex.axis
  )
  mtext(
    side = 3 ,
    outer = TRUE,
    text = 'unstructured large-area effect',
    xpd = NA,
    cex = 1.25,
    font = 2
  )
  
  legend(
    'topleft',
    legend = c("Pr(eta<0) > 0.8", "Pr(eta>0) > 0.8", 'zero'),
    pch = c(0, 0, NA),
    lty = c(1, 1, 2),
    col = c("skyblue", 'lightcoral', 'black'),
    lwd = 2,
    bg = "transparent",
    box.lty = 0
  )
  
  
}
# plot spatial effects on a map
map.param <- function(data.list, extracted, param, area = "small", shape_area_id){

  if(area == "small") { shape = data.list$shape_small; area_id = data.list$small_area_names[data.list$small_area_id]; area.names = data.list$small_area_names}
  if(area == "large") { shape = data.list$shape_large; area_id = data.list$large_area_names[data.list$large_area_id]; area.names = data.list$large_area_names}
  
  if(param == "residual"){
  extracted$residual = t(sapply(1:dim(extracted$mu)[1],function(s){t(data.list$y) - extracted$y_gen[s,]}))
  }
  
  if(param == "eta"){
    extracted$eta <- extracted$eta*array(extracted$sigma_eta,c(dim( extracted$eta)[1],dim( extracted$eta)[2]))
  }
  if(param == "psi"){
    extracted$psi <- extracted$psi*array(extracted$sigma_gamma,c(dim( extracted$psi)[1],dim( extracted$psi)[2]))
  }
  if(param == "phi"){
    extracted$phi <- extracted$phi*array(extracted$sigma_gamma,c(dim( extracted$phi)[1],dim( extracted$phi)[2]))
  }

  # if a param is at the individual level aggregate to the relevant area
  if(param == "residual" | param == "y" | param == "mu"){
    # put y in the parm list
    if(param=="y"){extracted$y = t(data.list$y) 
    # aggregate
    outcome = data.table(t(extracted[param][[1]]),area_id = area_id)
    outcome = outcome[,lapply(.SD,sum),by = "area_id"]
    outcome.area_id = outcome$area_id
    outcome  = outcome [,!"area_id"]
    }else{
    # aggregate
    outcome = data.table(t(extracted[param][[1]]),area_id = area_id)
    outcome = outcome[,lapply(.SD,mean),by = "area_id"]
    outcome.area_id = outcome$area_id
    outcome = apply(outcome[,!"area_id"],1,mean)}
  }else{
    # if it's at the area-level, no need to aggregate
    outcome = apply(t(extracted[param][[1]]),1,mean)
    }

  # if at the individual level will still have that area_id 
  if(param == "residual" | param == "y" | param == "mu"){
  shape$outcome = 0
  if(param=="y"){
  shape$outcome[match(outcome.area_id,shape[shape_area_id][[1]])] = as.numeric(as.character(unlist(outcome)))
  }else{
  shape$outcome[match(outcome.area_id,shape[shape_area_id][[1]])] = as.numeric(as.character(unlist(outcome)))
  } }else{
    # otherwie, make sure the levels match and stick the parm estimate in the shapefile 
    if(all(shape[shape_area_id][[1]]==area.names)){
    shape$outcome  = outcome
    }else{stop("shapefile order does not match area-parameter order")}
    }
  

  if(max(abs(shape$outcome))<=1){
    if(min(shape$outcome)<0){
  pal.seq = seq(-1, 1, length.out = 11) 
  col.scale = "Spectral"
    }else{
  pal.seq = seq(0, 1, length.out = 9)
  col.scale = "YlOrRd"
    } } else{
  
      if(max(shape$outcome)>1 & min(shape$outcome)>=0){
  pal.seq = seq(0, max(shape$outcome), length.out = 9)
  col.scale = "YlOrRd"    
      }else{
  pal.seq = seq(-max(abs(apply(extracted[param][[1]],2,mean))), max(abs(apply(extracted[param][[1]],2,mean))), length.out = 11)
  col.scale = "Spectral"  
  } 
      
      }
  
  if(param == "y"){param = "outcome"}
  if(param == "gamma"){param = "BYM2 small-area effect"}
  if(param == "psi"){param = "ICAR small-area effect"}
  if(param == "phi"){param = "unstructured small-area effect"}
  if(param == "eta"){param = "unstructured large-area effect"}
  
  tmp.main = str_to_title(paste(param,"spatial distribution"))
  tmp.main = gsub("Icar","ICAR",tmp.main)
  tmp.main = gsub("Bym2","BYM2",tmp.main)
  
    plot(
    shape["outcome"],
    breaks = pal.seq,
    pal = brewer.pal(length( pal.seq), name = col.scale)[-1] ,
    main = tmp.main,
    border = adjustcolor('black', 0.1)
  )
}

# plot ranked profiles
profile_ranks_logit <- function(sims ,
                                sims.quantiles ,
                                xlim ,
                                ylim ,
                                outcome ,
                                labels){
  
  
  
  plot(x = 1:dim(sims)[1],
       y = sims.quantiles ["50%",][order(sims.quantiles  ["50%",])],
       ylim = ylim,
       pch= NA,
       xlim = xlim,
       xaxt = "n",
       xlab = paste(outcome,'rank'),
       ylab = paste("log-odds of",outcome),
       main = str_to_title(paste('ranked',outcome,'profiles')))
  abline(h = 0,lty = 2)
  abline(h = mean(sims.quantiles ["50%",]),lty = 3,col = 'purple',lwd = 2)
  axis(side = 1,at =seq(1,dim(sims)[1] ,length.out = 30),
       labels =round(rev(seq(1,dim(sims)[1] ,length.out = 30))),
       cex.axis = 0.75)
  segments(x0 = 1:dim(sims)[1],
           x1 = 1:dim(sims)[1],
           y0 = sims.quantiles ["5%",][order(sims.quantiles ["50%",])],
           y1 = sims.quantiles ["95%",][order(sims.quantiles ["50%",])],
           lty = 1,
           col = adjustcolor('black',0.5))
  points(x = 1:dim(sims)[1],
         y = sims.quantiles ["50%",][order(sims.quantiles ["50%",])],
         pch = 0,cex = 0.5)
  points(x = c(1:dim(sims)[1])[round(seq(1,dim(sims)[1],length.out = 10))],
         y = sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
         pch = 15,cex = 0.75,col = 15)
  abline(v = round(seq(1,dim(sims)[1],length.out = 10)),lty = 3,col = 'orange')
  segments(x0 = round(seq(1,dim(sims)[1],length.out = 10)),
           x1 = rep(c(10,-10),5) + round(seq(1,dim(sims)[1],length.out = 10)),
           y0 = sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           y1 = rep(c(-1,1),5)*7.5 + sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           lty = 1,lwd = 2,
           col = adjustcolor('orange',0.75))
  points(x = rep(c(10,-10),5) +round(seq(1,dim(sims)[1],length.out = 10)),
         y = rep(c(-1,1),5)*7.5 + sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
         pch = 1,
         cex = 0.65,col = adjustcolor('orange',0.75))
  for(i in 1:10){
    text(x = rep(c(10,-10),5)[i] +round(seq(1,dim(sims)[1],length.out = 10))[i],
         y = rep(c(-1,1),5)[i]*11.5 + sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))][i], 
         labels = labels[order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))][i],
         cex = 0.65)
  }
}
plot.ranks.quantiles = function(sims,
                                sims.quantiles,
                                labels ,
                                outcome ,unit,
                                ylim.point,ylim.min,ylim.max,
                                xlim.min,xlim.max){
  par(mfrow =c (2,2))
  x.ratio = 10/sum(abs(xlim.min ),xlim.max )
  # 90% interval
  plot(x = 1:dim(sims)[1],
       y = sims.quantiles ["50%",][order(sims.quantiles  ["50%",])],
       ylim = c(min(sims.quantiles),max(sims.quantiles)),
       pch= NA,
       xlim = c(xlim.min,xlim.max),
       xaxt = "n",
       xlab = paste(outcome,'rank'),
       ylab = paste(unit,"of",outcome),
       main = str_to_title(paste('ranked',outcome,'profiles\n90% interval')))
  points(x = c(1:dim(sims)[1]),
         y = sims.quantiles ["50%",][order(sims.quantiles ["50%",])],
         pch = 0,cex = 0.25)
  axis(side = 1,at =seq(1,dim(sims)[1] ,length.out = 30),
       labels =round(rev(seq(1,dim(sims)[1] ,length.out = 30))),
       cex.axis = 0.75)
  segments(x0 = 1:dim(sims)[1],
           x1 = 1:dim(sims)[1],
           y0 = sims.quantiles ["5%",][order(sims.quantiles ["50%",])],
           y1 = sims.quantiles ["95%",][order(sims.quantiles ["50%",])],
           lty = 1,
           col = adjustcolor('black',0.5))
  # point estimates
  plot(x = 1:dim(sims)[1],
       y = sims.quantiles ["50%",][order(sims.quantiles  ["50%",])],
       pch= 1,
       ylim = c(0,ylim.point),
       xlim = c(xlim.min,xlim.max),
       xaxt = "n",
       xlab = paste(outcome,'rank'),
       ylab = paste(unit,"of",outcome),
       main = str_to_title(paste('ranked',outcome,'profiles\npoint estimates')))
  points(x = 1:dim(sims)[1],
         y = sims.quantiles ["50%",][order(sims.quantiles  ["50%",])],
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  axis(side = 1,at =seq(1,dim(sims)[1] ,length.out = 30),
       labels =round(rev(seq(1,dim(sims)[1] ,length.out = 30))),
       cex.axis = 0.75)
  offset.x = -c(rep(c(x.ratio*sum(abs(xlim.min ),xlim.max ),
                      2*x.ratio*sum(abs(xlim.min ),xlim.max )),5)[-10],
                x.ratio*sum(abs(xlim.min ),xlim.max )) 
  offset.y = c(rep(c(ylim.point/10,ylim.point/5),5)[-10],ylim.point/10)
  segments(x0 = round(seq(1,dim(sims)[1],length.out = 10)),
           x1 =  1.35*offset.x + round(seq(1,dim(sims)[1],length.out = 10)),
           y0 = sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           y1 = 1.35*offset.y + sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           lty = 1,lwd = 2,
           col = adjustcolor('orange',0.75))
  points(x =  1.35*offset.x  + round(seq(1,dim(sims)[1],length.out = 10)),
         y = 1.35*offset.y  + sims.quantiles ["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  text(x = c(rep(c(2,1.5),5)[-10],2)*offset.x + round(seq(1,dim(sims)[1],length.out = 10)) ,
       y = 2*offset.y +sims.quantiles["50%",][order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
       labels = labels[order(sims.quantiles ["50%",])][round(seq(1,dim(sims)[1],length.out = 10))],
       cex = 0.75,srt  = -42.5)
  # 5th percentile
  plot(x = 1:dim(sims)[1],
       y = sims.quantiles ["5%",][order(sims.quantiles  ["5%",])],
       pch= 1,
       ylim = c(0,ylim.min),
       xlim = c(xlim.min,xlim.max),
       xaxt = "n",
       xlab = paste(outcome,'rank'),
       ylab = paste(unit,"of",outcome),
       main = str_to_title(paste('ranked',outcome,'profiles\nlowest plausible estimate (5th percentile)')))
  points(x = 1:dim(sims)[1],
         y = sims.quantiles ["5%",][order(sims.quantiles  ["5%",])],
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  axis(side = 1,at =seq(1,dim(sims)[1] ,length.out = 30),
       labels =round(rev(seq(1,dim(sims)[1] ,length.out = 30))),
       cex.axis = 0.75)
  
  offset.x = -c(rep(c(x.ratio*sum(abs(xlim.min ),xlim.max ),
                      2*x.ratio*sum(abs(xlim.min ),xlim.max )),5)[-10],
                x.ratio*sum(abs(xlim.min ),xlim.max )) 
  offset.y = c(rep(c(ylim.min/10,ylim.min/5),5)[-10],ylim.min/10)
  
  segments(x0 = round(seq(1,dim(sims)[1],length.out = 10)),
           x1 = 1.2*offset.x + round(seq(1,dim(sims)[1],length.out = 10)),
           y0 = sims.quantiles ["5%",][order(sims.quantiles ["5%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           y1 = 1.2*offset.y + sims.quantiles ["5%",][order(sims.quantiles ["5%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           lty = 1,lwd = 2,
           col = adjustcolor('orange',0.75))
  points(x = 1.2*offset.x + round(seq(1,dim(sims)[1],length.out = 10)),
         y = 1.2*offset.y + sims.quantiles ["5%",][order(sims.quantiles ["5%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  text(x = c(rep(c(1.9,1.4),5)[-10],1.9)*offset.x + round(seq(1,dim(sims)[1],length.out = 10)) ,
       y = 1.75*offset.y + sims.quantiles["5%",][order(sims.quantiles ["5%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
       labels = labels[order(sims.quantiles ["5%",])][round(seq(1,dim(sims)[1],length.out = 10))],
       cex = 0.75,srt  = -42.5)
  
  # 95th percentile
  plot(x = 1:dim(sims)[1],
       y = sims.quantiles ["95%",][order(sims.quantiles  ["95%",])],
       pch= 1,
       ylim = c(0,ylim.max),
       xlim = c(xlim.min,xlim.max),
       xaxt = "n",
       xlab = paste(outcome,'rank'),
       ylab = paste(unit,"of",outcome),
       main = str_to_title(paste('ranked',outcome,'profiles\nlargest plausible estimate (95th percentile)')))
  points(x = 1:dim(sims)[1],
         y = sims.quantiles ["95%",][order(sims.quantiles  ["95%",])],
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  axis(side = 1,at =seq(1,dim(sims)[1] ,length.out = 30),
       labels =round(rev(seq(1,dim(sims)[1] ,length.out = 30))),
       cex.axis = 0.75)
  
  offset.x = -c(rep(c(x.ratio*sum(abs(xlim.min ),xlim.max ),
                      2*x.ratio*sum(abs(xlim.min ),xlim.max )),5)[-10],
                x.ratio*sum(abs(xlim.min ),xlim.max )) 
  offset.y = c(rep(c(ylim.max/10,ylim.max/5),5)[-10],ylim.max/10)
  
  segments(x0 = round(seq(1,dim(sims)[1],length.out = 10)),
           x1 = 1.1*offset.x + round(seq(1,dim(sims)[1],length.out = 10)),
           y0 = sims.quantiles ["95%",][order(sims.quantiles ["95%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           y1 = 1.1*offset.y + sims.quantiles ["95%",][order(sims.quantiles ["95%",])][round(seq(1,dim(sims)[1],length.out = 10))],
           lty = 1,lwd = 2,
           col = adjustcolor('orange',0.75))
  points(x = 1.1*offset.x + round(seq(1,dim(sims)[1],length.out = 10)),
         y = 1.1*offset.y + sims.quantiles ["95%",][order(sims.quantiles ["95%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
         pch = 16,
         cex = 0.65,col = adjustcolor('orange',0.75))
  text(x = c(rep(c(1.75,1.2),5)[-10],1.75)*offset.x + round(seq(1,dim(sims)[1],length.out = 10)) ,
       y = 1.5*offset.y +sims.quantiles["95%",][order(sims.quantiles ["95%",])][round(seq(1,dim(sims)[1],length.out = 10))], 
       labels = labels[order(sims.quantiles ["95%",])][round(seq(1,dim(sims)[1],length.out = 10))],
       cex = 0.75,srt  = -42.5)
}








