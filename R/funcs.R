# get logistic model from input
# dat_in is data frame of abundnace, site estimates by date
# resp is name of response variable, usually sg_prp
# new_val is number of obs for predicting from fitted mod
logis_est <- function(dat_in, resp = 'meanabu', new_vals = 100){

  pts <- dat_in 
  
  # logistic growth
  Asym <- max(pts[, resp], na.rm = T)
  xmid <- median(pts$Site, na.rm = T)
  scal <- quantile(pts$Site, 0.75, na.rm = T) - xmid
  form_in <- substitute(x ~ SSlogis(Site, Asym,  xmid, scal), 
                        list(x = as.name(resp)))
  
  # model
  mod <- try({nls(form_in, data = pts, na.action = na.exclude)}, silent = T)
  
  # values for prediction
  dep_rng <- range(pts[, 'Site'], na.rm = T)
  new.x <- seq(dep_rng[1], dep_rng[2], length = new_vals)
  
  # return NAs if model fail, else get model predictions
  if('try-error' %in% class(mod)) {
    pred <- rep(NA_real_, length = length(new.x))
    asym <- NA_real_
    logis_mod <- NA
    message('Cannot fit curve to seagrass data\n')
  } else {
    pred <- as.numeric(predict(mod, 
                               newdata = data.frame(Site = new.x)))
    pred <- data.frame(new.x, pred)
    names(pred) <- c('Site', resp)
    asym <- summary(mod)$coefficients['Asym', 'Estimate']
    logis_mod <- mod
  }
  
  # return output
  out <- list(pred = pred, asym = asym, logis_mod = logis_mod)
  return(out)
  
}


# get doc ests from fitted logistic regression curve
# dat_in is two column data frame first is depth, second is predicted proportion of points occupied by seagrass
# asym is asymptote estimate from logistic mod
get_ests <- function(dat_in, asym){
  
  if(class(dat_in) != 'data.frame') stop('Input must be data frame')
  
  # exp and resp columns from dat_in
  x_val <- dat_in[, 1]
  y_val <- dat_in[, 2]
  
  # first deriv
  inflect <- diff(y_val)/diff(x_val)
  ind_min <- which.min(inflect)
  
  est_fun <- NA
  z_cmin <- NA
  z_cmed <- NA
  z_cmax <- NA
  
  # get curve estimate if the minimum slope is not the last value
  if(ind_min != (nrow(dat_in) - 1)){
    
    inflect_val <- dat_in[ind_min + 1, ]
    slope_val <- inflect[ind_min]
    int_val <- inflect_val[, 2] - slope_val * inflect_val[, 1]
    est_fun <- function(x) slope_val * x + int_val
    z_cmax <- -1 * int_val / slope_val
    
    # get z_cmed, halfway between z_cmin and z_cmax
    # z_cmin is based on asymptote intercept with linear reg
    # z_cmin defaults to zero if value is extrapolated
    z_cmin <- max(c(0, (asym - int_val)/slope_val))
    z_cmed  <- z_cmin + ((z_cmax - z_cmin)/2)
    
  }
  
  # output
  out <- list(preds = dat_in, est_fun = est_fun, z_cmin = z_cmin, z_cmed = z_cmed, z_cmax = z_cmax)
  return(out)
  
}

# doc_in doc input object
# level percent level of prediction intervals
# nsim number of monte carlo simulations for each value in newdata
# trace logical for counter output
sens <- function(doc_in, level = 0.05, trace = T, remzero = T, ...){
  
  # get model from doc_in for vcov matrix
  mod <- attr(doc_in, 'logis_mod')
  if(!'nls' %in% class(mod)){
    message('No model in object')
    return(doc_in)
  }

  # get linear model from doc_in
  est_fun <- attr(doc_in, 'est_fun')
  if(!'function' %in% class(est_fun)){
    message('Inadequate logistic model')
    return(doc_in)
  }
  
  # predictor value
  Depth <- coefficients(mod)['xmid']
  names(Depth) <- 'Depth'
  
  if(trace) cat("\n")
  
  ### get lower and upper bounds on doc estimates
  # var for zcmin, zcmed, zcmax are all different
  # all depend on variance/covariance of beta/gamma 
  
  vcovmod <- vcov(mod)
  betavar <- vcovmod['xmid', 'xmid']
  gammavar <- vcovmod['scal', 'scal']
  covar <- vcovmod['xmid', 'scal']
  
  # quantile to eval given level
  zquant <- qnorm(1 - (level/2))
  
  # variance estimates for zcmed and zcmax change if beta - 2gamma (zcmin) is less than zero
  zc_min <- attr(doc_in, 'z_cmin')
  if(zc_min > 0){
    
    ## 
    # zcmin = beta - 2gamma
    zcmin_var <- betavar + 4 * gammavar - 4 * covar
    
    ## 
    # zcmed = beta
    zcmed_var <- betavar
    
  } else {
    
    ## 
    # zcmin = 0
    zcmin_var <- 0
    
    ## 
    # zcmed = (beta + 2gamma)/2
    zcmed_var <- (betavar + 4 * gammavar + 4 * covar) / 4
    
  }
  
  ## 
  # zcmax = beta + 2gamma
  # variance is same regardless of conditions above
  zcmax_var <- betavar + 4 * gammavar + 4 * covar
  
  # prediction intervals
  zcmin_print <- zquant * sqrt(zcmin_var)
  zcmed_print <- zquant * sqrt(zcmed_var)
  zcmax_print <- zquant * sqrt(zcmax_var)
  
  ##
  # lower prediction intervals
  z_cmin <- attr(doc_in, 'z_cmin') - zcmin_print
  z_cmed <- attr(doc_in, 'z_cmed') - zcmed_print
  z_cmax <- attr(doc_in, 'z_cmax') - zcmax_print
  lower_est <- list(
    lower_shift = -1 * c(zcmin_print, zcmed_print, zcmax_print), 
    z_cmin = z_cmin, z_cmed = z_cmed, z_cmax = z_cmax
  )
  
  # upper estimates based on uncertainty
  z_cmin <- attr(doc_in, 'z_cmin') + zcmin_print
  z_cmed <- attr(doc_in, 'z_cmed') + zcmed_print
  z_cmax <- attr(doc_in, 'z_cmax') + zcmax_print
  upper_est <- list(
    upper_shift = c(zcmin_print, zcmed_print, zcmax_print), 
    z_cmin = z_cmin, z_cmed = z_cmed, z_cmax = z_cmax
  )
  
  # replace all estimates with NA if z_cmax prediction interval includes zero
  if(remzero & lower_est$z_cmax <= 0){
    
    attr(doc_in, 'z_cmax') <- NA
    attr(doc_in, 'z_cmin') <- NA
    attr(doc_in, 'z_cmed') <- NA
    
    return(doc_in)
    
  }
  
  # output
  # fill lower_est, upper_est attributes
  attr(doc_in, 'lower_est') <- lower_est
  attr(doc_in, 'upper_est') <- upper_est
  return(doc_in)
  
}

#' Create a doc object
#' 
#' Wrapper for creating a doc object
#' 
#' @param  ls_in list input created internally within \code{\link{doc_est}}
#' 
#' @export doc
#' 
#' @return Returns a doc object to be used with S3 methods
#' 
#' @details 
#' This function is a simple wrapper to \code{\link[base]{structure}} that is used internally within other functions to create a doc object.  The function does not have to be used explicitly.    
#' 
doc <- function(ls_in){
  
  # sanity check
  if(any(!names(ls_in) %in% c('data', 'preds', 'logis_mod', 'est_fun', 'z_cmin', 'z_cmed', 'z_cmax', 'lower_est', 'upper_est'))) stop('Incorrect input for doc object')
  
  # create class, with multiple attributes
  structure(
    .Data = ls_in[['data']], 
    class = c('doc', 'data.frame'), 
    preds = ls_in$preds,
    logis_mod = ls_in$logis_mod, 
    est_fun = ls_in$est_fun,
    z_cmin = ls_in$z_cmin,
    z_cmed = ls_in$z_cmed, 
    z_cmax = ls_in$z_cmax,
    lower_est = NA,
    upper_est = NA
  )
  
}

# function for estimating depth of colonization
# also used for plots
# 'dat_in' is data from 'buff_ext'
# 'depth_var' is name of depth column in input data
# 'sg_var' is name of seagrass column in input data
# 'maxbin' maximum bin size
# 'minpts' minimum points in a bin
doc_est <- function(dat_in, depth_var = 'Site', sg_var = 'meanabu', maxbin = 0.5, minpts = 50){

  # sanity check
  if(!'data.frame' %in% class(dat_in)) dat_in <- data.frame(dat_in)
  
  # order by depth
  dat_in <- dat_in[order(dat_in[[depth_var]]), ]
  dat_in$depth <- dat_in[[depth_var]]

  # stop function if no seagrass found
  if(sum(dat_in[[sg_var]], na.rm = T) == 0) stop('No seagrass present')
  
  # create bins from shallow to inreasing depth
  # each is a maximum width, unless there are at least a minimum number of points in the bin that define the bin width
  dat_sub <- dat_in
  pts <- NULL
  while(1 <= nrow(dat_sub)){
    
    dep_rng <- with(dat_sub, Site <= (Site[1] + maxbin))
    dep_rng <- dat_sub[dep_rng, , drop = F]
    dep_rng <- dep_rng[1:min(minpts, nrow(dep_rng)), , drop = F]
    dep_bin <- with(dep_rng, data.frame(
      mean(Site, na.rm = T), 
      length(Site), 
      mean(meanabu, na.rm = T)
    )
    )
    
    pts <- rbind(pts, dep_bin)
    dat_sub <- dat_sub[-c(1:nrow(dep_rng)), , drop = F]
    
  }
  
  names(pts) <- c('Site', 'n', 'sg_prp')
  
  # remove Site bins with zero seagrass
  pts <- pts[pts$sg_prp > 0, ]
  
  ##
  # estimate a logistic growth function for the data
  
  mod_est <- logis_est(pts, 'sg_prp')
  
  # output
  preds <- mod_est$pred
  asym <- mod_est$asym
  logis_mod <- mod_est$logis_mod
  
  ##
  # calculate depth of col using get_est
  
  # if no curve fit
  if(any(is.na(logis_mod))){
    
    # create doc output
    ls_in <- list(data = pts, preds = preds, logis_mod = logis_mod, est_fun = NA, 
                  z_cmin = NA, z_cmed = NA, z_cmax = NA, lower_est = NA, upper_est = NA)
    
    out <- doc(ls_in)
    return(out)
    
  }
  
  # check if curve is monotonic descending
  if(!with(preds, all(sg_prp == cummin(sg_prp)))){
    
    ls_in <- list(data = pts, preds = preds, logis_mod = logis_mod, est_fun = NA, 
                  z_cmin = NA, z_cmed = NA, z_cmax = NA, lower_est = NA, upper_est = NA)
    
    out <- doc(ls_in)
    return(out)
    
  }
  
  # get doc estimates using get_ests functions
  ests <- get_ests(preds, asym)
  preds <- ests[['preds']]
  est_fun <- ests[['est_fun']]
  z_cmin <- ests[['z_cmin']]
  z_cmed <- ests[['z_cmed']]
  z_cmax <- ests[['z_cmax']]
  
  # all output
  ls_in <- list(data = pts, preds = preds, logis_mod = logis_mod, est_fun = est_fun, 
                z_cmin = z_cmin, z_cmed = z_cmed, z_cmax = z_cmax, lower_est = NA, 
                upper_est = NA)
  
  out <- doc(ls_in)
  return(out)
  
}

# get depth from site
# 'res' output from sense
# 'dt' original data
dep_est <- function(res, dt){
 
  # site max
  z_cmed <- attr(res, 'z_cmed')
  lo <- attr(res, 'lower_est')$'z_cmed'
  hi <- attr(res, 'upper_est')$'z_cmed'
  
  approx(x = dt$Site, y = dt$Depth_dem, xout = c(lo, z_cmed, hi))
  dep <- Hmisc::approxExtrap(x = dt$Site, y = dt$Depth_dem, xout = c(lo, z_cmed, hi))
  dep <- dep$y

  out <- tibble(
    var = c('Depth_dem', 'Site'), 
    lo = c(dep[1], lo), 
    z_cmed = c(dep[2], z_cmed), 
    hi = c(dep[3], hi)
  ) 
  
  return(out)
  
}
  
