#' Fit Waiting Time Distribution with saturated sampling
#'
#' @param data data frame (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model.
#' @param form an object of class "formula" (or one that can be coerced to that
#' class): a symbolic description of the model to be fitted. 
#' @param parameters model formulae for distribution parameters
#' @param start 
#' @param end 
#' @param reverse logical; Fit the reverse waiting time distribution.
#' @param id the name of the variable that identifies distinct individuals
#' @param subset 
#' @param robust logical; compute a robust estimate of variance.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of options, and is
#' na.fail if that is unset.
#' @param init starting values for the parameters.
#' @param ... 
#'
#' @returns satwtdttt returns an object of class "wtd" inheriting from "mle".
#' @export
#'
#' @examples
satwtdttt <- function(data, form, parameters=NULL, start=NA, end=NA, reverse=F, id=NA,
                      subset=NULL, robust=T, na.action=na.omit, init=NULL, ...) {
  
  browser()
  
  
  if(is.null(data) || (nrow(data)<1)) {
    stop("data must be non-empty")
  }
  
  if(!inherits(form, "formula") || attr(terms(form), "response")==0) {
    stop("obstime variable must be specified in model formula")
  }
  
  data <- as.data.table(data)
  
  if(!is.null(substitute(subset))) {
    
    rows <- enquo(subset)
    rows_val <- eval_tidy(rows, data)
    data <- data[rows_val ,]
    
    if(nrow(data)<1) {
      stop("data must be non-empty")
    }
    
  }
  
  
  obs.name <- all.vars(form)[1]
  covar.names <- unique(unlist(lapply(parameters, function(x) all.vars(x)[-1])))
  
  if(!(obs.name %in% names(data))) {
    stop(paste0("'", obs.name, "'", "is not in data"))
  }
  
  data <- na.action(data, cols = c(obs.name, covar.names))
  
  
  # creation of shifted dates
  
  # if(!is(data[[obs.name]], "Date") || !is(start, "Date") || !is(end, "Date"))
  #   stop(paste0("variables start, end and '", obs.name, "' must be all of class Date"))
  
  delta <- as.numeric(end - start)
  
  if(is.null(id) || length(id)!=1 || is.na(id)) {
    stop("The id variable must be provided")
  }
  
  if(!(id %in% names(data))) {
    stop(paste0("'", id, "'", "is not in data"))
  }
  
  
  # implementing saturated sampling
  
  # keep only unique dispensing times
  data <- data[, .SD[!duplicated(obs.name)], by = id]
  
  # 1) between first dispensation and last dispensation
  
  # add an artificial last observation to include index dates after the last observed dispensation until the end of the window
  data_max <- data[!duplicated(id)]
  data_max <- data_max[,obs_max := end]
  
  data_max <- data_max[,(obs.name):=NULL]
  setnames(data_max, "obs_max", obs.name)
  
  data <- rbind(data, data_max)
  
  # order by dispensing time
  data <- data[order(id,get(obs.name))]
  
  data_post <- data[get(obs.name) >= start & get(obs.name) <= end, .SD]
  
  # compute time from consecutive dispensations
  data_post <- data_post[, time_to_next := shift(get(obs.name), type = "lead") - get(obs.name) , by = id]
  data_post <- data_post[, time_to_next := as.numeric(time_to_next)]
  
  # remove the last row as it doesn't have time to the next one
  data_post <- data_post[!is.na(time_to_next)]
  
  # expand obs.name time_to_next times
  data_e_post <- data_post[rep(1:.N, times = time_to_next)]
  
  # compute time since last dispensation, considering all the possible index dates until the next dispensation
  bycols <- c(id, obs.name)
  data_e_post <- data_e_post[, days_since_last := (seq(.N)-1)+0.5, by = bycols]
  
  first_disp <- data[get(obs.name) >= start & get(obs.name) <= end, .SD[which.min(get(obs.name))], by = id]
  
  # 2) before first dispensation
  
  # data_pre_first_disp <- data %>% filter(get(obs.name) < first_disp & get(obs.name) >= first_disp-(end-start)) %>% slice_max(get(obs.name))
  
  data <- first_disp[, .(id, first_disp = disp_time)][data, on = "id"]
  data_pre_first_disp <- data[get(obs.name) < first_disp & get(obs.name) >= first_disp-(end-start), .SD[which.max(get(obs.name))], by = id]
  
  # add first dispensing (within start-end window)
  data_pre_first_disp <- rbind(data_pre_first_disp[,.SD, .SDcols = setdiff(names(data_pre_first_disp), "first_disp")], first_disp)
  
  # compute time from consecutive dispensations
  data_pre_first_disp <- data_pre_first_disp[, time_to_next := shift(get(obs.name), type = "lead") - get(obs.name) , by = id]
  data_pre_first_disp <- data_pre_first_disp[, time_to_next := as.numeric(time_to_next)]
  
  # remove the last row as it doesn't have time to the next one
  data_pre_first_disp <- data_pre_first_disp[!is.na(time_to_next)]
  
  # expand obs.name time_to_next times
  data_e_pre_first_disp <- data_pre_first_disp[rep(1:.N, times = time_to_next)]
  
  # compute time since last dispensation, considering all the possible index dates until the next dispensation
  bycols <- c(id, obs.name)
  data_e_pre_first_disp <- data_e_pre_first_disp[, days_since_last := (seq(.N)-1)+0.5, by = bycols]
  
  # max_disp_last_year <- as.numeric(data %>% filter(get(obs.name) < start & get(obs.name) >= start-365.25) %>% summarise(last_disp = max(get(obs.name))))
  
  max_disp_last_year <- data[get(obs.name) < first_disp & get(obs.name) >= first_disp-365.25, .SD[which.max(get(obs.name))], by = id]
  
  max_disp_last_year <- max_disp_last_year[, min_diff := start - disp_time]
  
  data_e_pre_first_disp <- max_disp_last_year[, .(id, min_diff)][data_e_pre_first_disp, on = "id"]
  
  data_e_pre_first_disp <- data_e_pre_first_disp[days_since_last > min_diff, .SD]
  
  # 3) bind time window before and after the first dispensation within start-end window ----
  data_e <- rbind(data_e_pre_first_disp[,.SD, .SDcols = setdiff(names(data_e_pre_first_disp), "min_diff")], data_e_post)
  
  ############
  
  # XXXX do formula rewrite like in wtdttt() ?
  disttmp <- attr(terms(form, specials=c("dlnorm", "dweib", "dexp")), "specials")
  
  dist <- if(isTRUE(disttmp$dlnorm==2)) "lnorm" # need isTRUE() as value can be NULL
  else if(isTRUE(disttmp$dweib==2)) "weib"
  else if(isTRUE(disttmp$dexp==2)) "exp"
  else stop("model must use one of dlnorm, dweib or dexp")
  
  if(dist == "lnorm") {
    
    newform <- days_since_last ~ dlnorm(logitp, mu, lnsigma)
    
  } else if (dist == "weib") {
    
    newform <- days_since_last ~ dweib(logitp, lnalpha, lnbeta)
    
  } else if (dist == "exp") {
    
    newform <- days_since_last ~ dexp(logitp, lnbeta)
    
  }
  
  if(nrow(data_e)==0)
    stop("All dates are out of the window defined by start and end")
  
  # start and end need to be different compared to those provided in the call of the function (which define the sampling window)
  # they also need to be numerical now, as days_since_last is numerical
  start <- 0
  end <- 365
  
  out <- wtdttt(form = newform, parameters = parameters,
                start = start, end = end, reverse = reverse, id = id,
                preprocess = F, init = init, data = data_e)
  
  
  if (!robust) {
    
    out <- out
    
  } else {
    
    vcov_s <- sand_vcov(out)
    out@vcov <- vcov_s
    
  }
  
  
  return(out)
  
  
}

