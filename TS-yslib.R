library(stats)
library(itsmr)
library(magrittr)
#######################
    #out_error_plot
#######################
#out-of-sample forecasting error function & plot
out_error_plot = function(data, m, a_order, s_order = c(0L, 0L, 0L), period = NA){
  #data : time-series data
  #m : initial training set index = 1:m
  n = length(data)
  xhat = numeric(n-m)
  err = numeric(n -m)
  for (i in 0:(n-m-1)){
    fit = arima(data[1:(m+i)], order = a_order, seasonal = list(order = s_order, period = period),
                include.mean = F)
    xhat[i+1] = forecast::forecast(fit, h = 1)$mean
    err[i+1] = (data[m+i+1] - xhat[i+1])^2
  }
  plot.ts(data, ylim = c(min(min(data), min(xhat)), max(max(data), max(xhat))))
  lines(as.vector(time(data))[(m + 1):n], xhat, col = "red")
  return(mean(err))
}

#######################
      #mspe_plot
#######################
mspe_plot = function(data, m, a_order, s_order = c(0L, 0L, 0L), period = NA){
  n = length(data)
  fit = arima(data[1:m], order = a_order, seasonal = list(order = s_order, period =period),
              include.mean = F)
  xhat = forecast::forecast(fit, h = n-m)$mean
  mspe = (1 - mean(abs((xhat - data[(m+1):n])/data[(m+1):n])))*100
  plot.ts(data, ylim = c(min(min(data), min(xhat)), max(max(data), max(xhat))))
  lines(as.vector(time(data))[(m+1):n], xhat, col = "red")
  return(mspe)
}


#######################
      #test_coef
#######################
#test if all of coef away from zero
test_coef = function(model, p_val = 0.05){
  test = 2 * (1 - pnorm(abs(model$coef / sqrt(diag(model$var.coef))))) < p_val
  return(test)
}

#######################
    #compare_stats
#######################
#compare aic, bic, sigma, loglik
compare_stats = function(data, p = 0, d = 0, q = 0, P = 0, D = 0, Q = 0, period = NA){
  idx = matrix(0, ncol = 6)
  for (i in 0:P){
    for (j in 0:D){
      for (k in 0:Q){
        for (l in 0:p){
          for (m in 0:d){
            for (n in 0:q){
              idx = rbind(idx, c(l, m, n, i, j, k))
            }
          }
        }
      }
    }
  }
  idx = idx[c(-1, -2), ]
  out = data.frame()
  id = data.frame()
  if (!is.null(nrow(idx))){
    for (i in 1:nrow(idx)){
      possibleError <- tryCatch(
        if (sum(idx[i, 4:6]) == 0){
          fit = arima(data, order = idx[i, 1:3], include.mean = F)
          out = rbind(out, c(stringr::str_flatten(as.character(idx[i,1:3])), 
                             AIC(fit), BIC(fit), fit$sigma2, fit$loglik))
        } else{
          fit = arima(data, order = idx[i, 1:3], include.mean = F,
                      seasonal = list(order = idx[i, 4:6], period = period))
          out = rbind(out, c(stringr::str_flatten(as.character(idx[i,])), 
                             AIC(fit), BIC(fit), fit$sigma2, fit$loglik))
        },
        error=function(e) e
      )
      if(inherits(possibleError, "error")) next
    }
    colnames(out) = c("order" ,"AIC", "BIC", "sigma2", "loglik")
    return(out %>% 
             dplyr::mutate(AIC = as.double(AIC), BIC = as.double(BIC),
                           sigma2 = as.double(sigma2), loglik = as.double(loglik)) %>% 
             dplyr::arrange(AIC, BIC, sigma2, desc(loglik)))
  } else{
    if (sum(idx[4:6]) ==0){
      fit = arima(data, order = idx[1:3], include.mean = F)
      return(c(stringr::str_flatten(as.character(idx)), 
               AIC(fit), BIC(fit), fit$sigma2, fit$loglik))
    }else{
      fit = arima(data, order = idx[1:3], include.mean = F,
                  seasonal = list(order = idx[4:6], period = period))
      out = rbind(out, c(stringr::str_flatten(as.character(idx)), 
                         AIC(fit), BIC(fit), fit$sigma2, fit$loglik))
    }
  }
}
