# from https://github.com/haozhestat/RFIntervals
# modified so doesn't retrain random forest
#' @export
RFOOBInterval <- function(rf,
                          x0,
                          alpha = 0.05,
                          symmetry = TRUE){

  ntest = nrow(x0)
    
  test_pred <- predict(rf, x0)
  oob_abs_error = sort(abs(rf$y - rf$predicted))
  oob_error = sort(rf$y - rf$predicted)
  
  upper_pred = rep(NA, ntest)
  lower_pred = rep(NA, ntest)
  
  ## symmetry = TRUE leads to the symmetric OOB Intervals
  ## symmetry = FALSE leads to the standard OOB Intervals
  if(symmetry){
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_abs_error,1-alpha)
      lower_pred[i] = test_pred[i] - quantile(oob_abs_error,1-alpha)
    }
  }
  else{
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_error, 1-alpha/2)
      lower_pred[i] = test_pred[i] + quantile(oob_error, alpha/2)
    }
    
  }
  
  return(list(pred = matrix(test_pred,ntest,1),
              lo = matrix(lower_pred,ntest,1),
              up = matrix(upper_pred,ntest,1),
              fit = matrix(rf$predicted,length(rf$predicted),1)))
           #   fit = matrix(predict(rf,x),n,1)))
}
