
rm(list = ls())   # CLEAR VARIABLES
dev.off()         # REMOVE GRAPHICS

library(MASS)
library(caret)
library(ggplot2)
source("a3_q6_generate_the_data.R")
source("a3_q6_execute_lasso.R")
source("a3_q6_do_a_set.R")

set.seed(1977)
RHO = 0.7
SIGMA2 = 10
COEFFICIENTS = c(runif(10, 0.5, 1), rep(0, times=20))

gamma = 0.5 # -0.3 # the max lambda is 10^(gamma)
n_points = 149

arr = ABC(RHO_ = 0.1,
          SIGMA2_ = 1,
          COEFFICIENTS_ = COEFFICIENTS)

plot(log(1+arr[seq(1, length(arr), 32)]), abs(arr[seq(2, length(arr), 32)]), type="l", lwd=2, pch=19,
     col="#ffcc00aa",
     xlim = c(0, 1.0),
     ylim = c(0, 1.1),
     main=expression(paste('|', beta, '| vs  log(1+',lambda,') where ', rho, ' = 0.7 and ', sigma^2,' = 10')))


for (i in 3:32) {
  lines(log(1+arr[seq(1, length(arr), 32)]), abs(arr[seq(i, length(arr), 32)]), lwd=2, pch=19, col=(i*4))  
}

ABC = function(RHO_ = 0.1, SIGMA2_ = 1, COEFFICIENTS_=c(runif(20, 0.5, 1), rep(0, times=20))) {
  
  df = generate_the_data(RHO_, SIGMA2_, COEFFICIENTS_)
  index = createDataPartition(df$y, group=c(71), list=FALSE, times=1)
  N1_ = df[index,] # retain what appears in index
  N2_ = df[-index,] # retain only what does not appear in index
  
  arr = rep(0, times=32*n_points) # THINK OF IT LIKE UINT32 BUT DUMB
  count = 1
  lambda_arr = c(0, 10^seq(from=-3, to=gamma, length=(n_points-1)))
  for (i in lambda_arr) {
    
    ESTIMATE = DO_ONE_SET_AT_LAMBDA(TRAINING = N1_,
                                    N2 = N2_,
                                    COEFFICIENTS_ = COEFFICIENTS_,
                                    LAMBDA_ = i)  
      
      arr[count] = i
      count = count + 1
      for (j in 2:32) {
        arr[count] = ESTIMATE[j-1]
        count = count + 1
      }
    }

  return(arr)
}





DO_ONE_SET_AT_LAMBDA = function(TRAINING, N2, COEFFICIENTS_, LAMBDA_) {
  
  TRAINING.Y = TRAINING[,1]
  TRAINING.X1 = TRAINING[,2:11]
  TRAINING.X2 = TRAINING[,12:31]
  
  DF.GOOD = as.data.frame(cbind(TRAINING.Y, TRAINING.X1))
  DF.FULL = as.data.frame(cbind(TRAINING.Y, TRAINING.X1, TRAINING.X2))
  
  model.OLS.GOOD = lm(TRAINING.Y ~ . - 1, data=DF.GOOD) # the well-specified model
  model.OLS.FULL = lm(TRAINING.Y ~ . - 1, data=DF.FULL) # the one with everything
  
  model.LASSO = do_lasso_and_lasso_related_activities(lambda_vector_ = LAMBDA_, train_df_ = TRAINING)
  
  
  # PREDICTIONS ON N2
  
  predict.LASSO = predict(model.LASSO, newdata=N2)
  predict.OLS.GOOD = predict(model.OLS.GOOD, newdata=N2)
  predict.OLS.FULL = predict(model.OLS.FULL, newdata=N2)
  
  # CALCULATE MSE
  
  rmse.OLS.FULL = RMSE(predict.OLS.FULL, N2$y)
  MSE.OLS.FULL = (rmse.OLS.FULL)^2
  
  rmse.LASSO = RMSE(predict.LASSO, N2$y)
  MSE.LASSO = (rmse.LASSO)^2
  
  rmse.OLS.GOOD = RMSE(predict.OLS.GOOD, N2$y)
  MSE.OLS.GOOD = (rmse.OLS.GOOD)^2
  
  # CALCULATE THE SS(BIAS)
  
  coef.OLS.FULL = model.OLS.FULL$coefficients
  bias.OLS.FULL = COEFFICIENTS_ - coef.OLS.FULL
  SSB.OLS.FULL = sum(bias.OLS.FULL^2)
  
  coef.OLS.GOOD = c(model.OLS.GOOD$coefficients, rep(0, times=20))
  bias.OLS.GOOD = COEFFICIENTS_ - coef.OLS.GOOD
  SSB.OLS.GOOD = sum(bias.OLS.GOOD^2)
  
  # THE LASSO COEFFICIENTS : [1] is the intercept, remove it
  coef.LASSO = coef(model.LASSO$finalModel, model.LASSO$finalModel$lambdaOpt)[2:31]
  bias.LASSO = COEFFICIENTS_ - coef.LASSO
  SSB.LASSO = sum(bias.LASSO^2)
  
  return(coef(model.LASSO$finalModel, model.LASSO$finalModel$lambdaOpt))
}
