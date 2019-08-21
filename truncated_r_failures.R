simulate_fixed_censored_data <- function(r, pf1, shape, scale)
{
  u_vector <- runif(r)
  
  n <- r/pf1
  t_c <- qweibull(pf1, shape, scale)
  F_tc <- 1 - exp( -(t_c/scale)^shape )
  
  truncated_data <- scale * ( -log( 1-u_vector*F_tc ) )^( 1/shape )
  
  data <- list(Number_of_Failures = r,
               Censor_Time = t_c,
               Failure_Times = truncated_data,
               Total_Number = n)
  
  return(data)
}

