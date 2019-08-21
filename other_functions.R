conditional_prob <- function(beta, eta, t_c, t_w){
  return( (pweibull(t_w, beta, eta) - pweibull(t_c, beta, eta))/
            (1-pweibull(t_c, beta, eta)) )
}

compute_coverage_prob <- function(quan, num_at_risk, p_true){
  p4 <- pbinom(quan[4], num_at_risk, p_true)
  p3 <- pbinom(quan[3], num_at_risk, p_true)
  
  p2 <- dbinom(quan[2], num_at_risk, p_true)+pbinom(quan[2], num_at_risk, p_true, lower.tail = F)
  p1 <- dbinom(quan[1], num_at_risk, p_true)+pbinom(quan[1], num_at_risk, p_true, lower.tail = F)
  
  return(c(p1, p2, p3, p4))
}

compute_gpq_conditional_prob <- function(MLEs, BT_MLEs, t_c, t_w)
{
  MLE_mu <- log (MLEs[[2]])
  MLE_sigma <- 1/MLEs[[1]]
  
  BT_mu <- log (BT_MLEs[[2]])
  BT_sigma <- 1/BT_MLEs[[1]]
  
  GPQ_Beta <- 1 / ( MLE_sigma*MLE_sigma / BT_sigma )
  GPQ_Eta <- exp( MLE_mu + (MLE_mu - BT_mu)/BT_sigma*MLE_sigma )
  
  p_star_star <- conditional_prob(GPQ_Beta, GPQ_Eta, t_c, t_w)
  
  return (p_star_star)
}