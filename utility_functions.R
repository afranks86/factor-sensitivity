get_bounds <- function(df, data_list, treatment_index=1) {
  beta <- df %>% filter(M==1, J==treatment_index) %>% pull(beta)
  B <- df %>% filter(J==treatment_index) %>% pull(B) %>% matrix(dat=., nrow=data_list$K, ncol=data_list$M, byrow=TRUE)
  nrms <- apply(B, 1, function(x) sum(x^2))
  tibble(K=1:data_list$K, norm=nrms, beta=beta)
}





cal_rv_est <- function(beta, norm_beta) {
  w <- (sd(tr)*beta)^2 / norm_beta
  sign(beta)*w/(1+w)
  
}

pmist <- function(y, prob, mu=10, sigma) {
    prob * VGAM::pfoldnorm(y, mu, sigma) + (1-prob) * VGAM::pfoldnorm(y, -1*mu, sigma)
}

qmist2n = function(q, mu=10, sigma=1, prob){
  # the min minmax below is computed to supply a range to the solver
  # the solution must be between the min and max
  # quantile of the mixed distributions

  minmax <- range(VGAM::qfoldnorm(q, mu, sigma), VGAM::qfoldnorm(q, -mu, sigma))
  uniroot(function(x) pmist(x, prob, mu, sigma) - q,
                interval = minmax,
                tol = 10^{-16})$root  

} #qmist2n


## get odds ratio distribution
r2_to_odds <- function(r2, probs=NULL, alpha=0.95) {

  prob <- mean(probs)
  sigma2_t <- prob*(1-prob)

  r2 <- abs(r2)
  mu_lambda <- 1/2 * 1/(sigma2_t) * r2/(1-r2)
  sigma_lambda <- sqrt(1/sigma2_t * r2/(1-r2))

  qalpha <- sapply(1:length(sigma_lambda), function(i) {
    if(is.nan(mu_lambda[i]))
      NaN
    else {
      # browser()
      #qmist2n(alpha, mu_lambda[i], sigma_lambda[i], prob)
      VGAM::qfoldnorm(alpha, mu_lambda[i], sigma_lambda[i])
    }
  })
  
  odds_ratio <- exp(qalpha)
  odds_ratio
}

compute_bounds_and_robustness <- function(df, R2_no_nc=0.0, R2_nc=R2_no_nc, NC_index, data_list=data_list, 
                             treatment_index=1, 
                             marginal_sensitivity_param=FALSE, probs=NULL, 
                             full_rank=FALSE,
                             r2y=NULL, sdy=1) {


  sigma_t_hat <- sd(unlist(data_list$x[, treatment_index]))

  m <- ifelse(full_rank, data_list$K, data_list$M)
  beta_naive <- df %>% filter(M==1, J==treatment_index) %>% pull(beta) * sdy
  
  if(full_rank) {

    B <- df %>% filter(J==treatment_index) %>% pull(B) %>% 
      matrix(dat=., nrow=data_list$K, ncol=data_list$M, byrow=TRUE)
    lambda <- df %>% 
      filter(J==treatment_index) %>% 
      filter(M==1) %>% 
      pull(lambda) %>% diag * sdy^2
    
      BB <- B %*% t(B)
      if(!is.null(r2y)){

        lambda_scaler <- ((diag(BB) + diag(lambda)) * r2y - diag(BB))/diag(lambda)
        lambda_scaler <- pmax(lambda_scaler, rep(0, length(lambda_scaler)))
        lambda_scaler[is.na(lambda_scaler)] <- 0
      }

      eig_gamma <- eigen(B %*% t(B) + lambda_scaler*lambda)
      Gamma_tilde <- eig_gamma$vectors %*% diag(sqrt(abs(eig_gamma$values)))  * sdy
      lambda <- lambda - lambda_scaler*lambda
    
  } else {
    Gamma_tilde <- df %>% filter(J==treatment_index) %>% pull(B) %>% 
      matrix(dat=., nrow=data_list$K, ncol=m, byrow=TRUE) * sdy
    lambda <- df %>% 
      filter(J==treatment_index) %>% 
      filter(M==1) %>% 
      pull(lambda) %>% diag * sdy^2
  }

  var_y <- diag(Gamma_tilde %*% t(Gamma_tilde) + lambda)
  
  Gamma_c <- matrix(Gamma_tilde[NC_index,], ncol = m)
  Gamma_c_inv <- MASS::ginv(Gamma_c)
  
  beta_naive_nc <- beta_naive[NC_index]
  nrms <- apply(Gamma_tilde, 1, function(x) sum(x^2))
  w <- (sigma_t_hat*beta_naive)^2 / nrms
  rv_init <- sign(beta_naive)*w/(1+w)
  
  # check compatibility of negative control assumptions #
  if (!all.equal(as.numeric(Gamma_c %*% Gamma_c_inv %*% beta_naive_nc),
                 as.numeric(beta_naive_nc), check.names = FALSE)) {
    stop('Negative control assumptions are not compatible.')
  }
  
  R2_min_nc <- sigma_t_hat^2*sum((Gamma_c_inv %*% beta_naive_nc)^2) / 
    (1 + sigma_t_hat^2*sum((Gamma_c_inv %*% beta_naive_nc)^2))
  # beta_nc <- sigma_t_hat^2*sqrt(1-R2)*Gamma_c_inv %*% tau_naive_c

  bias_no_nc <- (sqrt(apply(Gamma_tilde, 1, function(x) sum(x^2))) / sigma_t_hat) * sqrt(R2_no_nc/(1-R2_no_nc))
  bias_nc <- c(Gamma_tilde %*% Gamma_c_inv %*% beta_naive_nc)
  r2_less_than_min <- FALSE

  if (any(R2_nc < R2_min_nc)) {
    R2_nc <- R2_min_nc + 1e-9
    r2_less_than_min <- TRUE
  }
  
  
  if(ncol(MASS::Null(t(Gamma_c))) == 0){
    cat("PATEs are identifiable.")
    est_df <- cbind(beta_naive, 
                    beta_naive - bias_no_nc, beta_naive + bias_no_nc, 
                    beta_naive - bias_nc)
    colnames(est_df) <- c("R2_0", "R2_lwr", "R2_upr", "R2_nc_lwr", "R2_nc_upr")
    est_df$rv <- 1
    est_df$K <- 1:data_list$K
    
  } else {
    bias_ncperp <- apply((diag(m) - Gamma_c_inv %*% Gamma_c) %*% t(Gamma_tilde), 
                         2, function(x) sqrt(sum(x^2))) %o%
      (rep(sqrt(R2_nc / (sigma_t_hat^2*(1-R2_nc)) - sum((Gamma_c_inv %*% beta_naive_nc)^2)), each=2) *
         rep(c(-1, 1), length(R2_nc)))
    if(any(is.nan(bias_ncperp)))
      browser()
    
    bias_ncperp[abs(bias_ncperp) < 1e-12] <- 0

    est_df <- cbind(beta_naive, beta_naive - bias_no_nc, beta_naive + bias_no_nc, 
                    beta_naive - bias_nc, beta_naive - bias_nc + bias_ncperp)
    colnames(est_df) <- c("Naive", "R2_lwr", "R2_upr", "R2_nc_0", "R2_nc_lwr", "R2_nc_upr")

    est_df <- as_tibble(est_df)
    debiased <- beta_naive - bias_nc
    debiased[abs(debiased) < 1e-12] <- 0
    rv_mediate <- (debiased^2/
                     apply((diag(m) - Gamma_c_inv %*% Gamma_c) %*% t(Gamma_tilde), 2, function(x) sum(x^2)) + 
                     sum((Gamma_c_inv %*% beta_naive_nc)^2)) * sigma_t_hat^2
    
    rv <- rv_mediate / (1 + rv_mediate)
    est_df$rv_init <- rv_init
    est_df$rv_nc <- sign(debiased+1e-16)*rv
    est_df$rv_nc_cond <- (est_df$rv_nc - rv_init)/(1-rv_init)
    
    est_df$K <- 1:data_list$K
    est_df$r2_less_than_min <- r2_less_than_min

    ## r2y stats
    est_df$r2y_gc <- diag(Gamma_tilde %*% t(Gamma_tilde)) / (diag(Gamma_tilde %*% t(Gamma_tilde)) + diag(lambda))
    est_df$lambda <- diag(lambda)
    est_df$total_var <- diag(Gamma_tilde %*% t(Gamma_tilde)) + diag(lambda)
    
    ## C&H RV and XRV
    b <- beta_naive^2 * (sigma_t_hat^2/var_y)
    c <- -b
    a <- 1
    rv_ch <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    xrv <- b / (1+b)
    
    est_df$rv_ch <- rv_ch
    est_df$xrv <- xrv
    
    if(marginal_sensitivity_param) {
      if(is.null(probs))
        stop("Must specify probs when marginal_sensitivity_param is TRUE.")

      
      est_df$odds_init <- r2_to_odds(est_df$rv_init, probs=probs)
      est_df$odds_nc <- r2_to_odds(est_df$rv_nc, probs=probs)
      est_df$odds_ch <- r2_to_odds(est_df$rv_ch, probs=probs)
      est_df$odds_xrv <- r2_to_odds(est_df$xrv, probs=probs)
    }

    est_df
  }
 

}

compute_r2y <- function(stan_results) {
  
  stan_results %>% spread_draws(B[K, M], lambda[K]) %>% 
    nest(data = c(B, lambda, K, M, .chain)) %>% 
    { map(.$data, function(x) {
      lambda <- x %>% filter(M==1) %>% pull(lambda)
      B <- matrix(x$B, nrow=max(x$K), ncol=max(x$M), byrow=TRUE)
      diag(B %*% t(B)) / (diag(B %*% (t(B))) + lambda)
    }) } -> R2list
  
  R2list    
}

get_beta_w_norm <- function(stan_results, data_list, NC_index = 1,
                treatment_index=1, marginal_sensitivity_param = TRUE,
                probs=NULL, full_rank=FALSE, R2_no_nc=0.0, R2_nc=R2_no_nc, r2y=NULL, sdy=1) {

  stan_results %>% spread_draws(B[K, M], lambda[K], beta[K, J]) %>% 
    group_by(.draw) %>% 
    nest() -> nest_df
  
  beta_w_norm <- nest_df %>% 
    mutate(data = purrr::map(data, function(
      d) compute_bounds_and_robustness(d, R2_no_nc=R2_no_nc, R2_nc=R2_nc, NC_index=NC_index, 
                          treatment_index = treatment_index, data_list=data_list, 
                          marginal_sensitivity_param = marginal_sensitivity_param, 
                          probs=probs, full_rank=full_rank, r2y=r2y, sdy=sdy))) %>%
    unnest(cols = c(data))
  
}

compute_effect_intervals <- function(beta_w_norm) {
  beta_w_norm %>% 
    group_by(K) %>% summarize(lower_naive=quantile(Naive, 0.025, na.rm=TRUE),
                              upper_naive=quantile(Naive, 0.975, na.rm=TRUE),
                              lower_init=quantile(R2_lwr, 0.025, na.rm=TRUE),
                              upper_init=quantile(R2_upr, 0.975, na.rm=TRUE),
                              lower_nc=quantile(R2_nc_lwr, 0.025, na.rm=TRUE),
                              upper_nc=quantile(R2_nc_upr, 0.975, na.rm=TRUE),
                              lower_rv_init = quantile(abs(rv_init), 0.025, na.rm=TRUE),
                              lower_rv_nc = quantile(abs(rv_nc), 0.025, na.rm=TRUE),
                              lower_odds_init = quantile(abs(odds_init), 0.025, na.rm=TRUE),
                              lower_odds_nc = quantile(abs(odds_nc), 0.025, na.rm=TRUE),
                              lower_odds_ch = quantile(abs(odds_ch), 0.025, na.rm=TRUE),
                              lower_odds_xrv = quantile(abs(odds_xrv), 0.025, na.rm=TRUE),
    ) %>% 
    mutate(Metabolite = colnames(Y)) %>%
    mutate(Metabolite = fct_relevel(Metabolite, 
                                    "HDL", "Methylmercury", "Lead", "Iron", "Potassium", "LDL", 
                                    "Cadmium", "Sodium", "Triglycerides", "Glucose")) %>% 
    arrange(Metabolite)
}

make_interval_plot <- function(effect_intervals, nc_plot_index=1, title="",
                               labels=c(expression(R[paste(T,'~',U,'|',X)]^2~"=0"), 
                                        expression(R[paste(T,'~',U, '|', X)]^2~"="~R[paste("min")]^2))) {
  orig_color = "#023FA5"
  nc_color = "#CC1C2F"
  effect_intervals %>%
    ggplot() + 
    geom_segment(aes(x=forcats::fct_inorder(Metabolite), xend=forcats::fct_inorder(Metabolite), 
                     y=lower_nc, yend=upper_nc, col='null'), 
                 size=3, alpha=0.5) +
    geom_segment(aes(x=forcats::fct_inorder(Metabolite), xend=forcats::fct_inorder(Metabolite), 
                     y=lower_naive, yend=upper_naive, col='nuc'), size=1.5, alpha=0.75) +
    geom_rect(aes(xmin=as.numeric(forcats::fct_inorder(Metabolite))-0.1, xmax=as.numeric(forcats::fct_inorder(Metabolite))+0.1, 
                     ymin=lower_init, ymax=upper_init, col="rect"), fill=alpha("white", 0)) +
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(2, 10),
              label=ifelse(sign(lower_naive) == sign(upper_naive), format(round(lower_odds_ch, digits=1), nsmall=1), "X")), fontface="bold") +    
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1.8, 10),
              label=ifelse(sign(lower_naive) == sign(upper_naive), format(round(lower_odds_xrv, digits=1), nsmall=1), "X")), fontface="bold") +
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1.6, 10),
                  label=ifelse(sign(lower_naive) == sign(upper_naive), format(round(lower_odds_init, digits=1), nsmall=1), "X")), fontface="bold") +
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1.4, 10),
                  label=ifelse(sign(lower_nc) == sign(upper_nc), format(round(lower_odds_nc, digits=1), nsmall=1), "")), col=nc_color, fontface="bold", nudge_x=.3) +
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1.4, 10),
                  label=ifelse(sign(lower_nc) == sign(upper_nc), format(round(lower_odds_nc / lower_odds_init[nc_plot_index], digits=1), nsmall=1), "")), 
              col=nc_color, fontface="bold", nudge_x=-.25) +
    geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1.4, 10),
              label=ifelse(sign(lower_nc) == sign(upper_nc), "/", "X")), col=nc_color, fontface="bold") +
    
    geom_text(data=tibble(), aes(label=c(as.character(expression(RV["NC"]^Gamma)), as.character(expression(RV^Gamma)), "XRV", as.character(expression(RV^1))), 
                                 x=rep(-.2, 4), y=seq(1.4, 2, by=0.2)), 
              fontface="bold", hjust=0, col="dark blue", parse=TRUE) +
    geom_point(aes(x=nc_plot_index, y=rep(0, length(nc_plot_index))), size=3, colour=nc_color, alpha=0.5, data=tibble()) +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    xlab("") + ylab("Treatment Effect") +  expand_limits(x=-.35) +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_color_manual(name="", labels=labels, 
    #                    values=c(nuc="black", 'null'=nc_color)) + 
    ggtitle(title) + theme(legend.text.align = 0)
}

