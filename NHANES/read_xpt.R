library(tidyverse)
library(rstan)

xpt_files  <- list.files(pattern="*\\.xpt", ignore.case=TRUE)

df_list  <- lapply(xpt_files, function(x) haven::read_xpt(x))
names(df_list) <- xpt_files

df <- df_list %>% reduce(left_join, by="SEQN")

dim(df)

## URXUMS: Albumin, urine (mg/L)
## URXUCR: Creatinine, urine (mg/dL)
## URXUAS: Arsenic, Total urine (mg/L)
## LBXGLU: Fasting glucose (mg/DL)
## LBDHDD: Direct HDL-Chol (mg/DL)
## LBDLDLN LDL Chol NIN Eqn 2
## URXUCM: Chromium (urine, ug/L)


df %>% select(ends_with("SI"), URXUAS, URXUCM, URXUHG) %>%
  as.matrix -> mat

hclust(cr)$ord

df %>% select(ends_with("SI"), URXUAS, URXUCM, URXUHG) %>%
  cor(method="spearman", use="pairwise.complete.obs") -> cr

corrplot::corrplot(cr)

cr %>% eigen %>% .$values %>% plot


## w/ arsenic and chromium
df %>% select(ends_with("SI"), URXUAS, URXUCM, URXUHG) %>%

  colnames


colnames(df)


## Fish

## Alcohol
df_list[["P_ALQ.XPT"]] %>% summarize_all(function(x) mean(is.na(x)))


df_list[["P_ALQ.XPT"]] %>% pull(ALQ130) %>% table()
## ALQ130 - Avg # alcoholic drinks/day - past 12 mos

## Smoking
df_list[["P_SMQRTU.XPT"]]

## Background covariates
## RIDAGEYR - Age in years at screening - 0 to 24 mos
## RIAGENDR - Gender: 1 male, 2 female
## DMDEDUC2 - Education 1: less than 9th grade, 5 college or above
df %>% dplyr::select(RIDAGEYR, RIAGENDR, DMDEDUC2)

## ALQ121: 1 every day or 2 nearly every day a drink

## ``The one treated individual consumed between one and three drinks on most days. The controls drank once a week or less, and most did not drink at all.'' (Rosenbaum, 2021)
df %>% select(SEQN, ALQ121, ALQ130) %>% filter(!is.na(ALQ121) & !is.na(ALQ130)) %>%
  filter(!(ALQ121 %in% c(1, 2) & ALQ130 > 3)) %>%
  mutate(light_drinking = ALQ121 %in% c(1, 2)) ->
  drinking_df

mean(drinking_df$light_drinking)

drinking_SEQN <- drinking_df$SEQN
demo_df <- df %>% dplyr::select(SEQN, RIDAGEYR, RIAGENDR, DMDEDUC2) %>% filter(SEQN %in% drinking_SEQN)

## LBDBGMSI - Mercury, methyl (nmol/L)
## LBDHDDSI HDL cholesterol
#dplyr::select(SEQN, LBDBGMSI, LBDHDDSI, URXUAS, URXUCM, URXUHG)
outcome_df <- df %>% filter(SEQN %in% drinking_SEQN)



## HDL
summary(lm(outcome_df$LBDHDDSI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))
## Methylmercury
summary(lm(outcome_df$LBDBGMSI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))

## GLUCOSE
summary(lm(outcome_df$LBXSGL ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))

## Potasium
summary(lm(outcome_df$LBXSKSI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))

## Sodium
summary(lm(outcome_df$LBXSNASI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))



outcomes <- outcome_df %>% 
  dplyr::select(SEQN, LBDHDDSI, LBDBGMSI, LBXSGL, LBXSKSI, 
                LBXSNASI, LBDSIRSI, LBDBPBSI) %>%
  drop_na %>% arrange(SEQN)
predictors <- left_join(demo_df, drinking_df) %>%
  filter(SEQN %in% outcomes$SEQN) %>%
  arrange(SEQN) %>%
  dplyr::select(SEQN, light_drinking, RIDAGEYR, RIAGENDR)


X <- predictors %>%  mutate(light_drinking = ifelse(light_drinking, 1, 0),
                            gender = ifelse(RIAGENDR==2, 1, 0)) %>%
  select(-one_of("SEQN", "RIAGENDR")) %>% as.matrix

matched_predictors <- MatchIt::match.data(MatchIt::matchit(light_drinking ~ RIAGENDR + RIDAGEYR, data=predictors))
X <- matched_predictors %>% mutate(light_drinking = ifelse(light_drinking, 1, 0),
             gender = ifelse(RIAGENDR==2, 1, 0),
             age = RIDAGEYR) %>%
  arrange(SEQN) %>%
  select(light_drinking, gender, age) %>% as.matrix
Y <- outcomes %>% 
  filter(SEQN %in% matched_predictors$SEQN) %>% 
  arrange(SEQN) %>% 
  dplyr::select(-SEQN) %>% 
  mutate_all(log) %>%
  as.matrix
colnames(Y) <- c("HDL", "Methylmercury", "Glucose", "Potasium", "Sodium", "Iron", "Lead")

sm <- stan_model("../metabolomics/multivariate_regression.stan")


data_list <- list(K=ncol(Y), J=ncol(X), N=nrow(Y), M=2, x=X, y=scale(Y))

stan_results <- sampling(sm, data_list)

stan_results %>% spread_draws(beta[K, J]) %>% 
  group_by(J, K) %>% 
  summarize(q025=quantile(beta, 0.025), 
            q975=quantile(beta, 0.975)) %>%
  mutate(sig = sign(q025*q975)==1) %>% 
  mutate(feature = colnames(Y)) %>% 
  filter(J==1) %>% 
  arrange(desc(q025)) %>% relocate(feature, q025, q975)

get_bounds <- function(df, treatment_index=1) {
  beta <- df %>% filter(M==1, J==treatment_index) %>% pull(beta)
  B <- df %>% pull(B) %>% matrix(dat=., nrow=data_list$K, ncol=data_list$M, byrow=TRUE)
  nrms <- apply(B, 1, function(x) sum(x^2))
  tibble(K=1:data_list$K, norm=nrms, beta=beta)
}

stan_results %>% spread_draws(B[K, M], lambda[K], beta[K, J]) %>% 
  group_by(.draw) %>% 
  nest() -> nest_df

beta_w_norm <- nest_df %>% mutate(data = purrr::map(data, function(d) get_bounds(d))) %>%
  unnest(cols = c(data))

cal_rv_est <- function(beta, norm_beta) {
  w <- (sd(tr)*beta)^2 / norm_beta
  sign(beta)*w/(1+w)
  
}

R2 <- 0.0
tr <- X[, 1]
beta_w_norm %>% 
  group_by(K) %>% 
  mutate(beta_age_quantile=ifelse(abs(quantile(beta, 0.025)) < abs(quantile(beta, 0.975)), quantile(beta, 0.025), quantile(beta, 0.975))) %>% 
  ungroup() %>% 
  mutate(upper = beta + sqrt(R2/(1-R2)) * sqrt(norm)/sd(tr), 
         lower = beta - sqrt(R2/(1-R2)) * sqrt(norm)/sd(tr), 
         RV = cal_rv_est(beta_age_quantile, norm)) %>%
  group_by(K) %>% 
  summarize(mean_upper = quantile(upper, 0.975),                                    
            mean_lower = quantile(lower, 0.025),
            lower_rv = quantile(RV, 0.05),
            median_rv = median(RV),
            signif=sign(mean_lower*mean_upper)==1) %>%
  mutate(feature=colnames(Y)) %>% 
  filter(signif) %>% arrange(desc(signif), desc(lower_rv))

## Null controls

cal_tau_cali_wNC <- function(df, R2=0.5, NC_index, treatment_index=1) {
  
  sigma_t_hat <- sd(data_list$x[, treatment_index])
  R2_orig <- R2  
  m <- data_list$M
  beta_naive <- df %>% filter(M==1, J==treatment_index) %>% pull(beta)
  Gamma_tilde <- df %>% pull(B) %>% matrix(dat=., nrow=data_list$K, ncol=m, byrow=TRUE)
  
  
  Gamma_c <- matrix(Gamma_tilde[NC_index,], ncol = data_list$M)
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
  
  bias_no_nc <- sqrt(apply(Gamma_tilde, 1, function(x) sum(x^2)) / sigma_t_hat^2 * R2/(1-R2))
  bias_nc <- c(Gamma_tilde %*% Gamma_c_inv %*% beta_naive_nc)
  r2_less_than_min <- FALSE
  if (any(R2 < R2_min_nc)) {
    # warning(paste0("R^2 is less than ", 
    #             round(R2_min_nc, 2), ""))
    R2 <- R2_min_nc + 1e-9
    r2_less_than_min <- TRUE
  }
  
  if(ncol(MASS::Null(t(Gamma_c))) == 0){
    cat("PATEs are identifiable.")
    est_df <- cbind(beta_naive, beta_naive - bias_no_nc, beta_naive + bias_no_nc, beta_naive - bias_nc)
    browser()
    colnames(est_df) <- c("R2_0", "R2_lwr", "R2_upr", "R2_nc_lwr", "R2_nc_upr")
    est_df$rv <- 1
    est_df$K <- 1:data_list$K
    
  } else {
    bias_ncperp <- apply((diag(m) - Gamma_c_inv %*% Gamma_c) %*% t(Gamma_tilde), 
                         2, function(x) sqrt(sum(x^2))) %o%
      (rep(sqrt(R2 / (sigma_t_hat^2*(1-R2)) - sum((Gamma_c_inv %*% beta_naive_nc)^2)), each=2) *
         rep(c(-1, 1), length(R2)))
    if(any(is.nan(bias_ncperp)))
      browser()

    bias_ncperp[abs(bias_ncperp) < 1e-16] <- 0
    
    est_df <- cbind(beta_naive, beta_naive - bias_no_nc, beta_naive + bias_no_nc, beta_naive - bias_nc, beta_naive - bias_nc + bias_ncperp)
    
    colnames(est_df) <- c("Naive", "R2_lwr", "R2_upr", "R2_nc_0", "R2_nc_lwr", "R2_nc_upr")
    
    est_df <- as_tibble(est_df)
    debiased <- beta_naive - bias_nc
    debiased[abs(debiased) < 1e-16] <- 0
    rv_mediate <- (debiased^2/
                     apply((diag(m) - Gamma_c_inv %*% Gamma_c) %*% t(Gamma_tilde), 2, function(x) sum(x^2)) + 
                     sum((Gamma_c_inv %*% beta_naive_nc)^2)) * sigma_t_hat^2
    
    rv <- rv_mediate / (1 + rv_mediate)
    est_df$rv_init <- rv_init
    est_df$rv_nc <- sign(beta_naive)*rv
    est_df$K <- 1:data_list$K
    est_df$r2_less_than_min <- r2_less_than_min
    est_df
  }
}


stan_results %>% spread_draws(B[K, M], lambda[K], beta[K, J]) %>% 
  group_by(.draw) %>% 
  nest() -> nest_df

beta_w_norm <- nest_df %>% mutate(data = purrr::map(data, function(d) cal_tau_cali_wNC(d, R2=0.3, NC_index=2, treatment_index = 1))) %>%
  unnest(cols = c(data))

orig_color = "#023FA5"
nc_color = "#CC1C2F"
beta_w_norm %>% group_by(K) %>% summarize(lower_init=quantile(rv_init, 0.025, na.rm=TRUE),
                                          upper_init=quantile(rv_init, 0.975, na.rm=TRUE),
                                          lower_nc=quantile(rv_nc, 0.025, na.rm=TRUE),
                                          upper_nc=quantile(rv_nc, 0.975, na.rm=TRUE)) %>% 
  mutate(Feature = colnames(Y)) %>%
  filter(sign(lower_init) == sign(upper_init)) %>%
  mutate(rv_init = pmin(abs(lower_init), abs(upper_init)), rv_nc = pmin(abs(lower_nc), abs(upper_nc))) %>%
  arrange(desc(rv_init)) %>%
  select(-one_of("lower_init", "upper_init", "lower_nc", "upper_nc")) %>%
  pivot_longer(cols=c("rv_init", "rv_nc"), names_to="Type", values_to="rv") %>%
  ggplot() + geom_line(aes(x=forcats::fct_inorder(Feature), y=rv, col=Type, group=Type), size=1.4) + 
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  xlab("Feature") + ylab("Robustness Value") + 
  scale_color_manual(name="Type",
                     labels=c("No null controls", "Methylmercury is a null control"), values=c(orig_color, nc_color)) + 
  theme(legend.position = c(0.8, 0.8))




