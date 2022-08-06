library(tidyverse)
library(rstan)
library(tidybayes)
source("utility_functions.R")
options(mc.cores=10)

xpt_files  <- list.files(pattern="*\\.xpt", ignore.case=TRUE)

df_list  <- lapply(xpt_files, function(x) haven::read_xpt(x))
names(df_list) <- xpt_files

df <- df_list %>% reduce(left_join, by="SEQN")

## RIDAGEYR - Age in years at screening - 0 to 24 mos
## RIAGENDR - Gender: 1 male, 2 female
## DMDEDUC2 - Education 1: less than 9th grade, 5 college or above
df %>% dplyr::select(RIDAGEYR, RIAGENDR, DMDEDUC2)

## ALQ121: 1 every day or 2 nearly every day a drink

## ``The one treated individual consumed between one and three drinks on most days. The controls drank once a week or less, and most did not drink at all.'' (Rosenbaum, 2021)
df %>% select(SEQN, ALQ121, ALQ130) %>% filter(!is.na(ALQ121) & !is.na(ALQ130)) %>%
  filter(ALQ121 %in% c(1, 2) | ALQ121 >= 5) %>%
  filter(ALQ130 < 3) %>%
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

## potassium
summary(lm(outcome_df$LBXSKSI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))

## Sodium
summary(lm(outcome_df$LBXSNASI ~ drinking_df$light_drinking + demo_df$RIDAGEYR + demo_df$RIAGENDR))



outcomes <- outcome_df %>% 
  mutate(BPXODI = (BPXODI1 + BPXODI2 + BPXODI3)/3, BPXOSY = (BPXOSY1 + BPXOSY2 + BPXOSY3)/3) %>%
  dplyr::select(SEQN, LBDHDDSI, LBDBGMSI, LBXSGL, LBXSKSI, 
                LBXSNASI, LBDSIRSI, LBDBPBSI, LBDLDNSI, LBDTRSI, LBDBCDSI ) %>%
  drop_na %>% arrange(SEQN)
predictors <- left_join(demo_df, drinking_df) %>%
  filter(SEQN %in% outcomes$SEQN) %>%
  arrange(SEQN) %>%
  dplyr::select(SEQN, light_drinking, RIDAGEYR, RIAGENDR, DMDEDUC2) %>% 
  drop_na()


X <- predictors %>%  mutate(light_drinking = ifelse(light_drinking, 1, 0),
                            age = RIDAGEYR,
                            gender = ifelse(RIAGENDR==2, 1, 0),
                            educ = ifelse(DMDEDUC2 <=3, 0, 1)) %>%
  select(-one_of("SEQN", "RIDAGEYR", "RIAGENDR", "DMDEDUC2"))

## GLM calibration



preds_full <- predict(glm(light_drinking ~ RIDAGEYR + RIAGENDR + (DMDEDUC2 <= 3), data=predictors, family="binomial"), type="response")
preds_no_educ <- predict(glm(light_drinking ~ RIDAGEYR + RIAGENDR , data=predictors, family="binomial"), type="response")

preds_no_age <- predict(glm(light_drinking ~ RIAGENDR + (DMDEDUC2 <= 3), data=predictors, family="binomial"), type="response")
preds_no_gender <- predict(glm(light_drinking ~ RIDAGEYR + (DMDEDUC2 <= 3), data=predictors, family="binomial"), type="response")

or_educ <- ((preds_full/(1-preds_full)) / (preds_no_educ/(1-preds_no_educ)))
or_age <- ((preds_full/(1-preds_full)) / (preds_no_age/(1-preds_no_age)))
or_gender <- ((preds_full/(1-preds_full)) / (preds_no_gender/(1-preds_no_gender)))
or_educ <- ((preds_full/(1-preds_full)) / (preds_no_educ/(1-preds_no_educ)))

or_all <- ((preds_full/(1-preds_full)) / (mean(predictors$light_drinking)/(1-mean(predictors$light_drinking))))

tibble(age_or = or_age) %>% 
  ggplot() + geom_histogram(aes(x=age_or)) + 
  theme_bw(base_size=16) + ylab("Count") + xlab("Odds Ratio") +
  ggtitle("Change in odds ratio - Age")

tibble(gender_or = or_gender) %>% 
  ggplot() + geom_histogram(aes(x=gender_or)) + 
  theme_bw(base_size=16) + ylab("Count") + xlab("Odds Ratio") +
  ggtitle("Change in odds ratio - Gender")


tibble(both_or = pmax(1/or_both, or_both)) %>% 
  ggplot() + geom_histogram(aes(x=both_or)) + 
  theme_bw(base_size=16) + ylab("Count") + xlab("Odds Ratio") +
  ggtitle("Change in odds ratio - Gender") + xlim(c(1, 2))

max(1/min(or_age), max(or_age))
quantile(pmax(1/or_age, or_age), 0.95)

max(1/min(or_gender), max(or_gender))
quantile(pmax(1/or_gender, or_gender), 0.95)

max(1/min(or_educ), max(or_educ))
quantile(pmax(1/or_educ, or_educ), 0.95)


max(1/min(or_both), max(or_both))
quantile(pmax(1/or_both, or_both), 0.95)


summary(glm(light_drinking ~ RIDAGEYR + RIAGENDR + (DMDEDUC2<=3), data=predictors, family="binomial"))

hist(predict(glm(light_drinking ~ RIDAGEYR + RIAGENDR + as.factor(DMDEDUC2), data=predictors, family="binomial"), type="response"), breaks=50)


## treatment variable
tr <- X[, "light_drinking"]

Y <- outcomes %>% 
  filter(SEQN %in% predictors$SEQN) %>% 
  arrange(SEQN) %>% 
  dplyr::select(-SEQN) %>% 
  mutate_all(log) %>%
  as.matrix
colnames(Y) <- c("HDL", "Methylmercury", "Glucose", 
                 "Potassium", "Sodium", "Iron", "Lead", "LDL", "Triglycerides", "Cadmium")

Y[Y[, "Methylmercury"]==log(0.9), "Methylmercury"] <- NA
Y[Y[, "Cadmium"]==log(0.632), "Cadmium"] <- NA
bounds_mat <- matrix(c(2, -Inf, min(Y[, "Methylmercury"], na.rm=TRUE), 10, -Inf, min(Y[, "Cadmium"], na.rm=TRUE)), nrow=2, byrow=TRUE)
Y <- Amelia::amelia(Y, bounds=bounds_mat)$imp[[1]]

bart_res_hdl <- bcf::bcf(y=Y[, "HDL"], z=tr$light_drinking, x_control=as.matrix(X[,-1]), 
                     pihat = predict(glm(light_drinking ~ ., data=X, family=binomial("logit")), type="response"),
                     nburn=1000, nsim=2000)
bart_res_hg <- bcf::bcf(y=Y[, "Methylmercury"], z=tr$light_drinking, x_control=as.matrix(X[,-1]), 
                     pihat = predict(glm(light_drinking ~ ., data=X, family=binomial("logit")), type="response"),
                     nburn=1000, nsim=2000)

res <- rbart::rbart(x.train = as.matrix(X), y.train = Y[, "HDL"])



X %>% mutate(preds=apply(bart_res$tau, 2, mean), 
             lower=apply(bart_res$tau, 2, function(x) quantile(x, 0.025)),
             upper=apply(bart_res$tau, 2, function(x) quantile(x, 0.975))) %>% 
  mutate(gender=as.factor(ifelse(gender==1, "male", "female")), educ=as.factor(ifelse(educ==1, "ged", "no ged"))) %>%
  ggplot() + 
  geom_line(aes(x=age, y=preds)) + 
  geom_ribbon(aes(x=age, ymin=lower, ymax=upper), alpha=0.25) + 
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_wrap(~educ+gender)




yres <- (sapply(1:10, function(x) lm(scale(Y[, x]) ~ as.matrix(X))$residuals))
colnames(yres) <- colnames(Y)

as_tibble(yres) %>% pivot_longer(cols=1:10, values_to="value", names_to="outcome") %>% 
  ggplot() + geom_histogram(aes(x=value)) + facet_wrap(~outcome) + theme_bw()


pvals <- sapply(colnames(Y), function(x) coef(summary(lm(as.formula(sprintf("%s ~ light_drinking + age + gender + educ", x)), cbind(Y, X))))[2:5, 4]) < 0.05

rsq_full <- sapply(colnames(Y), function(x) summary(lm(as.formula(sprintf("%s ~ light_drinking + age + gender + educ", x)), data=cbind(Y, X)))$r.squared)
rsq_no_age <- sapply(colnames(Y), function(x) summary(lm(as.formula(sprintf("%s ~ light_drinking + gender + educ", x)), data=cbind(Y, X)))$r.squared)
rsq_no_gender <- sapply(colnames(Y), function(x) summary(lm(as.formula(sprintf("%s ~ light_drinking + age + educ", x)), data=cbind(Y, X)))$r.squared)
rsq_no_educ <- sapply(colnames(Y), function(x) summary(lm(as.formula(sprintf("%s ~ light_drinking + age + gender", x)), data=cbind(Y, X)))$r.squared)

cbind(rsq_no_age, rsq_no_gender, rsq_no_educ)


partial_gender <- (rsq_full - rsq_no_gender) / (1 -rsq_no_gender)
partial_age <- (rsq_full - rsq_no_age) / (1 -rsq_no_age)
partial_educ <- (rsq_full - rsq_no_educ) / (1 -rsq_no_educ)



sm <- stan_model("../multivariate_regression.stan")

for(m in 1:2) {
  data_list <- list(K=ncol(Y), J=ncol(X), N=nrow(Y), M=m, x=X, y=scale(Y))
  stan_results <- sampling(sm, data_list)
  saveRDS(stan_results, file=sprintf("drinking_samples_weduc_rank%i.RDS", m))
}



stan_results <- readRDS("drinking_samples_weduc_rank5.RDS")

naive_estimates <- stan_results %>% spread_draws(beta[K, J]) %>% 
  filter(.chain != 3) %>% 
  group_by(J, K) %>% 
  summarize(q025=quantile(beta, 0.025), 
            q975=quantile(beta, 0.975)) %>%
  mutate(sig = sign(q025*q975)==1) %>% 
  mutate(feature = colnames(Y)) %>% 
  filter(J==1) %>% 
  arrange(desc(q025)) %>% relocate(feature, q025, q975)
naive_estimates




############################################
###### This takes a few minutes ############
############################################


  
########### 
stan_results <- readRDS(file=sprintf("drinking_samples_weduc_rank%i.RDS", m))
data_list$M <- m
NC_index <- 2

r2seq <- seq(0, 0.1, by=0.0001)
odds <- r2_to_odds(r2seq, probs=preds_full)
r2_benchmark <- r2seq[which.min(abs(odds - 3.5))]

## No single outcome confounding
r2y <- rep(NA, 10)

beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=NC_index,
                               probs=preds_full, full_rank=FALSE,
                               R2_no_nc = r2_benchmark, R2_nc = 0.0,
                               r2y = r2y)
effect_intervals <- compute_effect_intervals(beta_w_norm)

nc_plot_index <- which(effect_intervals$K %in% c(NC_index))  
make_interval_plot(effect_intervals, nc_plot_index=nc_plot_index,
                   title="Factor Confounding",
                   labels=c("NUC", expression(Lambda[0.95]~"="~Lambda[min]))) -> plt

plt + 
  scale_color_manual(name="", labels=c("NUC", expression(Lambda[0.95]~"= 3.5"), expression(Lambda[0.95]~"="~Lambda[min])), 
                    values=c(nuc="black", rect="black", 'null'="#CC1C2F")) + 
  theme(legend.key = element_rect(color = c(NA, "black", NA), fill=c(NA, NA, NA)))

ggsave(file=sprintf("effect_intervals_weduc_%s_odds_m%i.pdf", paste(NC_index, collapse="_"), m), width=10, height=5)

beta_w_norm %>% ggplot() + geom_point(aes(x=rv_init, y=rv_nc, col=as.factor(sign(rv_init)*sign(rv_nc)))) + facet_wrap(~K) + theme_bw() +
  geom_abline(slope=1, linetype="dashed") +
  geom_abline(slope=-1, linetype="dashed")

beta_w_norm %>% ggplot() + geom_density_2d(aes(x=abs(rv_init), y=abs(rv_nc)), contour_var="ndensity") + 
    geom_abline(slope=1, linetype="dashed") + facet_wrap(~K)



## R2y=1 for mercury
NC_index <- 2
r2y <- rep(NA, 10)
r2y[2] <- 1

beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=NC_index,
                               probs=preds_full, full_rank=TRUE,
                               R2_no_nc = r2_benchmark, R2_nc = 0.0,
                               r2y = r2y)
effect_intervals <- compute_effect_intervals(beta_w_norm)
nc_plot_index <- which(effect_intervals$K %in% c(NC_index))  
make_interval_plot(effect_intervals, nc_plot_index=nc_plot_index, 
                   title=expression(R[Y~"~"~"U|X"]^2~"=1 for methylmercury"),
                   labels=c("NUC", expression(Lambda[0.95]~"="~Lambda[min]))) -> plt

plt + 
  scale_color_manual(name="", labels=c("NUC", expression(Lambda[0.95]~"= 3.5"), expression(Lambda[0.95]~"="~Lambda[min])), 
                     values=c(nuc="black", rect="black", 'null'="#CC1C2F")) + 
  theme(legend.key = element_rect(color = c(NA, "black", NA), fill=c(NA, NA, NA)))

ggsave(file=sprintf("effect_intervals_%s_odds_m%i_null_r2y1.pdf", paste(NC_index, collapse="_"), m), width=10, height=5)

## R2y=1 for all
r2y <- rep(NA, 10)
#r2y[2] <- 0.2
r2y[4] <- 0.2
NC_index <- c(2)
beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=NC_index,
                               probs=preds_full, full_rank=TRUE,
                               r2y = r2y)
effect_intervals <- compute_effect_intervals(beta_w_norm)
nc_plot_index <- which(effect_intervals$K %in% c(NC_index))  
make_interval_plot(effect_intervals, nc_plot_index=nc_plot_index)

ggsave(file=sprintf("effect_intervals_%s_odds_m%i_all_blah.pdf", 
                    paste(NC_index, collapse="_"), m), width=10, height=5)

effect_intervals %>% filter(K %in% c(2, 4)) %>% glimpse

## All except HDL are null, r2y=1 for all
r2y <- rep(1, 10)
beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=2:10,
                               probs=preds_full, full_rank=TRUE,
                               r2y = r2y)
effect_intervals <- compute_effect_intervals(beta_w_norm)
nc_plot_index <- which(effect_intervals$K %in% 2:10)  
make_interval_plot(effect_intervals, nc_plot_index=nc_plot_index)

ggsave(file=sprintf("effect_intervals_%s_odds_m%i_full_rank_all_null.pdf", 
                    paste(NC_index, collapse="_"), m), width=10, height=5)

## All except HDL are null, r2y=1 for all
r2y <- rep(NA, 10)
beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=2:10,
                               probs=preds_full, full_rank=TRUE,
                               r2y = r2y)
effect_intervals <- compute_effect_intervals(beta_w_norm)
nc_plot_index <- which(effect_intervals$K %in% 2:10)  
make_interval_plot(effect_intervals, nc_plot_index=nc_plot_index)







###

NC_index <- 2
r2y_seq <- c(seq(1, 0.2, by=-0.1), seq(0.2, 0.1, by=-0.02))
effect_intervals_joined <- purrr::map_dfr(r2y_seq, function(x) {
  print(x)
  r2y <- rep(NA, 10)
  r2y[2] <- x
  beta_w_norm <- get_beta_w_norm(stan_results, data_list, NC_index=NC_index,
                               probs=preds_full, full_rank=TRUE,
                               r2y = r2y)
  effect_intervals <- compute_effect_intervals(beta_w_norm)
  effect_intervals$r2ch3hg <- x
  effect_intervals
  }
)

nc_color = "#CC1C2F"
effect_intervals_joined %>% filter(K==1) %>% 
  ggplot(aes(x=r2ch3hg)) + geom_ribbon(aes(ymin=lower_nc, ymax=upper_nc), fill=nc_color, alpha=0.5) + 
  theme_bw() + geom_hline(yintercept=0, linetype='dashed') + 
  xlim(c(0, 1))

orig_color = "#023FA5"
nc_color = "#CC1C2F"
effect_intervals_joined %>%
  filter(r2ch3hg <= 0.35) %>% 
  ggplot() + 
  geom_segment(aes(x=forcats::fct_inorder(Metabolite), xend=forcats::fct_inorder(Metabolite), 
                   y=lower_nc, yend=upper_nc, col='null'), 
               size=3, alpha=0.5) +
  geom_segment(aes(x=forcats::fct_inorder(Metabolite), xend=forcats::fct_inorder(Metabolite), 
                   y=lower_init, yend=upper_init, col='nuc'), size=1.5) +
  transition_states(as.factor(r2ch3hg), state_length=0) +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  geom_hline(yintercept=0, linetype="dashed") +  
  geom_point(aes(x=nc_plot_index, y=rep(0, length(nc_plot_index))), size=3, colour=nc_color, alpha=0.5, data=tibble()) +
  geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(1, 60),
                label=ifelse(sign(lower_naive) == sign(upper_naive), format(round(lower_odds_init, digits=1), nsmall=1), "X")), fontface="bold") +
  geom_text(aes(x=forcats::fct_inorder(Metabolite), y = rep(0.9, 60),
                label=ifelse(sign(lower_nc) == sign(upper_nc), format(round(lower_odds_nc, digits=1), nsmall=1), "X")), col=nc_color, fontface="bold") +
  xlab("") + ylab("Treatment Effect") + 
  scale_color_manual(name="", labels=c(expression(R[paste(T,'~',U,'|',X)]^2~"=0"), expression(R[paste(T,'~',U, '|', X)]^2~"="~R[paste("min")]^2)), 
                     values=c(nuc="black", 'null'=nc_color)) + 
  theme(legend.text.align = 0) +
  labs(title = "R2_Y Ch3Hg: {closest_state}") -> intervals_animation

anim_save("test.gif", animate(intervals_animation, width=10, height=5, units = "in", res = 150, rewind=TRUE))



# beta_w_norm %>% ggplot() +
# ggridges::geom_density_ridges(aes(x=log(abs(odds_init)), y=as.factor(K), fill=as.factor(K)), alpha=0.5,
#                               quantile_lines=TRUE,
#                               quantile_fun=function(x,...) median(x, na.rm=TRUE)) + theme_bw(base_size=16) +
#   theme(legend.position="none") + xlim(c(0, 5))
#   xlab("Robustness Value") + ylab("") + ggtitle("Posterior Robustness")

effect_intervals_joined %>% filter(K==1) %>%
  mutate(covers=ifelse(sign(lower_nc)!=sign(upper_nc), "YES", "NO")) %>%
  mutate(Effect=(lower_nc+(upper_nc-lower_nc)/2)) %>%
  filter(covers=="NO") %>% 
  ggplot(aes(y=lower_odds_nc, x=r2ch3hg)) + 
  geom_line() + 
  geom_rect(xmin=0, xmax=0.13, ymin=0, ymax=10, fill="red", alpha=0.02) +  
  geom_point(aes(col=Effect), size=3) +
  ylim(c(1, 8)) + xlim(c(0, 1)) + 
  ylab(expression(Lambda[0.95]~" for HDL")) + xlab(expression(R[Y~"~"~"U|X"]^2~" for methylmercury"))+
  theme_bw(base_size=16) +

  ggrepel::geom_label_repel(aes(x=x, y=y, label=label), min.segment.length = 0, seed = 42, box.padding = 0.5, 
                            data=tibble(x=c(partial_age[2], partial_gender[2], partial_educ[2]), 
                                        y=c(3.5, 1.5, 1.5), 
                                        label=c("Age", "Gender", "Education"))) +
  colorspace::scale_color_continuous_sequential(palette="Viridis")



effect_intervals_joined %>% filter(K==4) %>%
  mutate(covers=ifelse(sign(lower_nc)!=sign(upper_nc), "YES", "NO")) %>%
  mutate(midpoint=(lower_nc+(upper_nc-lower_nc)/2)) %>%
  filter(covers=="NO") %>% 
  ggplot(aes(y=lower_odds_nc, x=r2ch3hg)) + 
  geom_line() + 
  geom_point(aes(col=midpoint), size=3) +
  ylim(c(1, 8)) + xlim(c(0, 1)) + 
  ylab(expression(Lambda[0.95]~" for HDL")) + xlab(expression(R[Y~"~"~"U|X"]^2~" for methylmercury"))+
  theme_bw(base_size=16) +
  ggrepel::geom_label_repel(aes(x=x, y=y, label=label), min.segment.length = 0, seed = 42, box.padding = 0.5, 
                            data=tibble(x=c(partial_age[2], partial_gender[2], partial_educ[2]), 
                                        y=c(3.5, 1.5, 1.5), 
                                        label=c("Age", "Gender", "Education"))) +
  colorspace::scale_color_continuous_sequential(palette="Viridis",)








NC_index <- 2
r2y_seq <- c(seq(0.05, 1, by=0.05))

r2y_grid <- as.matrix(expand.grid(x=r2y_seq, y=r2y_seq))

## Faster for run time.  Use posterior mean of B and lambda and 0.05 quantile of beta for rv


stan_results %>% spread_draws(B[K, M], lambda[K]) %>% 
  nest(data = c(B, K, M, lambda)) %>% 
  { map(.$data, function(x) {
    
    B <- matrix(x$B, nrow=max(x$K), ncol=max(x$M), byrow=TRUE)
    lambda <- x %>% filter(M==1) %>% pull(lambda) %>% diag
    
    #B %*% t(B)
    
    B %*% t(B) #+ diag(lambda)
    
  }) } -> BBlist

Bstar <- t(chol(Reduce('+', BBlist)/4000))[, 1:data_list$M]

df_mean <- stan_results %>% spread_draws(B[K, M], lambda[K], beta[K, J]) %>%
  group_by(K, M, J) %>%
  summarize(B = mean(B), lambda=mean(lambda), beta_lower=quantile(beta, 0.025), beta_upper=quantile(beta, 0.975)) %>%
  filter(J==1)
df_mean$B <- as.numeric(t(Bstar))

k_odds <- numeric(nrow(r2y_grid))
for(i in 1:nrow(r2y_grid)) {

  x <- r2y_grid[i, ]

  r2y <- rep(NA, 10)
  r2y[2] <- x[1]
  r2y[1] <- x[2]
  
  print(c(r2y[2], r2y[4]))
    
  upper <- cal_tau_cali_wNC(df_mean %>%  mutate(beta=beta_lower), R2=0.0, NC_index=NC_index, 
                   treatment_index = 1, data_list=data_list, 
                   marginal_sensitivity_param = TRUE, 
                   probs=preds_full, full_rank=TRUE, r2y=r2y)
  lower <- cal_tau_cali_wNC(df_mean %>%  mutate(beta=beta_upper), R2=0.0, NC_index=NC_index, 
                   treatment_index = 1, data_list=data_list, 
                   marginal_sensitivity_param = TRUE, 
                   probs=preds_full, full_rank=TRUE, r2y=r2y)
  
  upper4 <- upper %>% filter(K==1) %>% select(R2_nc_lwr, odds_nc) 
  lower4 <- lower %>% filter(K==1) %>% select(R2_nc_lwr, odds_nc)

  k_odds[i] <- ifelse(sign(upper4$R2_nc_lwr) != sign(lower4$R2_nc_lwr), NA, min(upper4$odds_nc, lower4$odds_nc))
}


tibble(r2ch3hg=r2y_grid[, 1], r2k=r2y_grid[, 2], odds=k_odds) %>% 
  ggplot(aes(x=r2ch3hg, y=r2k)) + 
  geom_tile(aes(fill=odds), color="black") +
  geom_text(aes(label=round(odds, 1)), col="black") +   
  colorspace::scale_fill_continuous_sequential(palette="ag_GrnYl", na.value="white") + theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank()) + 
  xlab(expression(R["Y~X|U"]^2~" for methylmercury")) + ylab(R["Y~X|U"]^2~" for potassium") + ggtitle(expression(Lambda)) +
  theme(text = element_text(size=20))


tibble(r2ch3hg=r2y_grid[, 1], r2k=r2y_grid[, 2], odds=k_odds) %>% 
  ggplot(aes(x=r2ch3hg, y=r2k)) + 
  stat_contour(aes(z=odds)) +
  metR::geom_text_contour(aes(z=odds)) +
   theme_bw() +
  xlab(expression(R["Y~X|U"]^2~" for methylmercury")) + ylab(R["Y~X|U"]^2~" for potassium") + ggtitle(expression(Lambda)) +
  theme(text = element_text(size=20))

  
  

#####################################################################

rank1 <- readRDS(file="drinking_samples_weduc_rank1.RDS")
loo1 <- loo(rank1)


rank2 <- readRDS(file="drinking_samples_weduc_rank2.RDS")
loo2 <- loo(rank2)

rank3 <- readRDS(file="drinking_samples_weduc_rank3.RDS")
loo3 <- loo(rank3)
  
rank4 <- readRDS(file="drinking_samples_weduc_rank4.RDS")
loo4 <- loo(rank4)

rank5 <- readRDS(file="drinking_samples_weduc_rank5.RDS")
loo5 <- loo(rank5)

rank6 <- readRDS(file="drinking_samples_weduc_rank6.RDS")
loo6 <- loo(rank6)

rank7 <- readRDS(file="drinking_samples_weduc_rank7.RDS")
loo7 <- loo(rank7)
 
rank8 <- readRDS(file="drinking_samples_weduc_rank8.RDS")
loo8 <- loo(rank8)

rank10 <- readRDS(file="drinking_samples_weduc_rank10.RDS")
loo10 <- loo(rank10)

r2y_mat <- sapply(list("rank4"=rank4, "rank5"=rank5, "rank6"=rank6, "rank7"=rank7, "rank8"=rank8), 
       function(x) apply(matrix(unlist(compute_r2y(x)), nrow=10), 1, median))

rownames(r2y_mat) <- colnames(Y)
as_tibble(r2y_mat, rownames="outcome") %>% pivot_longer(cols=-outcome, names_to="rank") %>%
  ggplot() + geom_line(aes(x=as.numeric(as.factor(rank))+3, y=value, col=as.factor(outcome))) +
  xlab("Rank") + ylab("r2y") + theme_bw() + ylim(c(0,1))

r2y_nsoc <- apply(matrix(unlist(compute_r2y(stan_results)), nrow=10), 1, median)



rsq_x_mat <- apply(round(cbind(partial_age, partial_gender, partial_educ), digits=3), 2, as.character)
for(i in 1:3) {
  sig_indices <- which(!t(pvals)[, i+1])
  rsq_x_mat[sig_indices, i] <- "NS"
}
rsq_x_mat <- cbind(rsq_x_mat, round(r2y_nsoc, 3))
rownames(rsq_x_mat) <- colnames(Y)
xtable::xtable(rsq_x_mat)








comp <- loo::loo_compare(x=list("rank1"=loo1, "rank2"=loo2, "rank3"=loo3, "rank4"=loo4, "rank5"=loo5, 
                                "rank6"=loo6, "rank7"=loo7, "rank8"=loo8, "rank10"=loo10))
print(comp, digits=2)
xtable::xtable(comp[, 1:2])

loo::loo_compare(x=list("rank7"=loo7, "rank8"=loo8))




stan_results <- readRDS(file=sprintf("drinking_samples_full_rank%i.RDS", data_list$M))
stan_results %>% spread_draws(B[K, M], lambda[K]) %>% 
  nest(data = c(B, K, M, lambda)) %>% 
  { map(.$data, function(x) {

    B <- matrix(x$B, nrow=max(x$K), ncol=max(x$M), byrow=TRUE)
    lambda <- x %>% filter(M==1) %>% pull(lambda) %>% diag

    #B %*% t(B)

    B %*% t(B) #+ diag(lambda)
    
  }) } -> BBlist

eig <- eigen(Reduce('+', BBlist)/4000)

summary(sapply(BBlist, function(x) cov2cor(x)[1, 2]))

gamma_cols <- data_list$M #data_list$K
gamma_init <- eig$vectors[, 1:gamma_cols] %*% sqrt(diag(eig$values[1:gamma_cols]))
gamma_init <- varimax(gamma_init)$loadings

g2 <- gamma_init[2, ]
g2 <- g2/sqrt(sum(g2^2))
n2 <- rstiefel::NullC(g2)

g1 <-  n2 %*% t(n2) %*% gamma_init[1, ]
g1 <- g1/sqrt(sum(g1^2))

n21 <- rstiefel::NullC(cbind(g2, g1))

g8 <-  n21 %*% t(n21) %*% gamma_init[8, ]
g8 <- g8/sqrt(sum(g8^2))

n218 <- rstiefel::NullC(cbind(g2, g1, g8))

g4 <- n218 %*% t(n218) %*% gamma_init[4, ]
g4 <- g4/sqrt(sum(g4^2))

n2184 <- rstiefel::NullC(cbind(g2, g1, g8, g4))

g3 <- n2184 %*% t(n2184) %*% gamma_init[3, ]
g3 <- g3/sqrt(sum(g3^2))
nrest <- rstiefel::NullC(cbind(g2, g1, g8, g4, g3))

#gamma <- gamma_init %*% cbind(g2, g1, g8, g4, g3, nrest)
gamma <- gamma_init %*% cbind(g2, g1, g8, g4, g3)

ord <- c(2, 1, 8, 4, 3, 6, 9, 7, 10, 5)
gamma <- gamma[ord, ]



tibble(expand.grid(outcome=colnames(Y)[ord], col=1:gamma_cols)) %>% 
  mutate(Value = as.numeric(gamma)) %>%  
  ggplot(aes(x=col, y=forcats::fct_rev(as.factor(outcome)))) + 
  geom_tile(aes(fill=Value), col="black") + 
  geom_text(aes(label=round(Value, 2)), col="black") +   
  colorspace::scale_fill_continuous_diverging() + theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + ylab("") + ggtitle(expression(Gamma)) +
  theme(text = element_text(size=20)) -> gamma_fig

ggsave(gamma_fig, filename=sprintf("nhanes_heat_rank_%i.pdf", data_list$M), width=6, height=5)

stan_results %>% spread_draws(lambda[K]) %>% 
  nest(data = c(K, lambda)) %>% 
  { map(.$data, function(x) {

    lambda <- x %>% pull(lambda)

  }) } -> Lambda_list

lambda_mean <- Reduce('+', Lambda_list)/4000
lambda_mat <- diag(lambda_mean)

tibble(expand.grid(outcome=colnames(Y)[ord], col=1)) %>% 
  mutate(Value = as.numeric(lambda_mean[ord])) %>%  
  ggplot(aes(x=col, y=forcats::fct_rev(as.factor(outcome)))) + 
  geom_tile(aes(fill=Value), col="black") + 
  geom_text(aes(label=round(Value, 2)), col="black") +   
  colorspace::scale_fill_continuous_diverging() + theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + ylab("") + ggtitle(expression(Lambda)) +
  theme(text = element_text(size=20)) -> lambda_fig


ggsave(lambda_fig, filename="nhanes_lambda_5.pdf", width=3.5, height=5)

library(patchwork)
ggsave(gamma_fig + lambda_fig, filename="nhanes_gamma_lambda_5.pdf", width=10, height=5)
