# Model Comparation

library(loo)
library(ggplot2)

exp_model <- readRDS("exp_PH_model_draws.rds")

weibull_model <- readRDS("weibull_PH_model_draws.rds")

semiparametric_model <- readRDS("semiparametric_model_fit.rds")
gc()
# Extract log-likelihoods
log_lik_exp_draws <- exp_model$loo(cores = 4)
log_lik_weibull_draws <- weibull_model$loo(cores = 4)
log_lik_semiparametric_draws <- semiparametric_model$loo(cores = 4)

# Compare models using LOO

comp <- loo_compare(
  log_lik_exp_draws,
  log_lik_weibull_draws,
  log_lik_semiparametric_draws)

comp <- as.data.frame(comp)

write.csv(comp, "model_comparation.csv")

# Plot LOO comparison
rownames(comp) <- c("Semiparametric", "Weibull", "Exponential")

ggplot(comp, aes(x = 1:3, y = looic, 
                 fill =c("Semiparametric", "Weibull", "Exponential"))) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = looic - se_looic, ymax = looic + se_looic), width = 0.2) +
  labs(x = "", 
       y = "LOOIC",
       fill = "Model",
       ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
