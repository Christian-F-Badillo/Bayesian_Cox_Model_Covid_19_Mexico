# Librerias

## Bayes libreria
library(cmdstanr)
library(posterior)
options(mc.cores = parallel::detectCores())

cmdstan_path()
# Debes ser "C:/Users/crish/.cmdstan/cmdstan-2.36.0"

# cmdstan version
cmdstan_version()

# Survival analysis
library(survival)

# Gráficos
library(repr)
library(bayesplot)
library(ggplot2)
library(survminer)
color_scheme_set("teal")


## Menjo de Datos
library(dplyr)
library(tidyr)
library(reshape2)

#------------------------------------------------------------------
# Arreglo de Datos
datos <- read.csv("COVID19MEXICO2022_muestra_clean.csv")

cat(names(datos), sep = " + " )

# Quitamos lo que no se va a usar
data_ord <- datos %>%
    select(-c(fecha_ingreso,fecha_sintomas,fecha_def,fin_ano, tipo_paciente))

head(data_ord)
tail(data_ord)

# Verificamos NA
print(apply(data_ord, 2, function(x) sum(is.na(x)) ))

# Quitmaos los na
data_ord <- data_ord %>%
    drop_na()

# Verificamos NA
print(apply(data_ord, 2, function(x) sum(is.na(x)) ))

# TIEMPO y EVENTOS
time.data <- data_ord$tiempo
event.data <- as.numeric(as.logical(data_ord$event))

# Convertimos a Factor
data_ord <- data_ord %>%
    mutate(across(
        c(sexo, neumonia, diabetes, epoc, asma,
          inmusupr, hipertension, otra_com, cardiovascular,
          obesidad, renal_cronica, tabaquismo, otro_caso),
        ~ as.factor(.x)
    ))


#--------------------------------------------------------------------------------
# Analísis Inferencial de los Datos Frecuentista
surv_model <- Surv(time.data, event.data) ~ sexo + neumonia + edad + diabetes + epoc + asma + inmusupr + hipertension + otra_com + cardiovascular + obesidad + renal_cronica + tabaquismo + otro_caso

# Ajustamos el modelo
fit <- coxph(surv_model, data = data_ord)

summary(fit)

curva <- survfit(fit)

ggsurvplot(curva, data = data_ord,
           xlab = "Tiempo",
           ylab = "Supervivencia estimada",
           conf.int = TRUE,
           risk.table = TRUE,
           ggtheme = theme_minimal())

# Guardamos el modelo
saveRDS(curva, file = "curva_fit_frec.rds")
saveRDS(fit, file = "modelo_cox_frec.rds")

# ---------------------------------------------------------------------------------
# Bayesian Cox Model

# Definimos el modelo Exponencial de Riesgos Proporcionales.

## Modelo en Stan

write(
"
data {
  int<lower=1> N;                       // Número de observaciones
  int<lower=1> K;                       // Número de covariables
  vector[N] time;                       // Tiempos observados (evento o censura)
  array[N] int<lower=0,upper=1> event;  // Indicador de evento (1=evento, 0=censura)
  matrix[N, K] X;                       // Matriz de covariables (sin intercepto)
}

parameters {
  vector[K] beta;                // Coeficientes de las covariables
  real<lower=0> lambda;          // Parámetro de riesgo base
}

model {
  // Priors
  beta ~ normal(0, 2);           // Prior para coeficientes
  lambda ~ gamma(1, 4);          // Prior para el riesgo base

  // Vectorizada (más eficiente)
  vector[N] rate = lambda * exp(X * beta);
  vector[N] log_rate = log(rate);
  vector[N] rate_time = rate .* time;

  vector[N] event_vector = to_vector(event); // ✅ Conversión explícita
  vector[N] event_complement = 1 - event_vector; // ✅ Operación válida

  target += dot_product(event_vector, log_rate - rate_time) + 
            dot_product(event_complement, -rate_time);
  
}

generated quantities {
  vector[N] log_lik;  // Log-verosimilitud

  // Convertir 'event' a vector
  vector[N] event_vec = to_vector(event);
  vector[N] event_comp = 1 - event_vec;

  // Términos vectorizados
  vector[N] rate = lambda * exp(X * beta);
  vector[N] log_rate = log(rate);
  vector[N] rate_time = rate .* time;

  // Log-verosimilitud
  log_lik = event_vec .* (log_rate - rate_time) + event_comp .* (-rate_time);
}
", 
    file = "exp_PH_model.stan"
)
rm(fit)
rm(curva)

# Definimos los datos para el modelo
stan_data <- list(
    N = nrow(data_ord),
    K = ncol(data_ord) - 3, # Sin tiempo y evento
    time = time.data,
    event = event.data,
    X = data_ord[, -c(15, 16, 17)] # Sin tiempo y evento
)

# Definimos el modelo
exp_PH_model <- cmdstan_model("exp_PH_model.stan")
exp_PH_model$print(line_numbers = TRUE)

# Ajustamos el modelo
fit_bayes <- exp_PH_model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    save_warmup = T,
    iter_sampling = 1000,
    refresh = 100,
    seed = 26534167
)


# Guardamos el modelo
fit_bayes$save_object(file = "exp_PH_model_draws.rds")

# Abrimos el modelo
fit_bayes <- readRDS("exp_PH_model_draws.rds")

# Diagnóstico
fit_bayes$cmdstan_diagnose()

# Resumen de los resultados
sum_draws <- fit_bayes$summary(
    variables = c("beta", "lambda"),
    posterior::default_summary_measures(),
    default_convergence_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
    )

write.csv(sum_draws, file = "summary_exp_PH_model.csv")

# Param Names
betas_str <- character(length =  ncol(data_ord) - 3)
for (i in 1:( ncol(data_ord) - 3)) {
    betas_str[i] <- paste("beta[", i, "]", sep = "")
}

betas_str

# Cadenas
p <- mcmc_trace(
    fit_bayes$draws(),
    regex_pars  = c("beta", "lambda"),
    n_warmup = 500,
    facet_args = list(nrow = 5, ncol = 3, labeller = label_parsed)
)
p + facet_text(size = 15)

# Gráficos Posteriores.
mcmc_areas_theme <- theme(
    text = element_text(size = 16),           # Cambia todo el texto
    axis.title = element_text(size = 18),     # Título de ejes
    axis.text = element_text(size = 14),      # Etiquetas de ticks
    axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5),   # Girar las etiquetas del eje Y
    legend.text = element_text(size = 14),    # Texto de leyenda
    legend.title = element_text(size = 16),   # Título de leyenda
    plot.title = element_text(size = 30, face = "bold")  # Título del gráfico
)

mcmc_areas(fit_bayes$draws(), 
           regex_pars = c("beta", "lambda"),
           prob_outer = 1, 
           point_est = "mean", 
           prob = 0.95,
           area_method = "equal height") +
    scale_y_discrete(labels = parse(text = c(betas_str, "lambda"))) +
    mcmc_areas_theme

## Correlación entre betas
samples_df <- fit_bayes$draws(variables = "beta", format = "df") # or format="array"

samples_df <- samples_df %>% select(-c(.chain, .iteration, .draw))

# Calcular la matriz de correlación
cor_matrix <- cor(samples_df)

# Convertir la matriz de correlación a formato largo
cor_melt <- melt(cor_matrix)

ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", 
                         mid = "white", midpoint = 0, 
                         limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.15),
          axis.title = element_text(size = 18),     # Título de ejes
          axis.text = element_text(size = 16)) +
    labs(title = "", x = "", y = "") + 
    scale_y_discrete(labels = parse(text = betas_str)) + 
    scale_x_discrete(labels = parse(text = betas_str))
rm(p)

# Curva de Sobrevivencia
         