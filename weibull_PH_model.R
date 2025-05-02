# Librerias

## Bayes libreria
library(cmdstanr)
library(posterior)
library(HDInterval)
options(mc.cores = parallel::detectCores())

cmdstan_path()
# Debes ser "C:/Users/crish/.cmdstan/cmdstan-2.36.0"
set_cmdstan_path("C:/Users/crish/.cmdstan/cmdstan-2.36.0")

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
        c(sexo, intubado, neumonia, diabetes, epoc, asma,
          inmusupr, hipertension, otra_com, cardiovascular,
          obesidad, renal_cronica, tabaquismo, otro_caso),
        ~ as.factor(.x -1)
    ))

summary(data_ord)

# Los tiempos deed 0 se ponen con 0.01
time.data[time.data == 0] <- 0.01
summary(time.data)

#---------------------------------------------------------------------------------
# Bayesian Weibull Proportional Hazards Model

## Modelo en Stan
write("
data {
  int<lower=1> N;                       // Número de observaciones
  int<lower=1> K;                       // Número de covariables
  vector[N] time;                       // Tiempos observados
  array[N] int<lower=0, upper=1> event; // Indicador de evento
  matrix[N, K] X;                       // Matriz de diseño (sin intercepto)
}

parameters {
  real lambda;               // Parámetro de escala (log(gamma))
  real<lower=0> alpha;       // Parámetro de forma
  vector[K] beta;            // Coeficientes de las covariables
}

transformed parameters {
  vector[N] sigma;           // Parámetro de escala de Stan por observación
  sigma = exp(-(lambda + X * beta)/alpha);
}

model {
  // Priors
  lambda ~ normal(10, 5);     // Prior normal para lambda
  alpha ~ gamma(1, 1);       // Prior gamma para alpha
  beta ~ normal(0, 3);       // Prior para coeficientes
  
  // Likelihood
  for (i in 1:N) {
    if (event[i]) {
      target += weibull_lpdf(time[i] | alpha, sigma[i]);
    } else {
      target += weibull_lccdf(time[i] | alpha, sigma[i]);
    }
  }
}

generated quantities {
  real gamma = exp(lambda);  // Parámetro de escala en parametrización de supervivencia
  vector[N] log_lik;         // Log-verosimilitud punto a punto
  
  for (i in 1:N) {
    // Usamos el sigma ya calculado en transformed parameters
    real sigma_i = sigma[i];
    
    if (event[i]) {
      log_lik[i] = weibull_lpdf(time[i] | alpha, sigma_i);
    } else {
      log_lik[i] = weibull_lccdf(time[i] | alpha, sigma_i);
    }
  }
}

", file = "weibull_PH_model.stan")

#---------------------------------------------------------------------------------
# Definimos los datos para el modelo
stan_data <- list(
    N = nrow(data_ord),
    K = ncol(data_ord) - 2, # Sin tiempo y evento
    time = time.data,
    event = event.data,
    X = data_ord[, -c(16, 17)] # Sin tiempo y evento
)

# Definimos el modelo
weibull_PH_model <- cmdstan_model("weibull_PH_model.stan")
weibull_PH_model$print(line_numbers = TRUE)

# Ajustamos el modelo
fit <- weibull_PH_model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    save_warmup = T,
    iter_sampling = 1000,
    refresh = 100,
    seed = 1234
)


# Guardamos el modelo
fit$save_object(file = "weibull_PH_model_draws.rds")

# Abrimos el modelo
fit <- readRDS("weibull_PH_model_draws.rds")

# Diagnóstico
fit$cmdstan_diagnose()
fit$print(
    c("lambda", "alpha", "beta"),
    digits = 4
)
  
# Resumen de los resultados
sum_draws <- fit$summary(
    variables = c("beta", "lambda", "alpha"),
    posterior::default_summary_measures(),
    default_convergence_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)

write.csv(sum_draws, file = "summary_weibull_PH_model.csv")

# Param Names
betas_str <- character(length =  ncol(data_ord) - 2)
for (i in 1:( ncol(data_ord) - 2)) {
    betas_str[i] <- paste("beta[", i, "]", sep = "")
}

betas_str

draws_parameters <- fit$draws(variables = c("beta", "lambda", "alpha"), 
                                    inc_warmup  = F,
                                    format = "df") # or format="array"

# Cadenas
p <- mcmc_trace(
    draws_parameters,
    regex_pars  = c("beta", "lambda", "alpha"),
    facet_args = list(nrow = 6, ncol = 3, labeller = label_parsed)
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

mcmc_areas(draws_parameters, 
           regex_pars = c("beta", "lambda", "alpha"),
           prob_outer = 1, 
           point_est = "mean", 
           prob = 0.95,
           area_method = "equal height") +
    scale_y_discrete(labels = parse(text = c(betas_str, "lambda", "alpha"))) +
    mcmc_areas_theme

# Calcular la matriz de correlación
cor_matrix <- cor(draws_parameters %>% select(-c(.chain, .iteration, .draw)))

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
    scale_y_discrete(labels = parse(text = c(betas_str, "lambda", "alpha"))) + 
    scale_x_discrete(labels = parse(text = c(betas_str, "lambda", "alpha")))
rm(cor_melt)



#--------------------------------------------------------------------------------
# Graficar Curvas de Supervivencia Posteriores

plot_posterior_survival_weibull <- function(posterior_draws, covariates, t_max = 365, n_curves = 100) {
  # Extraer parámetros del modelo Weibull
  lambda <- as.numeric(posterior_draws$lambda)
  alpha <- as.numeric(posterior_draws$alpha)
  betas <- as.matrix(posterior_draws %>% select(matches("^beta\\[.*\\]")))
  
  # Validar dimensiones
  if (ncol(betas) != length(covariates)) {
    stop("Número de betas y covariables no coincide")
  }
  
  # 1. Calcular predictor lineal (lambda + Xβ)
  linpred <- betas %*% covariates  # Vector de tamaño n_muestras
  
  # 2. Calcular parámetro de escala sigma = exp(-(lambda + Xβ)/alpha)
  sigma_vec <- exp(-(lambda + linpred)/alpha)  # Vector de tamaño n_muestras
  
  # 3. Secuencia de tiempo
  t_seq <- seq(0, t_max, length.out = 1000)
  
  # 4. Calcular matriz de supervivencia (n_muestras × n_tiempos)
  t_matrix <- matrix(t_seq, nrow = length(sigma_vec), ncol = length(t_seq), byrow = TRUE)
  sigma_matrix <- matrix(sigma_vec, nrow = length(sigma_vec), ncol = length(t_seq))
  alpha_matrix <- matrix(alpha, nrow = length(alpha), ncol = length(t_seq))
  
  S_matrix <- exp(-(t_matrix/sigma_matrix)^alpha_matrix)
  
  # 5. Estadísticos resumen
  median_S <- apply(S_matrix, 2, median)
  hdi_S <- apply(S_matrix, 2, function(x) HDInterval::hdi(x, credMass = 0.9))
  
  # 6. Dataframe para gráfico
  plot_df <- data.frame(
    time = t_seq,
    median = median_S,
    lower = hdi_S[1, ],
    upper = hdi_S[2, ]
  )
  
  # 7. Muestrear curvas posteriores
  set.seed(123)
  idx <- sample(nrow(S_matrix), min(n_curves, nrow(S_matrix)))
  posterior_curves <- data.frame(
    time = rep(t_seq, each = length(idx)),
    survival = as.vector(S_matrix[idx, ]),
    draw = rep(1:length(idx), times = length(t_seq))
  )
  
  # 8. Graficar
  ggplot(plot_df, aes(x = time)) +
    geom_line(
      data = posterior_curves,
      aes(y = survival, group = draw),
      color = "skyblue",
      alpha = 0.1,
      linewidth = 0.3
    ) +
    geom_line(aes(y = median), color = "blue", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
    labs(
      title = "Curvas de Supervivencia Posteriores (Weibull)",
      x = "Tiempo",
      y = "Probabilidad de Supervivencia"
    ) +
    ylim(0, 1) +
    theme_minimal()
}

# Definimos los valores de las covariables
data_summ <- data_ord %>%
    summarise(across(
        c(sexo, intubado, neumonia, edad, diabetes, epoc, asma,
          inmusupr, hipertension, otra_com, cardiovascular,
          obesidad, renal_cronica, tabaquismo, otro_caso),
        ~ mean(as.numeric(.x))))

# A excepción de la edad, todas las demás serán redondeadas a su valor entero
# más cercano.
data_summ <- data_summ %>%
    mutate(across(
        c(sexo, intubado, neumonia, diabetes, epoc, asma,
          inmusupr, hipertension, otra_com, cardiovascular,
          obesidad, renal_cronica, tabaquismo, otro_caso),
        ~ ifelse(.x > 1.5, 1, 0) ))

covariables <- as.vector(data_summ) %>% unlist()
summ_pars <- sum_draws %>%
    select(variable, mean, q5, q95) %>%
    rename(parameter = variable) %>%
    rename(hdi_lower = q5, hdi_upper = q95)

plot_posterior_survival_weibull(
    posterior_draws = draws_parameters,
    covariates = covariables,
    t_max = 365,
    n_curves = 200  # Número de curvas posteriores a mostrar
)