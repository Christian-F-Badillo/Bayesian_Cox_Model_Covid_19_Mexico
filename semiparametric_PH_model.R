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

event_times <- data_ord %>%
    mutate(event = as.logical(event)) %>%
    mutate(tiempo = as.numeric(tiempo)) %>%
    filter(event == T) %>%
    select(tiempo) %>%
    as.vector() %>%
    unlist()

# Los tiempos deed 0 se ponen con 0.01
time.data[time.data == 0] <- 0.01
summary(time.data)

#---------------------------------------------------------------------
# Modelo Semiparamétrico

# Definimos el modelo
write("
data {
  int<lower=1> N;               // Número de observaciones
  int<lower=1> J;               // Número de intervalos
  int<lower=1> K;               // Número de covariables
  vector<lower=0>[N] time;      // Tiempos observados
  array[N] int<lower=0, upper =1> event;// Indicador de evento
  matrix[N, K] X;               // Matriz de diseño
  vector<lower=0>[J+1] cuts;    // Puntos de corte (debe incluir 0 y max(time))
}

transformed data {
  vector[J+1] sorted_cuts = sort_asc(cuts);  // Ordenar los cortes
}

parameters {
  vector[J] log_lambda;         // Log-riesgo base por intervalo
  vector[K] beta;               // Coeficientes de covariables
  cholesky_factor_corr[J] L;    // Factor de Cholesky para correlaciones
  real<lower=0> sigma;          // Desviación estándar de los log-lambdas
}

transformed parameters {
  vector[J] lambda = exp(log_lambda);
  matrix[J,J] Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(rep_vector(sigma,J), L));
}

model {
  // Priors
  log_lambda ~ multi_normal_cholesky(rep_vector(0,J), diag_pre_multiply(rep_vector(sigma,J), L));
  sigma ~ exponential(1);
  beta ~ normal(0, 3);
  L ~ lkj_corr_cholesky(2);     // Prior para matriz de correlación
  
  // Likelihood
  for (i in 1:N) {
    real lin_pred = X[i] * beta;
    real cum_haz = 0;
    int current_interval = 1;
    
    // Calcular riesgo acumulado hasta el tiempo t
    while (current_interval <= J && time[i] > sorted_cuts[current_interval + 1]) {
      real width = sorted_cuts[current_interval + 1] - sorted_cuts[current_interval];
      cum_haz += lambda[current_interval] * width;
      current_interval += 1;
    }
    
    // Último intervalo parcial
    if (current_interval <= J) {
      real width = time[i] - sorted_cuts[current_interval];
      cum_haz += lambda[current_interval] * width;
    }
    
    // Contribución a la verosimilitud
    if (event[i]) {
      target += log(lambda[current_interval]) + lin_pred - exp(lin_pred) * cum_haz;
    } else {
      target += -exp(lin_pred) * cum_haz;
    }
  }
}

generated quantities {
  corr_matrix[J] Omega = multiply_lower_tri_self_transpose(L);
  vector[N] log_lik;
  
  // Cálculo de log-verosimilitud punto a punto
  for (i in 1:N) {
    real lin_pred = X[i] * beta;
    real cum_haz = 0;
    int current_interval = 1;
    
    while (current_interval <= J && time[i] > sorted_cuts[current_interval + 1]) {
      real width = sorted_cuts[current_interval + 1] - sorted_cuts[current_interval];
      cum_haz += lambda[current_interval] * width;
      current_interval += 1;
    }
    
    if (current_interval <= J) {
      real width = time[i] - sorted_cuts[current_interval];
      cum_haz += lambda[current_interval] * width;
    }
    
    log_lik[i] = event[i] 
      ? log(lambda[current_interval]) + lin_pred - exp(lin_pred) * cum_haz
      : -exp(lin_pred) * cum_haz;
  }
}
      "
,"semiparametric_PH_model.stan")

# Compilamos el modelo
semiparametric_model <- cmdstan_model("semiparametric_PH_model.stan",
                                      pedantic = TRUE)
#---------------------------------------------------------------------------
# Preparación de datos para el modelo

# Definimos los cortes para el tiempo
J <- 9  # Número de intervalos deseado
cuts <- quantile(event_times, probs = seq(0, 1, length.out = J + 1)) %>% 
    unname()

# Asegurar que el último corte cubra el tiempo máximo (364)
cuts[length(cuts)] <- 364 + 1e-5  # Reemplazar el último cuantil

table(cut(event_times, breaks = cuts)) # Contamos los eventos

# Definimos los datos para el modelo
data_list <- list(
    N = nrow(data_ord),
    J = J,
    K = ncol(data_ord) - 2, # Número de covariables
    time = time.data,
    event = event.data,
    X = data_ord[, -c(16, 17)],
    cuts = cuts
)

# Ajustamos el modelo
fit <- semiparametric_model$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    refresh = 100,
    seed = 14082001
)
#---------------------------------------------------------------------------
# Guardamos el modelo
fit$save_object(file = "semiparametric_model_fit.rds")

#---------------------------------------------------------------------------
# Resumen de resultados
#fit <- readRDS("semiparametric_model_fit.rds")

sum_draws <- fit$summary(
    variables = c(
        "beta",        # Coeficientes de las covariables (efectos proporcionales)
        "lambda"      # Riesgo base por intervalo (transformed parameter)
        #"sigma",       # Desviación estándar de los log-lambdas
        #"log_lambda"  # Log-riesgo base por intervalo (parámetro original)
        #"Omega"        # Matriz de correlación completa (generated quantity)
    ),
    posterior::default_summary_measures(),
    default_convergence_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975))
)

write.csv(
    sum_draws,
    file = "semiparametric_model_summary.csv",
    row.names = FALSE
)

#----------------------------------------------------------------------------
# Correlaciones entre intervalos

cor_matrix <- fit$draws("Omega") %>% posterior::as_draws_matrix()

# Extrae los nombres de las variables (asume que Omega es KxK)
K <- 9
var_names <- paste0("V", 1:K)  # Cambia por nombres reales si los tienes

# Convertimos muestras de Omega (array 3D: iter x Omega[i,j]) a tibble
cor_samples <- cor_matrix %>%
    as_tibble(rownames = "iter") %>%
    pivot_longer(cols = -iter, names_to = "param", values_to = "value") %>%
    mutate(
        i = str_extract(param, "(?<=\\[)\\d+") %>% as.integer(),
        j = str_extract(param, "(?<=,)\\d+") %>% as.integer()
    )

# Calculamos estadísticas por par (i,j)
cor_summary <- cor_samples %>%
    group_by(i, j) %>%
    summarise(median = median(value),
              lower = quantile(value, 0.05),
              upper = quantile(value, 0.95),
              .groups = "drop")

# Mapa de calor con medianas y opcionalmente los intervalos
ggplot(cor_summary, aes(x = factor(j), y = factor(i), fill = median)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", median)), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "", y = "", fill = "Correlación") +
    coord_fixed() +
    theme_minimal()
#---------------------------------------------------------------------------
cor_samples %>%
    filter(i < j) %>%
    ggplot(aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    facet_grid(i ~ j, labeller = label_both) +
    theme_minimal() +
    labs(title = "Posterior de correlaciones", x = "Valor", y = "Densidad")

#---------------------------------------------------------------------------
# Betas del Modelo

draws_parameters <- fit$draws("beta", format = "df")

# Param Names
betas_str <- character(length =  ncol(data_ord) - 2)
for (i in 1:( ncol(data_ord) - 2)) {
    betas_str[i] <- paste("beta[", i, "]", sep = "")
}

betas_str

# Cadenas
p <- mcmc_trace(
    draws_parameters,
    regex_pars  = c("beta"),
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

mcmc_areas(draws_parameters, 
           regex_pars = c("beta"),
           prob_outer = 1, 
           point_est = "mean", 
           prob = 0.95,
           area_method = "equal height") +
    scale_y_discrete(labels = parse(text = betas_str)) +
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
    scale_y_discrete(labels = parse(text = betas_str)) + 
    scale_x_discrete(labels = parse(text = betas_str))
rm(cor_melt)
gc()
#---------------------------------------------------------------------------
# Base Hazard Rates
draws_lambda <- fit$draws("lambda", format = "df")

# Param Names
lamda_str <- character(length =  length(cuts) - 1)
for (i in 1:( length(cuts) - 1 )) {
    lamda_str[i] <- paste("lambda[", i, "]", sep = "")
}

lamda_str

# Cadenas
p <- mcmc_trace(
    draws_lambda,
    regex_pars  = c("lambda"),
    facet_args = list(nrow = 3, ncol = 3, labeller = label_parsed)
)
p + facet_text(size = 15)

# Gráficos Posteriores.
mcmc_areas(draws_lambda, 
           regex_pars = c("lambda"),
           prob_outer = 1, 
           point_est = "median", 
           prob = 0.95,
           area_method = "equal height") +
    scale_y_discrete(labels = parse(text = lamda_str)) +
    mcmc_areas_theme

# Calcular la matriz de correlación
cor_matrix <- cor(draws_lambda %>% select(-c(.chain, .iteration, .draw)))

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
    scale_y_discrete(labels = parse(text = lamda_str)) + 
    scale_x_discrete(labels = parse(text = lamda_str))
rm(cor_melt)
gc()


#------------------------------------------------------------------------
# Posterior Survival Function

plot_posterior_survival_semiparametric <- function(posterior_draws, covariates, cuts, t_max = 365, n_curves = 100) {
    # Extraer parámetros del modelo
    log_lambda <- as.matrix(posterior_draws %>% select(matches("^log_lambda\\[.*\\]")))
    betas <- as.matrix(posterior_draws %>% select(matches("^beta\\[.*\\]")))
    
    # Validaciones
    if (ncol(betas) != length(covariates)) {
        stop("Número de betas y covariables no coincide")
    }
    if (length(cuts) != ncol(log_lambda) + 1) {
        stop("Número de cortes no coincide con los intervalos del modelo")
    }
    
    # 1. Calcular lambda y predictor lineal
    lambda <- exp(log_lambda)
    linpred <- exp(betas %*% covariates)  # Hazard ratio: exp(Xβ)
    
    # 2. Secuencia de tiempo
    t_seq <- seq(0, t_max, length.out = 1000)
    J <- ncol(lambda)
    
    # 3. Calcular matriz de supervivencia
    S_matrix <- matrix(nrow = nrow(lambda), ncol = length(t_seq))
    
    for (i in 1:nrow(lambda)) {
        for (t_idx in 1:length(t_seq)) {
            t <- t_seq[t_idx]
            cum_hazard <- 0
            current_interval <- 1
            
            # Calcular hazard acumulado
            while (current_interval <= J && t > cuts[current_interval + 1]) {
                width <- cuts[current_interval + 1] - cuts[current_interval]
                cum_hazard <- cum_hazard + lambda[i, current_interval] * width
                current_interval <- current_interval + 1
            }
            
            if (current_interval <= J) {
                width <- t - cuts[current_interval]
                cum_hazard <- cum_hazard + lambda[i, current_interval] * width
            }
            
            # Aplicar efecto de covariables
            S_matrix[i, t_idx] <- exp(-cum_hazard * linpred[i])
        }
    }
    
    # 4. Estadísticos resumen
    median_S <- apply(S_matrix, 2, median)
    hdi_S <- apply(S_matrix, 2, function(x) HDInterval::hdi(x, credMass = 0.9))
    
    # 5. Dataframe para gráfico
    plot_df <- data.frame(
        time = t_seq,
        median = median_S,
        lower = hdi_S[1, ],
        upper = hdi_S[2, ]
    )
    
    # 6. Muestrear curvas
    set.seed(123)
    idx <- sample(nrow(S_matrix), min(n_curves, nrow(S_matrix)))
    posterior_curves <- data.frame(
        time = rep(t_seq, each = length(idx)),
        survival = as.vector(S_matrix[idx, ]),
        draw = rep(1:length(idx), times = length(t_seq))
    )
    
    # 7. Graficar
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
            title = "Curvas de Supervivencia Posteriores (Modelo Semi-Paramétrico)",
            x = "Tiempo",
            y = "Probabilidad de Supervivencia"
        ) +
        ylim(0, 1) +
        geom_vline(xintercept = cuts, linetype = "dashed", alpha = 0.3) + # Líneas de corte
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

# Graficamos la función de supervivencia posterior
plot_posterior_survival_semiparametric(
    posterior_draws = fit$draws(format = "df"),
    covariates = covariables,
    cuts = cuts,
    t_max = 365,
    n_curves = 100
)
gc()
