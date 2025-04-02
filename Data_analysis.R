## Author: Christian Badillo
## Fecha: 2025-04-02
## Modelo de Riesgos proporcionales para datos de Covid-19 en México en el año del 2022 usando
## estimación bayesiana.

# librerias
library(rstan)
library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# Setup

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#rstan_options(threads_per_chain = 1)

# --------------------------------------------------------------
# Preparación de Datos

data <- read.csv("COVID19MEXICO2022_muestra_clean.csv")

View(data)

names(data)

summary(data)

data$event <- ifelse(data$event == "True", 1, 0)

# --------------------------------------------------------------
# Vars to Factors
data$tipo_paciente <- factor(data$tipo_paciente, labels = c("Ambulatorio", "Hospitalizado"))
data$sexo <- factor(data$sexo, labels = c("Hombre", "Mujer"))

data <- data %>% select(-c(nacionalidad))

data <- data %>% 
    mutate(across(c("neumonia", "indigena", "diabetes", 
                  "epoc", "asma", "inmusupr", "hipertension", 
                  "otra_com", "cardiovascular", "obesidad", 
                  "renal_cronica", "tabaquismo", "otro_caso"),
                  ~ factor(.x, labels = c("Si", "No"))))

# --------------------------------------------------------------
# Estimación de la curva de supervivencia usando Kaplan-Meier

sur_obj <- Surv(time = data$tiempo, event = data$event)

## General Kaplan Meier
sur_obj <- Surv(time = data$tiempo, event = data$event)

km_fit <- survfit(sur_obj ~ 1)

ggsurvplot(km_fit, data = data, conf.int = T, ylim=c(0.80, 1), legend = "none") + 
    labs(
        title = "Curva de sobrevivencia para Covid-19 en México en el año 2022",
        x = "Tiempo (días)",
        y = "Probabilidad de Sobrevivencia")

## KM por demás variables

KM_vars <- function(vars, data, sur_obj) {
    for (var in vars) {
        # Crear la fórmula correctamente
        form <<- formula(paste0(deparse(substitute(sur_obj)), " ~ ", var))
        
        # Ajustar el modelo de Kaplan-Meier
        km_fit_var <- survfit(as.formula(form), data = data)
        
        # Graficar la curva de sobrevivencia y asegurar que se imprime
        print(
            ggsurvplot(km_fit_var, data = data, conf.int = TRUE, ylim = c(0, 1)) +
                labs(
                    title = paste("Curva de sobrevivencia para Covid-19 por", var),
                    x = "Tiempo (días)",
                    y = "Probabilidad de Sobrevivencia"
                )
        )
    }
}


vars_to_KM <- c("sexo", "tipo_paciente" , "neumonia", "diabetes", "epoc", "hipertension", 
                "cardiovascular", "obesidad", "tabaquismo", "asma", 
                "inmusupr", "otra_com", "renal_cronica", "otro_caso")

KM_vars(vars_to_KM, data, sur_obj)
