# Librerias

## Bayes libreria
library(cmdstanr)
library(posterior)
library(bridgesampling)
options(mc.cores = parallel::detectCores())

# Gráficos
library(repr)
library(bayesplot)
library(ggplot2)
color_scheme_set("teal")


## Menjo de Datos
library(dplyr)
library(tidyr)
library(reshape2)