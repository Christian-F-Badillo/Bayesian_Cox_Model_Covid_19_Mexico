# Usar Python 3.11 como base
FROM python:3.11

# Instalar herramientas de compilación necesarias
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    python3-dev \
    libatlas-base-dev \
    libhdf5-dev \
    gfortran \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev

# Actualizar pip
RUN pip install --upgrade pip

# Instalar PyMC y dependencias necesarias
RUN pip install pymc3 numpy pandas matplotlib seaborn arviz[all] jupyter scikit-survival rpy2

# Instalar paquetes de R esenciales
RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'data.table', 'rstan', 'lme4'), repos='http://cran.r-project.org')"

# Establecer el directorio de trabajo
WORKDIR /app

# Copiar archivos del proyecto
COPY . .

# Exponer el puerto de Jupyter Notebook
EXPOSE 8888

# Comando por defecto: iniciar Jupyter Notebook
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
