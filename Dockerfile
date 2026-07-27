FROM rocker/r-ver:4.3.0
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libgit2-dev pandoc \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install --break-system-packages \
    pandas==2.1.4 numpy==1.26.4 scipy==1.11.4 \
    scikit-learn==1.4.2 statsmodels==0.14.1 \
    matplotlib==3.8.4 seaborn==0.13.2 \
    requests==2.31.0 tqdm==4.66.2 openpyxl==3.1.2
RUN install2.r --error --skipinstalled \
    data.table pROC ggplot2 matrixStats knitr rmarkdown \
    testthat ComplexHeatmap survminer survival
WORKDIR /workspace
COPY . /workspace/paralog_sl_predictor/
WORKDIR /workspace/paralog_sl_predictor
RUN R -e 'install.packages("R_package/", repos = NULL, type = "source")'
CMD ["R", "--no-save"]
