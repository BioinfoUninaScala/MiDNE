FROM rocker/shiny:latest

LABEL maintainer="Giovanni Scala giovanni.scala@unina.it"

# system libraries of general use
RUN apt-get update && apt-get install -y libx11-dev make zlib1g-dev libglpk-dev libxml2-dev cmake libfreetype6-dev libjpeg-dev libpng-dev libtiff-dev libcurl4-openssl-dev pandoc libfontconfig1-dev libfribidi-dev libharfbuzz-dev


# R packages

RUN R -e "install.packages(c('bigmemory', 'bigstatsr', 'crosstalk', 'data.table', 'dendextend', 'doFuture','doMC', 'doParallel', 'doSNOW', 'DT', 'dplyr', 'foreach', 'fpc', 'ggplot2', 'glue','graphics', 'igraph', 'kernlab', 'Matrix', 'methods', 'optparse', 'parallel','plotly', 'purrr', 'R.matlab', 'RColorBrewer', 'readr', 'reshape2', 'rstatix','Rtsne', 'shiny', 'shinyalert', 'shinycssloaders', 'shinydashboard', 'shinyjs','shinyMatrix', 'snow', 'stats', 'stringr', 'tibble', 'tidyverse', 'utils', 'uwot','visNetwork', 'vroom', 'zip', 'BiocManager', 'remotes'))"

RUN R -e "BiocManager::install('gprofiler2')"

RUN R -e "remotes::install_github('BioinfoUninaScala/MiDNE', build_vignettes=FALSE, repos=BiocManager::repositories(),dependencies=TRUE, type='source')"


COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(MiDNE::MiDNEshiny(), host = '0.0.0.0', port = 3838)"]