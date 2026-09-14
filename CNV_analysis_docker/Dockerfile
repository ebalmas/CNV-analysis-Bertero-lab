# RStudio Server + R 4.5.2, basato su Ubuntu 24.04 LTS
FROM rocker/rstudio:4.5.2

ENV DEBIAN_FRONTEND=noninteractive

RUN install2.r \
      --error \
      --skipinstalled \
      --ncpus -1 \
      dplyr \
      xml2 \
      ggplot2 \
      viridis \
rmarkdown \
knitr \
    && rm -rf /tmp/downloaded_packages

COPY --chown=rstudio:rstudio \
    script/ \
    /home/rstudio/project/script/

COPY --chown=rstudio:rstudio \
    examples/ \
    /home/rstudio/project/examples/

COPY --chown=rstudio:rstudio \
    README.md \
    run_CNV_analysis.Rmd \
    run_CNV_analysis_multi_XML.Rmd \
    /home/rstudio/project/

COPY --chown=rstudio:rstudio \
    SOP/ \
    /home/rstudio/project/SOP/

RUN mkdir -p \
      /home/rstudio/project/data \
      /home/rstudio/project/Output \
      /home/rstudio/project/scratch/input_files \
    && chown -R rstudio:rstudio /home/rstudio/project

WORKDIR /home/rstudio/project

EXPOSE 8787
CMD ["/init"]