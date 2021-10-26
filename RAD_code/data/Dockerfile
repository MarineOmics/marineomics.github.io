# Base image -  can use R, jupyter notebooks
FROM jupyter/r-notebook
USER root 

# Install modules needed
RUN apt-get update \
 && apt-get install -y --no-install-recommends vcftools \
 && rm -rf /var/lib/apt/lists/*
 
# Install the ggplot2 library and a few other dependencies we want to have available
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('pheatmap')"
RUN Rscript -e "install.packages('geosphere')"
RUN Rscript -e "install.packages('grImport2')"
RUN Rscript -e "install.packages('dendextend')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("SNPRelate");' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("SeqArray")'


CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
