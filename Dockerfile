FROM docker.io/rocker/rstudio:4.5.3 
LABEL Author Philip Smith
LABEL Version v1.2.0
LABEL org.opencontainers.image.title="CINSignatureQuantification"
LABEL org.opencontainers.image.description="CINSignatureQuantification package NON-COMMERCIAL USE ONLY."
LABEL org.opencontainers.image.licenses="LicenseRef-GAP-ASL-1.0"
LABEL org.opencontainers.image.source="https://github.com/markowetzlab/CINSignatureQuantification"
LABEL org.opencontainers.image.documentation="https://github.com/markowetzlab/CINSignatureQuantification"

ENV PATH=$PATH:/usr/local/lib/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

RUN apt-get -y update
RUN apt-get install -y unzip less libcurl4-openssl-dev git-all curl
RUN apt-get clean

COPY LICENSE /licenses/CINSignatureQuantification-LICENSE.txt
COPY THIRD_PARTY_NOTICES.txt /licenses/THIRD_PARTY_NOTICES.txt

RUN Rscript -e 'install.packages("pak", repos="https://cloud.r-project.org",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'BiocManager::install("Biobase",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'BiocManager::install("QDNAseq",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'install.packages("data.table",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'install.packages("limSolve",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'install.packages("stringr",lib="/usr/local/lib/R/site-library/")'
RUN Rscript -e 'pak::pkg_install("markowetzlab/CINSignatureQuantification",lib="/usr/local/lib/R/site-library/",dependencies=TRUE)'

RUN mkdir -p /licenses && \
    Rscript -e 'file.copy(system.file("LICENSE", package = "stringi"), "/licenses/stringi-LICENSE.txt")'

RUN echo "Rscript -e 'library(CINSignatureQuantification)'" >> /tests.sh
RUN echo "echo $LD_LIBRARY_PATH" >> /tests.sh
RUN echo "echo $PATH" >> /tests.sh
RUN chmod u+x /tests.sh

RUN echo '#!/bin/bash\n\
cat << "EOF" >&2\n\
================================================================\n\
NOTICE: This image includes CINSignatureQuantification, licensed\n\
under the GAP Available Source License v1.0 (ASL).\n\
NON-COMMERCIAL ACADEMIC USE ONLY.\n\
Full license: /licenses/CINSignatureQuantification-LICENSE.txt\n\
Repo: https://github.com/markowetzlab/CINSignatureQuantification\n\
================================================================\n\
EOF\n\
exec "$@"\n'\ >> /entrypoint.sh

run chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

