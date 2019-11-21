FROM continuumio/anaconda

RUN cd /home && \
    git clone https://github.com/jason-weirather/hla-polysolver.git && \
    cp /home/hla-polysolver/.condarc /root/.condarc && \
    conda install -c vacation hla-polysolver && \
    rm -r /home/hla-polysolver

RUN apt-get update && apt-get install -y \
    libncurses5 \
    && apt-get clean && apt-get purge

ENV CONDA_PREFIX /opt/conda

WORKDIR /home
