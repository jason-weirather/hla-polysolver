FROM continuumio/anaconda

RUN cd /home && \
    git clone https://github.com/jason-weirather/hla-polysolver.git && \
    cp /home/hla-polysolver/.condarc /root/.condarc && \
    conda install -c vacation hla-polysolver && \
    rm -r /home/hla-polysolver

RUN echo "export CONDA_PREFIX=/opt/conda" >> /root/.bashrc

WORKDIR /home
