FROM continuumio/miniconda3

ADD env.yaml /tmp/env.yaml

RUN conda env create -f /tmp/env.yaml -n ipd-ref-generator python=3.8
RUN echo "conda activate ipd-ref-generator" > ~/.bashrc
ENV PATH /opt/conda/envs/ipd-ref-generator/bin:$PATH
