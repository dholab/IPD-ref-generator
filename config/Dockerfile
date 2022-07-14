FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD env.yaml /tmp/env.yaml

RUN conda env create -f /tmp/env.yaml -n ipd-ref-generator python=3.8 \
	&& conda clean --all -f --yes \
	&& conda init bash \
	&& echo "conda activate ipd-ref-generator" > ~/.bashrc \
	&& /bin/bash -c 'source ~/.bashrc'

ENV PATH /opt/conda/envs/ipd-ref-generator/bin:$PATH

ENV CONDA_DEFAULT_ENV ipd-ref-generator

CMD ["/bin/bash"]
