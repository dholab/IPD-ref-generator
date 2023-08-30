FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD env.yaml /tmp/env.yaml

RUN conda env create -f /tmp/env.yaml -n ipd-ref-generator \
	&& conda clean --all -f --yes \
	&& conda init bash \
	&& echo "conda activate ipd-ref-generator" > ~/.bashrc

ENV PATH=/opt/conda/envs/ipd-ref-generator/bin:$PATH

ENV CONDA_DEFAULT_ENV=ipd-ref-generator

COPY src/downloadIPD.go .

RUN go build -o goDownloadIPD downloadIPD.go && \
	mv goDownloadIPD /usr/local/bin && \
	chmod +x /usr/local/bin/goDownloadIPD

CMD ["/bin/bash"]
