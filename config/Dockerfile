FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD env.yaml /tmp/env.yaml

# Install Mamba
RUN conda install mamba libarchive -c conda-forge --yes

RUN mamba env create -f /tmp/env.yaml -n ipd-ref-generator \
	&& mamba clean --all -f --yes \
	&& mamba init bash \
	&& echo "conda activate ipd-ref-generator" > ~/.bashrc

ENV PATH=/opt/conda/envs/ipd-ref-generator/bin:$PATH

ENV CONDA_DEFAULT_ENV=ipd-ref-generator

COPY src/downloadIPD.go .

RUN go mod init goDownloadIPD && \
	go get github.com/gofrs/flock && \
	go build -o goDownloadIPD downloadIPD.go && \
	mv goDownloadIPD /usr/local/bin && \
	chmod +x /usr/local/bin/goDownloadIPD

CMD ["/bin/bash"]
