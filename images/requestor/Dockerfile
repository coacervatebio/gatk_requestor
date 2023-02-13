FROM condaforge/mambaforge:4.13.0-1

COPY config /data/config
COPY resources /data/resources
COPY results /data/results
COPY workflow /data/workflow

WORKDIR /home
ARG YAG_VER=v0.12.0
ADD "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" .
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*

# Update base env, since that's where we get dropped into by default
# This installs GATK, snakemake and yapapi
RUN mamba env update --name base --file /data/workflow/envs/env.yml

CMD [ "snakemake", "-f", "/data/workflow/Snakefile", "-c", "4" ]
