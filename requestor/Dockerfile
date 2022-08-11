FROM condaforge/mambaforge

RUN apt-get update && apt-get install screen jq -y

VOLUME /yagna

COPY benchmarks /mnt/benchmarks
COPY config /mnt/config
COPY resources /mnt/resources
COPY results /mnt/results

# Install yagna requestor
ARG YAG_VER=v0.10.1
WORKDIR /requestor
RUN wget "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz"
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/bin/
RUN rm -rf /requestor/*

COPY workflow /mnt/workflow
# Update base env, since that's where we get dropped into by default
# This installs GATK, snakemake and yapapi
RUN mamba env update --name base --file /mnt/workflow/envs/env.yml

COPY start.sh /requestor/start.sh
ENTRYPOINT ["/bin/bash", "/requestor/start.sh"]
CMD [ "-m", "golem" ]
