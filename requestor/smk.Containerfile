FROM condaforge/mambaforge

VOLUME /yagna

COPY benchmarks /mnt/benchmarks
COPY config /mnt/config
# COPY resources /mnt/resources
COPY results /mnt/results
COPY workflow /mnt/workflow

# Install yagna requestor
ARG YAG_VER=v0.10.1
WORKDIR /requestor
RUN wget "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz"
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/bin/
RUN rm -rf /requestor/*

# Update base env, since that's where we get dropped into by default
RUN mamba env update --name base --file /mnt/workflow/envs/env.yml

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n=base"]
CMD ["snakemake", "-c=1", "-s=/mnt/workflow/rules/Snakefile"]