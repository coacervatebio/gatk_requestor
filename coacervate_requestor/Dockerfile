FROM condaforge/mambaforge:4.13.0-1

RUN apt-get update && apt-get install --no-install-recommends -y screen jq 
RUN useradd -u 1000 -s /bin/bash --create-home requestor

COPY --chown=1000:1000 config /data/config
COPY --chown=1000:1000 resources /data/resources
COPY --chown=1000:1000 results /data/results
COPY --chown=1000:1000 workflow /data/workflow

# Update base env, since that's where we get dropped into by default
# This installs GATK, snakemake and yapapi
RUN mamba env update --name base --file /data/workflow/envs/env.yml

# Install yagna requestor
WORKDIR /home/requestor
ARG YAG_VER=v0.10.1
ADD "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" .
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*

COPY --chown=1000:1000 start.sh /home/requestor/start.sh
USER 1000
ENTRYPOINT ["/bin/bash", "/home/requestor/start.sh"]
CMD [ "-m", "golem" ]
