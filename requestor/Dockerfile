FROM condaforge/mambaforge

RUN apt-get update && apt-get install screen jq -y
RUN useradd -u 5000 -s /bin/bash --create-home requestor

COPY config /data/config
COPY resources /data/resources
COPY results /data/results
COPY workflow /data/workflow
RUN chown -R 5000:5000 /data

# Install yagna requestor
WORKDIR /home/requestor
ARG YAG_VER=v0.10.1
RUN wget "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz"
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*

# Update base env, since that's where we get dropped into by default
# This installs GATK, snakemake and yapapi
RUN mamba env update --name base --file /data/workflow/envs/env.yml

USER 5000
COPY start.sh /home/requestor/start.sh
ENTRYPOINT ["/bin/bash", "/home/requestor/start.sh"]
CMD [ "-m", "golem" ]
