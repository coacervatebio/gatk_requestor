FROM python:3.9-bullseye

# Install yagna requestor
ARG YAG_VER=v0.12.0
ADD "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" .
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*

# Apt installs
ARG JRE_VER=17
RUN apt-get update && apt-get install -y samtools openjdk-${JRE_VER}-jre

# Install GATK4
WORKDIR /home
ARG GATK_VER=4.4.0.0
ADD https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip .
RUN unzip gatk-${GATK_VER}.zip && mv gatk-${GATK_VER}/gatk-package-${GATK_VER}-local.jar /usr/local/share/ && rm -rf gatk-${GATK_VER}*

# Pip installs
RUN pip install yapapi flytekit

# Copy workflow dependencies
COPY config /data/config
COPY reference /data/reference
COPY results /data/results
COPY agents /data/agents
