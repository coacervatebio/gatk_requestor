FROM python:3.9-slim-bullseye

WORKDIR /root
ENV VENV /opt/venv
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONPATH /root

# Apt installs
ARG JRE_VER=17
RUN apt-get update && apt-get install -y build-essential samtools openjdk-${JRE_VER}-jre unzip

# Install the AWS cli separately to prevent issues with boto being written over
RUN pip3 install awscli

# Similarly, if you're using GCP be sure to update this command to install gsutil
# RUN apt-get install -y curl
# RUN curl -sSL https://sdk.cloud.google.com | bash
# ENV PATH="$PATH:/root/google-cloud-sdk/bin"

ENV VENV /opt/venv
# Virtual environment
RUN python3 -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

# Install Python dependencies
COPY ./requirements.txt /root
RUN pip install -r /root/requirements.txt

# Copy the actual code
COPY . /root

# This tag is supplied by the build script and will be used to determine the version
# when registering tasks, workflows, and launch plans
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag

# Install yagna requestor
ARG YAG_VER=v0.12.0
ADD "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" .
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*


# Install GATK4
WORKDIR /home
ARG GATK_VER=4.4.0.0
ADD https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip .
RUN unzip gatk-${GATK_VER}.zip && mv gatk-${GATK_VER}/gatk-package-${GATK_VER}-local.jar /usr/local/share/ && rm -rf gatk-${GATK_VER}*
