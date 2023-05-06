FROM python:3.9-slim-bullseye

WORKDIR /root
ENV VENV /opt/venv
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONPATH /root

# Apt installs
ARG JRE_VER=17
RUN apt-get update && apt-get install --no-install-recommends -y \
 build-essential openjdk-${JRE_VER}-jre unzip \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 wget \
 make \
 perl \
 bzip2 \
 gnuplot \
 ca-certificates \
 gawk \
 python3 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install the AWS cli separately to prevent issues with boto being written over
RUN pip3 install awscli

ENV VENV /opt/venv
# Virtual environment
RUN python3 -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

# This tag is supplied by the build script and will be used to determine the version
# when registering tasks, workflows, and launch plans
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag

# Install yagna requestor
ARG YAG_VER=v0.12.0
RUN wget "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" && \
 tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz" && \
 mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/ && \
 rm -rf golem-requestor*

# Install GATK4
ARG GATK_VER=4.4.0.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip && \
 unzip gatk-${GATK_VER}.zip && \
 mv gatk-${GATK_VER}/gatk-package-${GATK_VER}-local.jar /usr/local/share/ && \
 rm -rf gatk-${GATK_VER}*

# Install samtools
ARG SAMTOOLS_VER=1.17
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
 rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
 cd samtools-${SAMTOOLS_VER} && \
 ./configure --prefix=/usr/local && \
 make && \
 make install

# Install Python dependencies
COPY ./requirements.txt /root
RUN pip install -r /root/requirements.txt

# Copy the actual code
COPY . /root