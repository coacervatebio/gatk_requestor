FROM alpine:latest

# Install yagna requestor
WORKDIR /home/coacervate
ARG YAG_VER=v0.12.0
ADD "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz" .
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/local/bin/
RUN rm -rf golem-requestor*

COPY --chown=1000:1000 start.sh /home/coacervate/start.sh
RUN mkdir -p /home/coacervate/.local/share/yagna && chown -R 1000:1000 /home/coacervate/.local/
USER 1000
ENTRYPOINT ["/bin/bash", "/home/coacervate/start.sh"]