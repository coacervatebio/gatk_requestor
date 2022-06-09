FROM ubuntu
RUN apt update && apt install -y python3 wget pip git
RUN pip install yapapi
WORKDIR /golem

ARG YAG_VER=v0.10.1

RUN wget "https://github.com/golemfactory/yagna/releases/download/${YAG_VER}/golem-requestor-linux-${YAG_VER}.tar.gz"
RUN tar -xzf "golem-requestor-linux-${YAG_VER}.tar.gz"
RUN mv golem-requestor-linux-${YAG_VER}/* /usr/bin/
RUN rm -rf /golem/*

COPY ./start.sh /start.sh
CMD ["/bin/bash", "-c", "/start.sh"]