FROM ubuntu
RUN apt update && apt install -y python3 wget pip git
RUN pip install yapapi
WORKDIR /requestor

COPY ./start.sh /start.sh
CMD ["/bin/bash", "-c", "/start.sh"]