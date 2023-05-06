#!/bin/zsh

docker build . -t coacervate/requestor:latest
docker save coacervate/requestor:latest -o ~/Downloads/coacervate-requestor-latest.tar
docker cp ~/Downloads/coacervate-requestor-latest.tar flyte-sandbox:/tmp
docker exec -t flyte-sandbox ctr image import /tmp/coacervate-requestor-latest.tar