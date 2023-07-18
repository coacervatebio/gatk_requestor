#!/bin/zsh

REPO=$1
IMAGE=$2
TAG=$3

printf "Building Docker image: %s/%s:%s\n" "$REPO" "$IMAGE" "$TAG"
docker build . -t "$REPO/$IMAGE:$TAG"

printf "Saving Docker image: %s/%s:%s\n" "$REPO" "$IMAGE" "$TAG"
docker save "$REPO/$IMAGE:$TAG" -o /tmp/"$REPO-$IMAGE-$TAG.tar"

printf "Copying Docker image to flyte-sandbox:/tmp\n"
docker cp /tmp/"$REPO-$IMAGE-$TAG.tar" flyte-sandbox:/tmp

printf "Importing Docker image on flyte-sandbox\n"
docker exec -t flyte-sandbox ctr image import /tmp/"$REPO-$IMAGE-$TAG.tar"
