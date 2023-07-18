#!/bin/zsh

REGISTRY=$1
REPO=$2
IMAGE=$3
TAG=$4

printf "Building Docker image: %s/%s/%s:%s\n" "$REGISTRY" "$REPO" "$IMAGE" "$TAG"
docker build . -t "$REGISTRY/$REPO/$IMAGE:$TAG"

printf "Saving Docker image: %s/%s/%s:%s\n" "$REGISTRY" "$REPO" "$IMAGE" "$TAG"
docker save "$REGISTRY/$REPO/$IMAGE:$TAG" -o /tmp/"$REGISTRY-$REPO-$IMAGE-$TAG.tar"

printf "Copying Docker image to flyte-sandbox:/tmp\n"
docker cp /tmp/"$REGISTRY-$REPO-$IMAGE-$TAG.tar" flyte-sandbox:/tmp

printf "Importing Docker image on flyte-sandbox\n"
docker exec -t flyte-sandbox ctr image import /tmp/"$REGISTRY-$REPO-$IMAGE-$TAG.tar"
