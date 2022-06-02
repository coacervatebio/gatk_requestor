#!/bin/zsh
echo "Building.."
podman build -f Containerfile -t coacervate:testing

echo ""
echo "Running new image.."
podman run --rm -it localhost/coacervate:testing snakemake --help

echo "Cleaning up.."
podman rmi localhost/coacervate:testing