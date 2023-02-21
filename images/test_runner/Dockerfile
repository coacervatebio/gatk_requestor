FROM debian:stable-slim

# Install python/pip
RUN apt-get update && apt-get install -y python3 python3-pip

# Install testing dependencies
RUN pip install pysam click

# Add checking logic
COPY checker /checker