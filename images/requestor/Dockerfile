FROM condaforge/mambaforge:4.13.0-1

RUN useradd -u 1000 -s /bin/bash --create-home coacervate

COPY --chown=1000:1000 config /data/config
COPY --chown=1000:1000 resources /data/resources
COPY --chown=1000:1000 results /data/results
COPY --chown=1000:1000 workflow /data/workflow

# Update base env, since that's where we get dropped into by default
# This installs GATK, snakemake and yapapi
RUN mamba env update --name base --file /data/workflow/envs/env.yml

CMD [ "-m", "golem", "-d", "on" ]
