FROM condaforge/mambaforge

# Create the environment:
WORKDIR /env
COPY snakemake/workflow/envs/env.yml .
RUN mamba env create -f env.yml

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n=coacervate_env"]
CMD ["snakemake", "-c=1", "-s=/config/Snakefile"]