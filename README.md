# Coacervate

## Introduction

Genomic data has been generated [faster and faster](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002195) since the first human genome was sequenced - with no signs of slowing down. The demand for expertise, and the raw compute power needed to turn those enormous datasets into actionable insights, has exploded alongside it.

Analyzing that data relies on a highly-specialized, [multi-disciplinary](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471940/) skillset. Tools exist that make this step easier, but they are often black-boxes and prohibitively expensive at scale - precluding access from aspiring researchers and stifling collaboration.

Usability aside, there's no getting around the sheer scale of computing resources needed. Coupled with the fact that we can no longer rely on [Moore's Law](https://www.technologyreview.com/2016/05/13/245938/moores-law-is-dead-now-what/) to pick up the slack, we need to get creative to maximize the pace of progress.

Leveraging an extremely low-cost [global supercomputer](golem.network) and incredible open-source [tools](https://gatk.broadinstitute.org/hc/en-us) and [frameworks](https://snakemake.github.io/) - this project aims to democratize access to the knowledge and infrastructure required to carry out these analyses. Coacervate is a free and open-source public good, built to empower every citizen-scientist and eek out every last drop of efficiency in the name of biomedical progess.


## Quickstart:
- `docker run --rm -it -v yagna_datadir:/home/coacervate/.local/share/yagna coacervate/requestor`

  This will pull and run the requestor component of Coacervate. A Yagna daemon is spun up in the container and the default workflow is executed on test data bundled with the image.

## Architecture
-

## Features:
- option to mount yagna datadir at /yagna to persist app key and funding across runs
- option to mount any of your own dirs (config, results, snakefile) under /mnt
  - TODO: add hierarchy here
- known issues

