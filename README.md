# Coacervate

## Introduction

Genomic data has been generated [faster and faster](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002195) since the first human genome was sequenced - with no signs of slowing down. The demand for expertise, and the raw compute power needed to turn those enormous datasets into actionable insights, has exploded alongside it.

Analyzing that data relies on a highly-specialized, [multi-disciplinary](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471940/) skillset. Tools exist that make this step easier, but they are often black-boxes and prohibitively expensive at scale; this precludes access from aspiring researchers and stifles collaboration.

Usability aside, there's no getting around the sheer scale of compute resources needed. Coupled with the fact that we can no longer rely on [Moore's Law](https://www.technologyreview.com/2016/05/13/245938/moores-law-is-dead-now-what/) to pick up the slack, we need to get creative to maximize the pace of progress.

Leveraging the fact that people are more likely to have access to [high-bandwidth connectivity](https://www.fiercetelecom.com/broadband/att-upgrades-its-fiber-network-offer-2-gig-5-gig-speeds) than to high-performance compute, Coacervate lets you run genomic analyses on an extremely low-cost [global supercomputer](https://www.golem.network). By using incredible open-source [tools](https://gatk.broadinstitute.org/hc/en-us) and [frameworks](https://flyte.org/), this project aims to democratize access to the knowledge _and_ infrastructure required to carry out groundbreaking research.

Coacervate is a free and open-source public good, built to empower every citizen-scientist and eek out every last drop of efficiency in the name of progess.

## Approach
Coacervate accepts genomic sequence alignments and produces [joint called](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants) VCF files ready for annotation and actionable interpretation. This is achieved by splitting the inputs to parallelize the most [computationally intensive](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0254363#pone-0254363-g003) step via the Golem Network. 

## Status
⚠️⚠️⚠️ This project is currently a **Proof of Concept** and undergoing a *significant* refactor. The Getting Started tutorial and testing are very much under-construction.

## Getting Started
Please follow the [Getting Started tutorial](getting_started.md) to run your first workflow! Broadly, you'll need to: 
- make sure docker and git are installed
- install flytectl and start the demo cluster
- git pull the project
- register the workflows

## Getting Started with Your Own Data
Coacervate relies on Flyte to run it's workflows, which in turn uses an object store to manage all the data being worked on. If you're following the tutorial then you'll just use the small alignments baked into the requestor image. However, if you'd like to use your own data then you'll need to upload it to the object store and specify that upload path in the launchplan when firing your workflow. The flyte sandbox will spin up a minio instance which you can conveniently access in your browser at `http://localhost:30080/minio/browser` (you can login in with `minio`/`miniostorage`). Flyte can obviously use any object store as a backend, like S3, but that setup is out of scope for this README. Please also note that currently inputs must be in the CRAM format. 

## Repo Tour
Most of the layout of this repo is fairly self-explanatory, but a brief tour can be helpful for those that want to adapt functionality to their needs. Please refer to a slightly redacted tree of salient parts below.
- The `reference` and `data` directories contain mostly static assets. Outputs from different steps are typically written to `data/output`.
- The `run` dir contains the config, tasks, and workflow that make up the GATK proof-of-concept.
   - Golem requestor agents are contained in `agents`, `hello_golem.py` exists as a simple test agent.
   - Advanced pod configuration for communicating with the yagna daemon deployed in the cluster hosting Flyte is captured in `yagna_template.py`.
   - Different Flyte tasks are separated by alignment mapping, variant calling and supporting utilities.
   - Tests covering these tasks are also available with the testing strategy described in [this gist](https://gist.github.com/pryce-turner/298ab8bb7f8bb7ee1b2507dc068c938e).
   - Finally, the POC workflow is available in it's eponymous directory, additional workflows should live alongside it.
```
src
├── requirements.txt
├── reference
│   └── resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
├── data
│   ├── alignments
│   │   ├── HG03633.cram
│   │   └── HG04149.cram
│   └── output
└── run
    ├── __init__.py
    ├── config.yaml
    ├── agents
    │   ├── haplotypecaller.py
    │   ├── hello_golem.py
    ├── pod
    │   └── yagna_template.py
    ├── tasks
    │   ├── calling.py
    │   ├── mapping.py
    │   └── utils.py
    ├── tests
    │   ├── helpers.py
    │   ├── register.py
    │   ├── test_combine_region.py
    │   ├── test_gather_vcfs.py
    │   ├── test_genotype.py
    │   ├── test_golem_call_variants.py
    │   ├── test_index_cram.py
    │   └── test_split_cram.py
    └── workflows
        ├── gatk_best_practices.py
```

## Future Plans
- Distributed data storage in IPFS to remove 1-to-many data transfer bottleneck
- Fully fog-based, only access a portal from your browser
- FHE for clinical samples / sensitive data using [Zama](https://www.zama.ai/)
- Verifiable compute using [Risc Zero](https://www.risczero.com/)

### References
  - Stephens, Zachary D., et al. “Big Data: Astronomical or Genomical?” PLOS Biology, Public Library of Science, https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.1002195. 
  - Heo, Go Eun, et al. “Analyzing the Field of Bioinformatics with the Multi-Faceted Topic Modeling Technique.” BMC Bioinformatics, BioMed Central, 31 May 2017, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471940/.
  - John, Aji, et al. “Evaluation of Serverless Computing for Scalable Execution of a Joint Variant Calling Workflow.” PLOS ONE, Public Library of Science, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0254363.
  - Nature Publishing Group. (2015, September 30). A global reference for human genetic variation. Nature News. https://www.nature.com/articles/nature15393 
