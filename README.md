Tractometry pipeline
====================

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/scilus/tractometry_flow)](https://github.com/scilus/tractometry_flow/releases)

[![Nextflow](https://img.shields.io/badge/nextflow-21.10.6-brightgreen.svg)](https://www.nextflow.io/)
[![Docker container badge](https://img.shields.io/docker/v/scilus/scilus?label=docker&logo=docker&logoColor=white)](https://hub.docker.com/r/scilus/scilus)

This pipeline allows you to extract tractometry information by combining
subjects's fiber bundles and diffusion MRI metrics.

Should you use this pipeline for your research, **please cite the following**

```
Cousineau, M., P-M. Jodoin, E. Garyfallidis, M-A. Cote, F.C. Morency, V. Rozanski, M. Grand'Maison, B.J. Bedell, and M. Descoteaux.
"A test-retest study on Parkinson's PPMI dataset yields statistically significant white matter fascicles."
NeuroImage: Clinical 16, 222-233 (2017) doi:10.1016/j.nicl.2017.07.020

Kurtzer GM, Sochat V, Bauer MW Singularity: Scientific containers for
mobility of compute. PLoS ONE 12(5): e0177459 (2017)
https://doi.org/10.1371/journal.pone.0177459

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows.
Nature Biotechnology 35, 316–319 (2017) https://doi.org/10.1038/nbt.3820
```

Requirements
------------

- [Nextflow](https://www.nextflow.io)

Please download [nextflow_21.10.6](https://github.com/nextflow-io/nextflow/releases/download/v21.10.6/nextflow-21.10.6-all) (or any version <= 21.10.6).

Installation
-----------

`ǹextflow pull scilus/tractometry_flow`

Usage
-----------

See *USAGE* or run `nextflow run scilus/tractometry_flow -r 1.2.0 --help`

Singularity/Docker
-----------
If you are on Linux, we recommend using the Singularity to run tractometry_flow pipeline.
If you have Apptainer (Singularity), launch your Nextflow command with:
`-with-singularity ABSOLUTE_PATH/scilus-2.1.0.sif`

Image is available [here](http://scil.dinf.usherbrooke.ca/en/containers_list/scilus-2.1.0.sif)

If you are on MacOS or Windows, we recommend using the Docker container to run tractometry_flow pipeline.
Launch your Nextflow command with:
`-with-docker scilus/scilus:2.1.0`

:warning: WARNING :warning:
---------
The official release 2.1.0 is **NOT** available now.

Please, either build the singularity container using this command: 

`singularity build scilus_latest.sif docker://scilus/scilus:latest` 

and then launch your Nextflow command with:
`-with-singularity ABSOLUTE_PATH/scilus_latest.sif`

Or launch your Nextflow command with docker:
`-with-docker scilus/scilus:latest`
