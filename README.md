Tractometry pipeline
====================

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/scilus/tractometry_flow)](https://github.com/scilus/tractometry_flow/releases)

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
- [scilpy](https://github.com/scilus/scilpy)

Singularity/Docker
-----------
If you are on Linux, we recommend using the Singularity to run tractometry_flow pipeline.
If you have Singularity == 3.*, launch your Nextflow command with:
`-with-singularity scilus/scilus:1.2.0_tractometryflow-1.0.0`

If you have rebuild singularity Singularity == 2.* image is available [here](http://scil.dinf.usherbrooke.ca/en/containers_list/scilus-1.2.0_tractometryflow-1.0.0.img)
Launch your Nextflow command with: `-with-singularity ABSOLUTE_PATH/scilus-1.2.0_tractometryflow-1.0.0.img`

If you are on MacOS or Windows, we recommend using the Docker container to run tractometry_flow pipeline.
Launch your Nextflow command with:
`-with-docker scilus/scilus:1.2.0_tractometryflow-1.0.0`

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`
