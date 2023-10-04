# eQTLGen data preprocessing for TOPMed RNA-seq paper

## Dependencies

* Singularity (v. >=3)
* NextFlow (v. >= 21.04.0)


## Setup

1. Pull a singularity container with some dependencies: `make singularity`
2. Fetch data: `make data`
3. Update nextflow.config as necessary for your compute environment.


## Running

1. Get the top hit per gene: `make top-hit-per-gene`
2. Lift and tabix: `make lift-and-tabix`