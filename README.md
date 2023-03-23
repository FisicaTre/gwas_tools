# gwas_tools [![PyPI version](https://badge.fury.io/py/gwas-tools.svg)](https://badge.fury.io/py/gwas-tools)
Utils functions to build a pipeline for scattered light noise hunting in gravitational waves detectors.

## Docs [![Documentation Status](https://readthedocs.org/projects/gwas-tools/badge/?version=latest)](https://gwas-tools.readthedocs.io/en/latest/?badge=latest)

[Code Documentation](https://gwas-tools.readthedocs.io/en/latest/)

##  Requirements

- Python 3.7+
- [requirements.txt](requirements.txt)
- Optional requirements
  - `gwdama`

## Description
- `algorithms`
  - main algorithms of the currently used pipelines.
- `automation`
  - functions to create the main scripts to be executed on distributed computing systems.
- `helpers`
  - functions to get inputs for the pipeline, and classes to build a `.sub` and a `.dag` file to submit the pipeline.
- `html`
  - class to build a html page to summarize the output of the pipeline.
- `summary_pages`
  - functions to create the summary pages of the currently used pipelines.
- `utils`
  - utils functions to process signals, manage output files, and plot pipeline's results.

## References
- [Application to LIGO data](https://iopscience.iop.org/article/10.1088/1361-6382/aa8e6b/meta) 
- [Application to Virgo data](https://iopscience.iop.org/article/10.1088/1361-6382/ab9719/meta) 
- [Omega algorithm](https://dspace.mit.edu/handle/1721.1/34388)
- [Time varying filter EMD](https://www.sciencedirect.com/science/article/pii/S0165168417301135?casa_token=e9Q5Bi85etAAAAAA:ow686chMeVLYYF4anHGXpMx_dNSzej0s3x9PJuCuyt1zYyyyYLUsOOw6VSWXQJgZQPgAUitW3IU)
- [Adaptive time series analysis](https://www.jstor.org/stable/pdf/53161.pdf?casa_token=ZqoSg2aXRR8AAAAA:c-vPcJu5-ymb9Z_zZmr3pD1twXy3pb7nBxyUN0oUUoJfKgVLX1MIQhGqovwLsNJFQSCDrXa3k7GFJPxfIJhkwAXO650sblUb3mnVphXSjg73yUpczlEj)

