## gwadaptive_scattering  
Utils functions to build a pipeline for scattered light noise hunting in gravitational waves detectors.

#### Version

v0.6.2

####  Requirements

- Python 3.6+
- [requirements.txt](requirements.txt)

#### Description
- `helpers`
  - functions to get inputs for the pipeline, and classes to build a `.sub` and a `.dag` file to submit the pipeline.
- `html`
  - class to build a html page to summarize the output of the pipeline.
- `utils`
  - utils functions to process signals, manage output files, and plot pipeline's results.

#### References
- Application to LIGO data: [Link](https://iopscience.iop.org/article/10.1088/1361-6382/aa8e6b/meta) 
- Application to Virgo data: [Link](https://iopscience.iop.org/article/10.1088/1361-6382/ab9719/meta) 
- Omega algorithm: [Link](https://dspace.mit.edu/handle/1721.1/34388)
- Time varying filter EMD methodological paper: [Link](https://www.sciencedirect.com/science/article/pii/S0165168417301135?casa_token=e9Q5Bi85etAAAAAA:ow686chMeVLYYF4anHGXpMx_dNSzej0s3x9PJuCuyt1zYyyyYLUsOOw6VSWXQJgZQPgAUitW3IU)
- Introduction to adaptive time series analysis: [Link](https://www.jstor.org/stable/pdf/53161.pdf?casa_token=ZqoSg2aXRR8AAAAA:c-vPcJu5-ymb9Z_zZmr3pD1twXy3pb7nBxyUN0oUUoJfKgVLX1MIQhGqovwLsNJFQSCDrXa3k7GFJPxfIJhkwAXO650sblUb3mnVphXSjg73yUpczlEj)

