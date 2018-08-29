# Cell Ranger Docker Image

This docker image installs Cell Ranger and all of its dependencies needed to run `cellranger count`.
We are using Cell Ranger from https://github.com/10XGenomics/cellranger. This is currently Cell Ranger 2.1.1:w

This is licensed with an MIT license and so is most compatible with HCA data.

Note the following system requirements for running Cell Ranger: https://support.10xgenomics.com/single-cell-gene-expression/software/overview/system-requirements


## Build docker image

Build container:

`docker build .`

## Interactively work with image

`docker run -it <image> /bin/bash`

## Verify with the following command:

`docker run <image>`

