# Steps for making Docker image:

We are using cell ranger from https://github.com/10XGenomics/cellranger . 
This is currently Cellranger 2.1.1:w

This is licensed with an MIT license and so is most compatible with HCA data.

## Build docker image

Build container:

`docker build .`

## Interactively work with image

`docker run -it <image> /bin/bash`

## Verify with the following command:

`docker run <image>`