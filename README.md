# Supplement: Pathogenic variant prediction in rare disease diagnostics using tissue specific and single cell gene expression data: VARPP v2.0 
Authors: Sebastian Rauschert, Denise Anderson, Timo Lassmann

## Content of this repo
This repository contains the code and description to create the pathogenic and benign variant files used in the paper.
Further, it contains the process as to how we get the HPO term associated genes, and how to prepare the input data for the run of `VARPP` and `RuleFit`.

### Sampling code
The code to sample down the number of benign variants can be found under `code/sample_benign.R`

### Custom function from the `yarn` package
To take advantage of the GTExV8, we customised the `downloadGTEx()` function for the version 8.  
The code for this can be found in `code/yarn_custom.R`
