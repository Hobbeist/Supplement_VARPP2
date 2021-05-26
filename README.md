# Supplement: Interpretable pathogenic variant prediction for rare disease diagnostics 
Authors: Sebastian Rauschert, Denise Anderson, Yiming Wang, Gareth Baynam, Timo Lassmann

## Content of this repo
This repository contains the code and description to create the pathogenic and benign variant files used in the paper.
Further, it contains the process as to how we get the HPO term associated genes, and how to prepare the input data for the run of `VARPP2` and `VARPP-RuleFit`.

### Sampling code
The code to sample down the number of benign variants can be found under `code/sample_benign.R`

### Custom function from the `yarn` package
To take advantage of the GTExV8, we customised the `downloadGTEx()` function for the version 8.  
The code for this can be found in `code/yarn_custom.R`
