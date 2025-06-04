# HMMintrogression-smk
### This is a snakemake pipeline for infering archaic introgression segments using the HMM model described in Seguin-Orlando et al. (2014) while also adapting the training functionality from Skov et al. (2018).

## Required software
#### Python (v3.10)
#### Snakemake (v7.2)
#### BCFtools (v1.21)
#### Numpy (v1.19.5) 
#### PyVCF (v0.6.8)
#### Pandas (v1.1.5)

## Input data
#### All required input data are listed in the config.yaml file.

## Usage
```
snakemake -j 10 -w 60 -kp
```
