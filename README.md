
# Ethraid

Characterize long-period companions with partial orbits.

[![Powered by radvel](https://img.shields.io/badge/powered_by-radvel-EB5368.svg?style=flat)](https://radvel.readthedocs.io/en/latest/)

## Environment
### Create new environment with python 3.7
- *\$ conda create --name trends_env python=3.7*
- *\$ conda activate trends_env*

## Download using pip
- *\$ pip install ethraid*
- You may have to first upgrade pip using *\$ curl https://bootstrap.pypa.io/get-pip.py | python*

## Example CLI usage
### Run orbit fits using parameters in configuration file
- *\$ ethraid run -cf ethraid/config_files/test1.py*
### Load and plot saved results
- *\$ ethraid plot -cf ethraid/config_files/test1.py -rfp results/test1/test1_raw.h5 -gn 100*
### Print 95\% confidence intervals of companion's semi-major axis and mass based on derived posterior
- *\$ ethraid less -rfp results/test1/test1_raw.h5 -gn 100*

## If downloading repo
### Install dependencies using requirements.txt 
- *\$ pip install -r requirements.txt*

### Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

## Use api_run.py as a reference
- *\$ python api_run.py*
