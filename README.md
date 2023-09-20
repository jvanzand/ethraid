
# Ethraid

Characterize long-period companions with partial orbits.

## Environment
### Create new environment with python 3.7
- *\$ conda create --name trends_env python=3.7*
- *\$ conda activate trends_env*

## Download using pip
- *\$ pip install ethraid*
- If the installation fails, try upgrading pip: *\$ curl https://bootstrap.pypa.io/get-pip.py | python*

## Example CLI usage
### Run orbit fits, plot results, and print 95\% confidence intervals all at once
- *\$ ethraid all -cf path/to/ethraid/example_config_files/config_191939.py* -rfp results/test1/test1_processed.h5 -t 1d 2d

### Or run each command separately:

### Run orbit fits using parameters in configuration file
- *\$ ethraid run -cf path/to/ethraid/example_config_files/config_191939.py*
### Load and plot saved results
- *\$ ethraid plot -cf ethraid/config_files/test1.py -rfp results/test1/test1_raw.h5 -gn 75 -t 1d 2d*
### Print 95\% mass and semi-major axis confidence intervals based on derived posterior
- *\$ ethraid lims -rfp results/test1/test1_raw.h5 -gn 100*

## Download repo from Github
### Install dependencies using requirements.txt 
- *\$ pip install -r requirements.txt*

### Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

### Use api_run.py as a reference
- *\$ python api_run.py*
