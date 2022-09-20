
# Trend Analysis


## Environment and dependencies
- Create new environment with python 3.9
- *\$ conda create --name trends_env python=3.9*
- *\$ conda activate trends_env*
- Install dependencies using requirements.txt
- *\$ pip install -r requirements.txt*

## Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

## To profile code or run without profiling
- *\$ kernprof -l -v runner.py*
- or
- *\$ python runner.py*
