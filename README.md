
# Trend Analysis


## Environment
- Create new environment with python 3.7
- *\$ conda create --name trends_env python=3.7*
- *\$ conda activate trends_env*

## Download using pip
- *\$ pip install ethraid*
- You may have to first upgrade pip using *\$ curl https://bootstrap.pypa.io/get-pip.py | python*

## Example CLI usage
- *\$ ethraid run -sn 191939_cli -ms 845.4 -ds 11088795.98 -gd 0.114 -gde 0.006 -gdd -0.00006 -gdde 0.000019 -rvb 430.25 -rvep 2458847.78 -dmu 0.1277 -dmue 0.0342 -vmag 8.97 -imwav 2.2 -cs 'path/to/contrast_curve.csv' -n 1000000 -gn 100 -v -s proc raw*

## If downloading repo
- Install dependencies using requirements.txt 
- *\$ pip install -r requirements.txt*

- Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

## Use api_run.py as a reference
- *\$ python api_run.py*
