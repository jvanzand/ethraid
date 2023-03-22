
# Trend Analysis


## Environment
- Create new environment with python 3.7
- *\$ conda create --name trends_env python=3.7*
- *\$ conda activate trends_env*

## Download from Test PyPI (regular pip download coming soon)
- *\$ python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ethraid*
- You may have to upgrade pip using *\$ curl https://bootstrap.pypa.io/get-pip.py | python*

## If downloading repo
- Install dependencies using requirements.txt 
- *\$ pip install -r requirements.txt*

- Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

## To profile code or run without profiling
- *\$ kernprof -l -v runner.py*
- or
- *\$ python runner.py*
