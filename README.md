
# Trend Analysis


## Installing dependencies
- Create a new environment (python=3.9.4)
- *\$ conda install -file requirements.txt*
- *\$ pip install radvel celerite*
- *\$ python setup_all.py build_ext --inplace*

## Profiling rv_post

- In root trends directory
- *\$ kernprof -l -v no_classes/runner.py*
