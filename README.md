# Ethraid

Characterize long-period companions with partial orbits.  
Please refer to Van Zandt \& Petigura (2023, in prep.) for details on ethraid's functionality.

## Environment
### Create new environment with python 3.7
- *\$ conda create --name ethraid_env python=3.7*
- *\$ conda activate ethraid_env*

## Download using pip
- *\$ pip install ethraid*
- If the installation fails, try upgrading pip: *\$ curl https://bootstrap.pypa.io/get-pip.py | python*

## Create a configuration file from the template provided and provide the required parameters and desired data
- *\$ cp template_config.py my_config.py*
- NOTE: ethraid uses AU for all distances and M_Jup for all masses. Access helpful conversion factors using e.g.

    ```
    from ethraid import Ms2Mj, pc_in_au
    ```
    which respectively convert solar masses to Jupiter masses and parsecs to AU, in your config file.

## Example usage
### CLI: Run orbit fits, plot results, and print 95\% confidence intervals all at once from the command line
- *\$ ethraid all -cf path/to/my_config.py* -rfp results/\{star_name\}/\{star_name\}_processed.h5 -t 1d 2d
    - Note that the *-rfp* (read file path) flag requires the path to the output directory where the fit results are stored. On a first run, this path *does not exist yet,* but it will be created after the fit and before plotting.

### Alternative run each command separately:

#### Run orbit fits using parameters in configuration file
- *\$ ethraid run -cf path/to/my_config.py*
#### Load and plot saved results
- *\$ ethraid plot -cf path/to/my_config.py -rfp results/\{star_name\}/\{star_name\}_raw.h5 -gn 75 -t 1d 2d*
#### Print 95\% mass and semi-major axis confidence intervals based on derived posterior
- *\$ ethraid lims -rfp results/\{star_name\}/\{star_name\}_raw.h5 -gn 100*

### Another alternative: use the api_run.py module to interface easily with the API

- Add the following code to the end of the module to run a fit, plot the results, and print a summary to the command line.

    ```
    if __name__ == "__main__":
    
        config_path = 'path/to/my_config.py
        read_file_path = 'results/\{star_name\}/\{star_name\}_processed.h5'
    
    
        plot=True
        verbose = True
    
        run(config_path, read_file_path,
            plot=plot, verbose=verbose)
    ```

## Download repo from Github
### Install dependencies using requirements.txt 
- *\$ pip install -r requirements.txt*

### Build code from top level of repo
- *\$ cd trends/*
- *\$ python setup.py build_ext --inplace*

### Use api_run.py as a reference
- *\$ python api_run.py*
