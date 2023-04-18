# Test module to run through some basic functionality to ensure code works as expected.

import subprocess
import api_run


def api_tester(calc=True, load=True, verbose=False):
    """
    This function tests the api_run.py module.
    First it accesses two different configuration files
    and calculates their posteriors.
    Then it accesses the resulting saved files, plots
    the posteriors, and prints the 2σ bounds.
    
    Note that if this function is first run with
    calc=False and load=True, the necessary files
    will not yet exist to load.
    """
    if calc:
        # Tests of full calculations
        # Two separate config files testing different options
        config_files = ['ethraid/config_files/test1.py',
                        'ethraid/config_files/test2.py',
                        'ethraid/config_files/test3.py'
                        ]
        calc_error_count = 0
        for cf in config_files:
            try:
                api_run.run(cf)
            except Exception as err:
                calc_error_count+=1
                cf_short = cf.split('/')[-1]
                print('test.api_tester: Failed to run config file {}'.format(cf_short))
                print(err)
                print('')
        
        print("{} errors encountered while performing full calculations from API.".format(calc_error_count))

    
    if load:
        # Tests of loaded data
        grid_num = 100
        plot=True
        scatter_plot=[10,10]
        outdir=''
        verbose=False
        
        read_file_paths = ['results/test1/test1_raw.h5',
                           'results/test1/test1_processed.h5',
                           'results/test2/test2_raw.h5',
                           'results/test2/test2_processed.h5',
                           'results/test3/test3_raw.h5',
                           'results/test3/test3_processed.h5'
                           ]

        load_error_count = 0
        for rfp in read_file_paths:
    
            try:
                # Note that I provide arguments after rfp, but they are only needed sometimes.
                # If the loaded file has scatter_plot or outdir, eg, then the loaded values supersede those given here.
                api_run.run(read_file_path=rfp,
                            grid_num=grid_num, plot=plot, 
                            scatter_plot=scatter_plot, 
                            outdir=outdir, verbose=verbose)
                
                if verbose:
                    print("test.api_tester: Successfully loaded file {}".format(rfp.split('/')[-1]))
                    print('')
            except Exception as err:
                load_error_count+=1
                rfp_short = rfp.split('/')[-1]
                print('test.api_tester: Error loading file {}'.format(rfp_short))
                print(err)
                
                import traceback
                traceback.print_exc()
                print('')
                

        print("{} errors encountered while loading/running saved arrays from API.".format(load_error_count))
    
    
    return calc_error_count, load_error_count

def cli_tester(calc=True, load=True, verbose=False):
    """
    This function tests the command line interface.
    First it accesses two different configuration files
    and calculates their posteriors using the CLI framework.
    Then it accesses the resulting saved files and runs
    the plot and less functions.
    
    Note that if this function is first run with
    calc=False and load=True, the necessary files
    will not yet exist to load.
    """
    if calc:
        # Tests of full calculations
        # Two separate config files testing different options
        config_files = ['ethraid/config_files/test1.py',
                        'ethraid/config_files/test2.py',
                        'ethraid/config_files/test3.py'
                        ]
                        
        calc_error_count = 0
        for cf in config_files:
            try:
                subprocess.run(["python", "ethraid/cli.py", "run", "-cf", cf])
            except Exception as err:
                calc_error_count+=1
                cf_short = cf.split('/')[-1]
                if verbose:
                    print('test.cli_tester: Error executing "run" on config file {}'.format(cf_short))
                    print('')
        
        print("{} errors encountered while performing full calculations from CLI.".format(calc_error_count))
    
    if load:
        # Tests of loaded data                             
        read_file_paths = ['results/test1/test1_raw.h5',
                           'results/test1/test1_processed.h5',
                           'results/test2/test2_raw.h5',
                           'results/test2/test2_processed.h5',
                           'results/test3/test3_raw.h5',
                           'results/test3/test3_processed.h5'
                           ]
                          
        load_error_count = 0    
        for rfp in read_file_paths:
            try:
                test_name = rfp.split('/')[-1].split('_')[0] # Eg 'test1'
                config_name = 'ethraid/config_files/{}.py'.format(test_name)
                
                func="plot"
                subprocess.run(["python", "ethraid/cli.py", func, "-cf", config_name,
                                "-t", "1d", "2d", "-rfp", rfp, "-gn", "100"], check=True)
                func="less"
                subprocess.run(["python", "ethraid/cli.py", func,
                                "-rfp", rfp, "-gn", "100"])
                                
            except Exception as err:
                load_error_count+=1
                print('')
                rfp_short = rfp.split('/')[-1] # Eg 'test_2_raw.h5'
                print('test.cli_tester: Error executing "{}" on arrays in file {}.'\
                      .format(func, rfp_short))
    
        print("{} errors encountered while loading/running saved arrays from CLI.".format(load_error_count))
    
    
    return calc_error_count, load_error_count

if __name__=="__main__":
    
    api_errs = api_tester()
    cli_errs = cli_tester()
    
    print("{} api, {} cli errors".format(api_errs, cli_errs))
    
    
    
    
    
    
    
    