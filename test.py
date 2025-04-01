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
    calc_error_count = None
    load_error_count = None
    
    if calc:
        # Tests of full calculations
        # Two separate config files testing different options
        config_files = ['test_config_files/test1.py',
                        'test_config_files/test2.py',
                        'test_config_files/test3.py'
                        ]
        calc_error_count = 0
        for cf in config_files:
            try:
                # import pdb; pdb.set_trace()
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
            
            test_num = rfp.split('_')[0][-1] # Get number 1, 2, or 3
            config_path = 'test_config_files/test{}.py'.format(test_num) # Get config path. A little ad hoc
    
            try:
                # Note that I provide arguments after rfp, but they are only needed sometimes.
                # If the loaded file has scatter_plot or outdir, eg, then the loaded values supersede those given here.
                api_run.run(config_path, read_file_path=rfp, plot=plot, verbose=verbose)
                
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

def cli_tester(calc=True, load=True, all_=True):
    """
    This function tests the command line interface.
    First it accesses two different configuration files
    and calculates their posteriors using the CLI framework.
    Then it accesses the resulting saved files and runs
    the plot and lims functions.
    Finally, it repeats the above by running the "all"
    command on each config file/output file pair.
    
    Note that if this function is first run with
    calc=False and load=True, the necessary files
    will not yet exist to load.
    """
    calc_error_count = None
    load_error_count = None
    all_error_count = None
    
    if calc:
        # Tests of full calculations
        # Two separate config files testing different options
        config_files = ['test_config_files/test1.py',
                        'test_config_files/test2.py',
                        'test_config_files/test3.py'
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
                config_name = 'test_config_files/{}.py'.format(test_name)
                
                # Run plot function
                func = "plot"
                subprocess.run(["python", "ethraid/cli.py", func, "-cf", config_name,
                                "-t", "1d", "2d", "-rfp", rfp], check=True)
                # Run lims function
                func = "lims"
                subprocess.run(["python", "ethraid/cli.py", func,
                                "-cf", config_name, "-rfp", rfp])
                                
            except Exception as err:
                print(err)
                load_error_count+=1
                print('')
                rfp_short = rfp.split('/')[-1] # Eg 'test_2_raw.h5'
                print('test.cli_tester: Error executing "{}" on arrays in file {}.'\
                      .format(func, rfp_short))
    
        print("{} errors encountered while loading/running saved arrays from CLI.".format(load_error_count))
        
       
    if all_: # Run the run, plot, and lims functions sequentially
        print("Running all")
        all_error_count = 0
        for test_num in ['1','2','3']:
            for file_type in ['raw', 'processed']:
            
                try:
                    config_name = 'test_config_files/test{}.py'.format(test_num) # Eg test1
                    rfp = 'results/test{0}/test{0}_{1}.h5'.format(test_num, file_type) # Eg test3_processed.h5
                    # Run all function
                    subprocess.run(["python", "ethraid/cli.py", "all", "-cf", config_name,
                                    "-t", "1d", "2d", "-rfp", rfp], check=True)
            
                except Exception as err:
                   all_error_count+=1
                   print('')
                   rfp_short = rfp.split('/')[-1] # Eg 'test_2_raw.h5'
                   print('test.cli_tester: Error executing "all" on arrays in file {}.'\
                         .format(rfp_short))
        
    
    
    return calc_error_count, load_error_count, all_error_count

if __name__=="__main__":
    
    api_errs = api_tester(calc=True, load=True)
    cli_errs = cli_tester(calc=True, load=True, all_=True)
    
    print("\n")
    print("Test complete:")
    type_list = ["API array calculation", "API array loading", 
                 "CLI array calculation", "CLI array loading", "CLI 'all' function"]
    for i, err_count in enumerate([*api_errs, *cli_errs]):
        use_type = type_list[i]
        if err_count==None:
            print("    Did not perform "+use_type)
        elif err_count==0:
            print("    0 errors encountered running "+use_type)
        elif err_count>0:
            print("    {} errors encountered running ".format(err_count)+use_type)
    
    
    
    
    
    
    
    