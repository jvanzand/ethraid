"""
Command Line Interface
"""

from argparse import ArgumentParser
import ethraid.driver


def main():
    psr = ArgumentParser(
              description='Ethraid: Find distant companions with partial orbits', 
              prog='Ethraid'
                         )
    psr.add_argument('-V', '--version',
                    action='version',
                    version="%(prog)s {}".format(ethraid.__version__),
                    help='Print version number and exit.'
                    )
                    
    subpsr = psr.add_subparsers(title='subcommands', dest='subcommand')
    
    
    ## Parent parser to define arguments common to all subcommands
    ###################################################################
    psr_parent = ArgumentParser(add_help=False)
    psr_parent.add_argument('-od', '--outdir', 
                            type=str,  
                            required=False, 
                            default='',
                            help='Path to directory where output files will be saved'
                            )
    psr_parent.add_argument('-v', '--verbose', 
                            action='store_true', 
                            required=False,
                            help='Print out warnings and (possibly) useful information'
                            )

    
    ## Run trend simulation
    #######################
    psr_run = subpsr.add_parser('run', parents=[psr_parent],
                                description="Derive a-m posterior by forward modelling orbits",
                                prefix_chars="-"
                                )
    psr_run.add_argument('-cf', '--config',
                         type=str,
                         required=True,
                         help='Relative path of configuration file.'
                         )

    psr_run.set_defaults(func=ethraid.driver.run)
    
    
    ## Plot results of trend calculations
    #####################################
    psr_plot = subpsr.add_parser('plot', parents=[psr_parent],
                                 description="Plot joint 2d and marginalized 1d posteriors",
                                 prefix_chars="-"
                                )
    psr_plot.add_argument('-cf', '--config',
                          type=str,
                          required=True,
                          help='Relative path of configuration file.'
                          )
    psr_plot.add_argument('-rfp', '--read_file_path',
                          type=str,
                          required=True,
                          help='File path to read in already-calculated posterior array'
                          )
    psr_plot.add_argument('-t', '--type',
                          type=str, nargs='+',
                          required=True,
                          choices=['1d', '2d'],
                          help='Generate either 2D or 1D posteriors'
                          )
    # grid_num not required for plot, but if you load *raw* arrays without it, you will get an error.
    psr_plot.add_argument('-gn', '--grid_num',
                          type=int,
                          required=False,
                          default=None,
                          help='Dimension of binned probability arrays. Required for raw  input arrays.'
                          )
    psr_plot.set_defaults(func=ethraid.driver.plot)
    
    
    ## Print only the 1D bounds calculated from the total posterior
    ###############################################################
    psr_less = subpsr.add_parser('less', parents=[psr_parent],
                                 description="Print 95% a-m confidence intervals from calculated posteriors",
                                 prefix_chars="-"
                                )
    psr_less.add_argument('-rfp', '--read_file_path', 
                        type=str,  
                        required=True,
                        help='File path to read in already-calculated posterior array')
    # grid_num not required for less, but if you load *raw* arrays without it, you will get an error.
    psr_less.add_argument('-gn', '--grid_num',
                          type=int,
                          required=False,
                          default=None,
                          help='Dimension of binned probability array. Required for raw  input arrays.'
                          )
    psr_less.set_defaults(func=ethraid.driver.less)

    
    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
