"""
Command Line Interface
"""

from argparse import ArgumentParser
import ethraid.driver


def main():
    psr = ArgumentParser(
              description='ethraid: Find distant companions with partial orbits', 
              prog='ethraid'
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

    
    ## run: Run trend simulation
    #######################
    psr_run = subpsr.add_parser('run', parents=[psr_parent],
                                description="Derive a-m posterior by forward modeling orbits",
                                prefix_chars="-"
                                )
    psr_run.add_argument('-cf', '--config',
                         type=str,
                         required=True,
                         help='Relative path of configuration file.'
                         )

    psr_run.set_defaults(func=ethraid.driver.run)
    
    
    ## plot: Plot results of trend calculations
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
    # grid_num not required for plot, but loading *raw* arrays without it will produce an error.
    # psr_plot.add_argument('-gn', '--grid_num',
    #                       type=int,
    #                       required=False,
    #                       default=100,
    #                       help='Dimension of binned probability arrays. Required for raw input arrays.'
    #                       )
    psr_plot.set_defaults(func=ethraid.driver.plot)
    
    
    ## lims: Print only the 1D limits calculated from the total posterior
    ###############################################################
    psr_lims = subpsr.add_parser('lims', parents=[psr_parent],
                                 description="Print 95% a-m confidence intervals from calculated posteriors",
                                 prefix_chars="-"
                                )
    psr_lims.add_argument('-rfp', '--read_file_path', 
                        type=str,  
                        required=True,
                        help='File path to read in already-calculated posterior array')
    psr_lims.add_argument('-cf', '--config',
                         type=str,
                         required=True,
                         help='Relative path of configuration file.'
                         )
    # grid_num not required for plot, but loading *raw* arrays without it will produce an error.
    # psr_lims.add_argument('-gn', '--grid_num',
    #                       type=int,
    #                       required=False,
    #                       default=100,
    #                       help='Dimension of binned probability array. Required for raw input arrays.'
    #                       )
    psr_lims.set_defaults(func=ethraid.driver.lims)
    
    ## all: Run the run, plot, and lims commands sequentially
    ###############################################################
    psr_all = subpsr.add_parser('all', parents=[psr_parent],
                                 description="Run the run, plot, and lims commands sequentially",
                                 prefix_chars="-"
                                )
    psr_all.add_argument('-cf', '--config',
                         type=str,
                         required=True,
                         help='Relative path of configuration file.'
                         )
    psr_all.add_argument('-rfp', '--read_file_path',
                          type=str,
                          required=True,
                          help='File path to read in already-calculated posterior array'
                          )
    psr_all.add_argument('-t', '--type',
                          type=str, nargs='+',
                          required=True,
                          choices=['1d', '2d'],
                          help='Generate either 2D or 1D posteriors'
                          )
    # # grid_num not required for plot, but loading *raw* arrays without it will produce an error.
    # psr_all.add_argument('-gn', '--grid_num',
    #                       type=int,
    #                       required=False,
    #                       default=100,
    #                       help='Dimension of binned probability arrays. Required for raw input arrays.'
    #                       )
    psr_all.set_defaults(func=ethraid.driver.all)

    
    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
