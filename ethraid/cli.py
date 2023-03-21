"""
Command Line Interface
"""

from argparse import ArgumentParser
import ethraid.driver

# ###############
# import ethraid.system_params as sp
# ###############

def main():
    psr = ArgumentParser(
              description='Ethraid: Find distant companions with partial orbits', 
              prog='Ethraid'
                         )
    psr.add_argument('--version',
                    action='version',
                    version="%(prog)s {}".format(ethraid.__version__),
                    help='Print version number and exit.'
                    )
                    
    subpsr = psr.add_subparsers(title='subcommands', dest='subcommand')
    
    
    ## Parent parser to define arguments common to multiple subcommands
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
                            help='Print out warnings and (possibly) useful information')

    
    ## Run trend simulation
    #######################
    psr_run = subpsr.add_parser('run', parents=[psr_parent],
                                description="Derive a-m posterior by forward modelling orbits",
                                prefix_chars="-"
                                )

    psr_run.add_argument('-sn', '--star_name', 
                        type=str, 
                        required=True,
                        help='Name of host star for file naming. Need not be official.')
    psr_run.add_argument('-ms', '--m_star', 
                        type=float,  
                        required=True,
                        help='Stellar mass in units of Jupiter masses')
    psr_run.add_argument('-ds', '--d_star', 
                        type=float, 
                        required=True,
                        help='Distance from Earth to the host star in AU')
            
    psr_run.add_argument('-gd', '--gammadot', 
                        type=float,  
                        required=True,
                        help='Linear trend in RVs in m/s/day')
    psr_run.add_argument('-gde', '--gammadot_err', 
                        type=float, 
                        required=True,
                        help='Error on gamma_dot')
                        
    # required=False so we can have targets with no curvature
    # default=0/1e8 so that if the flag isn't provided, we have values to pass to RV code
    # nargs='?' so I can pass a blank argument into Cadence (useful for running lists of targets)
    # const=0/1e8 so that if the flag IS passed but with NO args (like on Cadence), default values are filled in
    psr_run.add_argument('-gdd', '--gammaddot', 
                        type=float, 
                        required=False, 
                        default=0, 
                        nargs='?', 
                        const=0,
                        help='Curvature in RVs in m/s/day/day')
    psr_run.add_argument('-gdde', '--gammaddot_err', 
                        type=float, 
                        required=False, 
                        default=1e8, 
                        nargs='?', 
                        const=1e8,
                        help='Error on gamma_ddot')

    psr_run.add_argument('-rvb', '--rv_baseline', 
                        type=float, 
                        required=True,
                        help='Length of RV time baseline in days')
    psr_run.add_argument('-rvep', '--rv_epoch', 
                        type=float,  
                        required=True,
                        help='Epoch of RV timeseries in BJD, usually around baseline midpoint')
                
    # required=False for targets with no astrometry
    # nargs='?' so I can pass blank args for Cadence (because in a list of stars, some don't have dmu and some do; my .sh script has to provide a single command with which all stars are run)
    # default=None and const=None by default, no need to set
    psr_run.add_argument('-dmu', '--delta_mu', 
                        type=float,  
                        required=False, 
                        nargs='?',
                        help='Change in astrometric proper motion in milliarcseconds/yr')
    psr_run.add_argument('-dmue', '--delta_mu_err', 
                        type=float,  
                        required=False, 
                        nargs='?',
                        help='Error on dmu')

    # required=False for targets without imaging
    # nargs='?' to let me pass the flag with no args for Cadence
    # default=None and const=None by default, no need to set
    psr_run.add_argument('-vmag', '--vmag', 
                        type=float, 
                        required=False, 
                        nargs='?',
                        help='Apparent V-band magnitude of host star')
    psr_run.add_argument('-imwav', '--imag_wavelength', 
                        type=float,  
                        required=False, 
                        nargs='?',
                        help='Wavelength of imaging observations (μm)')
    psr_run.add_argument('-cs', '--contrast_str', 
                        type=str,  
                        required=False, 
                        nargs='?',
                        help='Path to dataframe containing contrast curve')            
    psr_run.add_argument('-n', '--num_points', 
                        type=int, 
                        required=False, 
                        default=1000000,
                        help='Number of orbit models to run')
    psr_run.add_argument('-gn', '--grid_num', 
                        type=int,  
                        required=False, 
                        default=100,
                        help='Dimension of binned probability array')

    # # The 'save' arg is special. Rather than having an assigned type, the "action='store_true'" argument assumes 1) boolean type and 2) default = False if the flag is not given. It also allows there to be no argument after the flag (in which case the value is set to True).
    psr_run.add_argument('-s', '--save',
                        nargs='+',
                        default='proc',
                        choices=['raw', 'proc'], 
                        required=False,
                        help = 'How to save probability arrays. \n'
                               'Pass no flag to save processed arrays. \n'
                               'Pass with "raw" or "proc" flags to save raw, processed, or both.'
                        )



    psr_run.set_defaults(func=ethraid.driver.run)
    
    
    ## Plot results of trend calculations
    #####################################
    psr_plot = subpsr.add_parser('plot', parents=[psr_parent],
                                 description="Plot joint 2d and marginalized 1d posteriors",
                                 prefix_chars="-"
                                )
    psr_plot.add_argument('-rfp', '--read_file_path', 
                          type=str,  
                          required=True,
                          help='File path to read in already-calculated posterior array'
                          )
    psr_plot.add_argument('-t', '--type',
                          type=str, nargs='+',
                          required=True,
                          choices=['2d', '1d'],
                          help='Generate either 2D or 1D posteriors'
                          )
    psr_plot.add_argument('-gn', '--grid_num', 
                          type=int,  
                          required=False, 
                          default=None,
                          help='Dimension of binned probability array. Required for raw  input arrays.'
                          )
    # nargs='+' so I can use the same flag to pass both semi_major axis and mass of a companion to plot
    psr_plot.add_argument('-sp', '--scatter_plot', 
                        type=float,  
                        required=False, 
                        nargs='+',
                        help='Semi_major axis (AU) and mass (M_J) of a known companion for 2d plot.\
                        Separate values by a space only.')
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
    psr_less.add_argument('-gn', '--grid_num',
                          type=int,
                          required=False,
                          default=100,
                          help='Dimension of binned probability array. Required for raw  input arrays.'
                          )
    psr_less.set_defaults(func=ethraid.driver.less)

    
    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
