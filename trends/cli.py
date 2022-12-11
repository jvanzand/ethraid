"""
Command Line Interface
"""

from argparse import ArgumentParser
from runner import run

###############
import trends.system_params as sp
###############

parser = ArgumentParser(description="Perform trend analysis")

parser.add_argument('-sn', '--star_name', type=str, metavar='\b', required=True,
                    help='Name of host star')
parser.add_argument('-ms', '--m_star', type=float, metavar='\b', required=True,
                    help='Stellar mass in units of Jupiter masses')
parser.add_argument('-ds', '--d_star', type=float, metavar='\b', required=True,
                    help='Distance to the host star in AU')
                    
parser.add_argument('-gd', '--gdot', type=float, metavar='\b', required=True,
                    help='Linear trend in RVs in m/s/day')
parser.add_argument('-gde', '--gdot_err', type=float, metavar='\b', required=True,
                    help='Error on gamma_dot')

parser.add_argument('-bl', '--baseline', type=float, metavar='\b', required=True,
                    help='Length of RV time baseline in days')
parser.add_argument('-rvep', '--rv_epoch', type=float, metavar='\b', required=True,
                    help='Epoch of RV timeseries in BJD, usually around baseline midpoint')
                    
# required=False so we can have targets with no curvature
# default=0/1e8 so that if the flag isn't passed, we have values to pass to RV code
# nargs='?' so I can pass a blank argument into Cadence (useful for running lists of targets)
# const=0/1e8 so that if the flag IS passed but with NO args (like on Cadence), default values are filled in
parser.add_argument('-gdd', '--gddot', type=float, metavar='\b', 
                    required=False, default=0, nargs='?', const=0,
                    help='Curvature in RVs in m/s/day/day')
parser.add_argument('-gdde', '--gddot_err', type=float, metavar='\b', 
                    required=False, default=1e8, nargs='?', const=1e8,
                    help='Error on gamma_ddot')
                    
# required=False for targets with no astrometry
# nargs='?' so I can pass blank args for Cadence (because in a list of stars, some don't have dmu and some do; my .sh script has to provide a single command with which all stars are run)
# default=None and const=None by default, no need to set
parser.add_argument('-dmu', '--delta_mu', type=float, metavar='\b', 
                    required=False, nargs='?',
                    help='Change in astrometric proper motion in milli-arcseconds/yr')
parser.add_argument('-dmue', '--delta_mu_err', type=float, metavar='', 
                    required=False, nargs='?',
                    help='Error on dmu')

# required=False for targets without imaging
# nargs='?' to let me pass the flag with no args for Cadence
# default=None and const=None by default, no need to set
parser.add_argument('-vmag', '--vmag', type=float, metavar='',
                    required=False, nargs='?',
                    help='Apparent V-band magnitude of host star')
parser.add_argument('-imwav', '--imag_wavelength', type=float, metavar='', 
                    required=False, nargs='?',
                    help='Wavelength of imaging observations (μm)')
parser.add_argument('-cs', '--contrast_str', type=str, metavar='', 
                    required=False, nargs='?',
                    help='Path to dataframe containing contrast curve')

# nargs='+' so I can pass both semi-major axis and mass with the same flag
parser.add_argument('-sp', '--scatter_plot', type=float, metavar='\b', required=False, nargs='+',
                    help='Semi-major axis (AU) and mass (M_J) of a known companion to plot. Separate values by a space only.')               
parser.add_argument('-n', '--num_points', type=int, metavar='\b', required=False,
                    help='Number of orbit models to run')
parser.add_argument('-gn', '--grid_num', type=int, metavar='\b', required=False,
                    help='Dimension of binned probability array')

# # The 'save' and 'plot' args are special. Rather than having an assigned type, the "action='store_true'" argument assumes 1) boolean type and 2) default = False if the flag is not given. It also allows there to be no argument after the flag (in which case the value is set to True).
parser.add_argument('-s', '--save', action='store_true', required=False,
                    help='Whether to save posterior files')
parser.add_argument('-p', '--plot', action='store_true', required=False,
                    help='Whether to plot joint (a,m) posterior')
parser.add_argument('-r', '--read', type=str, metavar='\b', required=False,
                    help='File path to read in already-calculated posterior array')

args = parser.parse_args()

if __name__=="__main__":
    
    run(args.star_name, args.m_star, args.d_star, args.gdot, args.gdot_err, 
        args.gddot, args.gddot_err, args.baseline, 
        args.rv_epoch, args.delta_mu, args.delta_mu_err,
        vmag=args.vmag, imag_wavelength=args.imag_wavelength, contrast_str=args.contrast_str,
        scatter_plot=args.scatter_plot,
        num_points=args.num_points, grid_num=args.grid_num, 
        save=args.save, plot=args.plot, read_file_path=args.read)
