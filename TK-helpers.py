#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
TK-helpers - a simple python3 script wrapping TKGWV2's helper functions.

    ## Usage:
    ```
    user@desktop:~$ ./TK-helpers.py <submodule> [args]
    ```

    ## Example:
    ```
    user@desktop:~$ ./TK-helpers.py downsampleBam --downsampleN 1800000    \
                                                  --extensionBam "\\.bam$" \
                                                  --suffixDownBam "_sub" 
    ```
"""

__author__  = 'Maël Lefeuvres'
__license__ = 'GPL-2.0'
__version__ = '1.0c'

import argparse
from subprocess      import Popen
from os.path         import isfile, dirname, realpath
from io              import TextIOWrapper
from distutils.spawn import find_executable

def is_available(software_name):
    """Check whether `name` is on PATH."""
    return find_executable(software_name) is not None


def check_plink_basename(
        basename,
        required_extensions=['.bed', '.bim', '.fam']
    ):
    """
    Check the user-provided plink basename, and ensure all the required 
    extensions are present within the filesystem 
    (i.e.: 'basename.bed', 'basename.bim', 'basename.fam').
    """
    if not all([isfile(basename + ext) for ext in required_extensions]):
        raise argparse.ArgumentTypeError(
            "Provided basename '{}' is missing one of their required "
            "file extension: {}".format(basename, required_extensions)
        )
    return basename


def format_args_to_R(name, value):
    """
    Convert from Python to Rscript type formatting. This is mainly intended
    for lists and booleans, and works best when manipulating namespaces or 
    dictionaries.

        ## Usage:
        ```
        args = dict(x=['a','b','c'], y=False)
        [format_args_to_R(key, value) for (key, value) in args.items()]
        ```

        ## Attributes:
          - name <str>          : Name of the command line argument.
          - value <data struct> : User-provided value of the argument.

        ## Notes:
        - integers,floats and strings are expected to pass through unmodified.
     
        ## Examples:
        |         Python         |         Rscript         |
        + ---------------------- + ----------------------- |
        |  list=["a", "b", "c"]  |  list=c("a", "b", "c")  |
        |  bool=True             |  bool=TRUE              |
        |  string="\\.bam$"      |  string="\\.bam$"       |
        |  value=None            |  value=NULL             |

    """
    def python_to_r(value):
        argtype=type(value)
        if argtype == list:
            r_value = r"c({})".format(', '.join(r'{0}'.format(
                python_to_r(item)) for item in value
            ))
        elif argtype == bool:
            r_value = str(value).upper()
        elif argtype == str:
            r_value = '"' + repr("'\0" + value)[6:]
        elif argtype == TextIOWrapper:
            r_value = r'"{}"'.format(value.name)
        elif argtype in [int, float]:
            r_value = value
        elif value is None:
            r_value = "NULL"

        else:
            raise RuntimeError("Failed to parse '{}' into a valid R value".format(value))
        return r_value

    return r"{}={}".format(name, python_to_r(value))


def format_R_command(args):
    """
    Generate a set of R instructions, given the provided user-input.
        ## Usage:
        ```
        args = Namespace(helper='downsampleBam',
                         downsampleN=1500000,
                         extensionBam='\\.bam$',
                         suffixDownBam='_subsampled'
        )
        subcommand = format_R_command(args)
        print(subcommand)
        ```

        ## Attributes:
          - args (Namespace): A set of command line arguments, generated using
                              argparser.
    """

    helper_script = "{pywd}/helpers/{helper}.R".format(
        pywd=dirname(realpath(__file__)),
        helper=args.helper
    )
    subcmd = r'source("{}"); '.format(helper_script)

    SKIPPED_PARAMS = ["helper", "verbose"]
    params = ', '.join([
        format_args_to_R(arg, getattr(args, arg)) 
        for arg in vars(args) 
        if arg not in SKIPPED_PARAMS
    ])

    subcmd += r"{}({})".format(args.helper, params)

    return subcmd

def call_helper(subcmd):
    R_CALLER = "/usr/bin/env Rscript"
    cmd = r"{} --vanilla -e '{}'".format(R_CALLER, subcmd)
    print("\nRunning: {}\n".format(cmd))
    Popen(cmd, shell=True).wait()


def build_parser():
    """
    Create the main command line parser + a subcommand for each provided TKGWV2
    helper script.
    """

    # Main parser
    parser = argparse.ArgumentParser(
        prog="TK-helpers.py",
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=25
        ),

        description="""TK-helpers - A simple utility script wrapping TKGWV2's
                       helper functions.
                    """
    )
    subparsers = parser.add_subparsers(dest='helper', required=True)

    
    # Build subparser for distSimulations.R
    dist_simulations_parser = subparsers.add_parser(
        'distSimulations', 
        help="""
             Generate simulated distribution curves and posterior probabilities
             for estimated HRC values based on an input plink frequencies file.
             """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    dist_simulations_parser.add_argument(
        "--sampleVec",
        help = """
               Vector of input frequency files to work on
               (e.g.: 'commSAMP1____SAMP2.frq')
               """,
        required=True,
        nargs='+',
        type=argparse.FileType('r'),
        metavar="file.frq"
    )

    dist_simulations_parser.add_argument(
        "--numSimPairs",
        help="""
             Number of simulated pairs of individuals per class to be",
             generated.
              """,
        type=int,
        default=2000,
        metavar="n"
    )

    dist_simulations_parser.add_argument(
        "--freqFileHeader",
        help="When false a header will be added to the loaded .frq file",
        action='store_true'
    )


    # Build subparser for downsampleBam.R
    downsample_bam_parser = subparsers.add_parser(
        'downsampleBam',
        help="Downsize BAM files for faster processing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    downsample_bam_parser.add_argument(
        '--downsampleN',
        help="Maximum number of reads to downsample BAM files to",
        type=int,
        default=1500000,
        metavar="n"
    )

    downsample_bam_parser.add_argument(
        '--extensionBam',
        help="Work on BAM files with this extension/suffix",
        type=str,
        default=r"\\.bam$",
        metavar="regex"
    )

    downsample_bam_parser.add_argument(
        '--suffixDownBam',
        help="Suffix for new downsampled BAM files",
        type=str,
        default="_subsampled",
        metavar="str"
    )

    downsample_bam_parser.add_argument(
        '--downsampleSeed',
        help="Optionally set seeding for samtool's subsampling.",
        type=int,
        default=None,
        metavar="n"
    )

    # Build subparser for downsamplePed.R
    downsample_ped_parser = subparsers.add_parser(
        'downsamplePed',
        help='Downsample text-PLINK files (ped/map) for faster processing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    downsample_ped_parser.add_argument(
        '--downsampleN',
        help="Maximum number of SNPs to downsample text-PLINK dataset to",
        type=int,
        default=80000,
        metavar="n"
    )

    downsample_ped_parser.add_argument(
        '--extensionPed',
        help="Work on PED/MAP files with this extension/suffix",
        type=str,
        default=r"\\.ped$",
        metavar="regex"
    )

    downsample_ped_parser.add_argument(
        '--suffixDownPed',
        help="Extenstion/suffix for new downsampled files",
        type=str,
        default="_subsampled",
        metavar="str"
    )

    # Build subparser for individualisePlinks.R
    individualise_plinks_parser = subparsers.add_parser(
        'individualisePlinks',
        help="""
             Extract and convert all samples from a binary PLINK dataset into
             individual text-PLINK sets (.ped/.map).
             """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    individualise_plinks_parser.add_argument(
        '--dataset',
        help="Target dataset name without extensions",
        required=True,
        type=check_plink_basename,
        metavar="<basename>"
    )


    # Build subparser for simsForErrorRates.R
    sims_for_error_rates_parser = subparsers.add_parser(
        'simsForErrorRates',
        help="""
             Calculate error rates for specific subsets of SNPs by
             generating simulated individuals based on a provided allele
             frequencies file (.frq)
             """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    sims_for_error_rates_parser.add_argument(
        '--thins',
        help="List of SNP subset sizes",
        type=int,
        nargs="+",
        default=[100,500,1000,2000,3000,4000],
        metavar="n"
    )
    
    sims_for_error_rates_parser.add_argument(
        '--simSeqs',
        help="List of the number Number of simulated relationships",
        type=int,
        nargs="+",
        default=[100,250,500,1000,2000,4000],
        metavar="n"
    )
        
    sims_for_error_rates_parser.add_argument(
        '--freqsFi',
        help="Frequencies file in plink format with header (.frq)",
        type=argparse.FileType('r'),
        required=True,
        metavar="file.frq"
    )

    return parser.parse_args()




def main():
    args = build_parser()

    if args.helper == "individualisePlinks" and not is_available("plink"):
        raise RuntimeError(
            "'plink' is a runtime dependency of individualisePlinks, but "
            "cannot be found within your path. Please ensure this tool is "
            "properly installed on your system and available from your PATH."
        )
        
    subcommand = format_R_command(args)
    call_helper(subcommand)

if __name__ == '__main__':
    main()
