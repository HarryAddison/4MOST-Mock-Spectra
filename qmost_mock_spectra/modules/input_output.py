import argparse
import textwrap
import yaml


def parser():
    '''
    Function that allows the user to parse a initfile path via the
    command line.
    '''

    help_msg = "A help message and description of the code"

    parser = argparse.ArgumentParser(prog="TiDES Spectral Analysis",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent(help_msg))

    parser.add_argument("-i",
                        "--initfile",
                        default=None,
                        help=("yaml file containing input frames "
                              "and other procedure arguments"))
    return parser


def load_config():
    '''
    Obtain the config parameters from the yaml file that was parsed via
    the command line.
    '''

    arg_parser = parser()
    args = arg_parser.parse_args()
    #print(f"\nInput parameter file: {args.initfile}\n\n")

    # If no parameter file was provided print help message
    if args.initfile is None:
        arg_parser.print_help()
        exit()

    # Open the initfile and read in the parameters
    with open(str(args.initfile), "r") as file:
        config = yaml.safe_load(file)

    return config
