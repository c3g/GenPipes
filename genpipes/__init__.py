"""
The top-level Genpipes command for bash interface
"""

import argparse
import re
import os
import sys
import importlib
from types import SimpleNamespace
import textwrap

from .utils import container_wrapper_argparse
from .version import __version__

PIPELINES = [
"AmpliconSeq",
"ChipSeq",
"CoVSeq",
"DnaSeq",
"MethylSeq",
"Nanopore",
"Nanopore_CoVSeq",
"RnaSeq",
"RnaSeq_DeNovo_Assembly",
"RnaSeq_Light"
]

def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """
    class DefaultCommand(object):
        @staticmethod
        def main(self):
            parser.print_help()
            return 0

    parser.set_defaults(__command__ = DefaultCommand)



def make_parser(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    container_wrapper_argparse("genpipes", argv)

    parser = argparse.ArgumentParser(
        prog        = "genpipes",
        description = "GenPipes consists of Python scripts which create a list of jobs running Bash commands. "
                      "Those scripts support dependencies between jobs and smart restart mechanism if some jobs "
                      "fail during pipeline execution. Jobs can be submitted in different ways: by being sent to "
                      "a PBS or a SLURM scheduler or by being run as a series of commands in batch through a Bash "
                      "script. Job commands and parameters can be modified through several configuration files.")
    
    parser.add_argument("-v", "--version", action="version", version='%(prog)s {version}'.format(version=__version__))

    subparsers = parser.add_subparsers()

    add_default_command(parser)

    for pipeline in PIPELINES:
        # Add a subparser for each command.
        module = pipeline.lower()
        classname = pipeline.replace('_','')
        p_name = classname.lower()
        p_module = importlib.import_module('.pipelines.' + module, package='genpipes')
        p_class = getattr(p_module, classname)

        epilog = p_class.process_help(argv)

        subparser = subparsers.add_parser(module,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            conflict_handler='resolve', epilog=epilog, help='\n'.join(p_class.__doc__.split()[0:2]))

        subparser.set_defaults(__command__ = p_module)

        # Let the class feed its subparser.
        subparser = p_class.argparser(subparser)

    return parser

def main(argv):
    parser = make_parser(argv)
    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)


    return parsed_args.__command__.main(parsed_args)
