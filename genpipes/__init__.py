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
import shtab

from .tools import tools

from .core.utils import container_wrapper_argparse, set_logger
from .__version__ import __version__

PIPELINES = {
    "AmpliconSeq": "AmpliconSeq pipeline",
    "ChipSeq": "ChipSeq pipeline",
    "CoVSeq": "CoVSeq pipeline",
    "DnaSeq": "DnaSeq pipeline",
    "LongRead_DnaSeq": "LongRead DnaSeq pipeline",
    "MethylSeq": "MethylSeq pipeline",
    "Nanopore_CoVSeq": "Nanopore CoVSeq pipeline",
    "RnaSeq": "RnaSeq pipeline",
    "RnaSeq_DeNovo_Assembly": "RnaSeq DeNovo Assembly pipeline",
    "RnaSeq_Light": "RnaSeq Light pipeline",
}

def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """
    class DefaultCommand():
        """
        The default command when no subcommand is provided.
        """
        @staticmethod
        def main(parsed_args):
            """
            The main function for the default command.
            """
            parser.print_help()
            return 0

    parser.set_defaults(__command__ = DefaultCommand)

def make_parser(argv=None):
    """
    Creates the GenPipes command-line parser.
    """
    if argv is None:
        argv = sys.argv[1:]

    container_wrapper_argparse("genpipes", argv)

    parser = argparse.ArgumentParser(
        prog="genpipes",
        description="GenPipes consists of Python scripts which create a list of jobs running Bash commands. Those scripts support dependencies between jobs and smart restart mechanism if some jobs fail during pipeline execution. Jobs can be submitted in different ways: by being sent to a PBS or a SLURM scheduler or by being run as a series of commands in batch through a Bash script. Job commands and parameters can be modified through several configuration files."
    )

    parser.add_argument("-v", "--version", action="version", version=f'{__version__}')
    shtab.add_argument_to(parser, ["-s", "--print-completion"])

    subparsers = parser.add_subparsers(dest='command')

    add_default_command(parser)

    for pipeline, description in PIPELINES.items():
        # Add a subparser for each command.
        module = pipeline.lower()
        classname = pipeline.replace('_', '')
        p_module = importlib.import_module('.pipelines.' + module, package='genpipes')
        p_class = getattr(p_module, classname)

        epilog = p_class.process_help(argv)

        subparser = subparsers.add_parser(
            module,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            conflict_handler='resolve',
            epilog=epilog,
            help=description,
        )

        subparser.set_defaults(__command__=p_module)

        # Let the class feed its subparser.
        subparser = p_class.argparser(subparser)

    # Create the parser for the "tools" command
    parser_tools = subparsers.add_parser('tools', help='GenPipes companion tools')
    tools.add_subcommands(parser_tools)

    return parser

def main(argv):
    """
    The main function for the GenPipes command.
    """
    parser = make_parser(argv)
    parsed_args = parser.parse_args(argv)

    if parsed_args.command == 'tools':
        if not hasattr(parsed_args, 'func'):
            subparser = parser._subparsers._group_actions[0].choices['tools']
            subparser.print_help()
            return
        if hasattr(parsed_args, 'func'):
            return parsed_args.func(parsed_args)


    sanity_check = getattr(parsed_args, 'sanity_check', None)
    loglevel = getattr(parsed_args, 'log', None)
    if loglevel:
        set_logger(loglevel, sanity_check=sanity_check)

    return parsed_args.__command__.main(parsed_args)
