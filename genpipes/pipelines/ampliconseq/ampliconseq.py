
import argparse
import sys

from ... import utils

from ..ampliconseq import AmpliconSeq

def main(argv=None):
    """
    Stub so people can still call the pipeline file as a script
    """

    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = AmpliconSeq.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = AmpliconSeq.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

if __name__ == '__main__':
    main()
