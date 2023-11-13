
import argparse
import os
import sys

from ... import utils

from ..ampliconseq import AmpliconSeq

# Append mugqic_pipelines directory to Python library path
# that is the old crappy setup!
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

def main(argv=None):

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

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    report = parsed_args.report
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file

    # Specific pipeline options
    protocol = parsed_args.protocol

    pipeline = AmpliconSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                           clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                           design_file=design_file, no_json=no_json, container=container,
                           protocol=protocol)

    pipeline.submit_jobs()


if __name__ == '__main__':
    main()
