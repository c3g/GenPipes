################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes Pipelines.
#
# GenPipes Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import logging
import os
import re
import sys
import subprocess
import datetime

log = logging.getLogger(__name__)


def set_logger(log_level, sanity_check = False):
    """
    Set the logger level and format
    :param log_level: log level
    :param sanity_check: if True, log level is set to WARNING
    :return: None
    """
    min_level = getattr(logging, log_level.upper())
    if sanity_check:
        level = max(logging.WARNING, min_level)
        logging.basicConfig(stream=sys.stdout, level=level, format='%(message)s')
        log.warning("Sanity Check Report :")
    else:
        logging.basicConfig(level=min_level)

def time_to_datetime(time):
    """
    From transform time str to delta time:
    Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
    In fact it will also work for pbs/torque.
    Note that is strip from the string everithing that is not the time
    :param time: sting containing a slurm "--time" format substring
    :return: timedelta object
    """

    time = re.search("([0-9]+-)?[0-9]+:[0-9]+(:[0-9]+)?", time).group()

    if '-' in time:
        days, rest = time.split('-')
        rest = rest.split(':')[0]
    else:
        rest = time.split(':')
        days = 0

    hours = int(rest[0]) + int(days * 24)
    minutes = int(rest[1])

    if len(rest) > 2:
        sec = int(rest[2])
    else:
        sec = 0

    return datetime.timedelta(hours=hours, minutes=minutes, seconds=sec)


def expandvars(path, skip_escaped=False):
    """Expand environment variables of form $var and ${var}.
       If parameter 'skip_escaped' is True, all escaped variable references
       (i.e. preceded by backslashes) are skipped.
       Unknown variables are set to '' like in a bash shell.
    """
    def replace_var(m):
        return os.environ.get(m.group(2) or m.group(1), '')
    re_var = (r'(?<!\\)' if skip_escaped else '') + r'\$(\w+|\{([^}]*)\})'
    return re.sub(re_var, replace_var, path)


def container_wrapper_argparse(script, argv):
    """
    This method starts the pipeline inside the system put in place by the genpipes_in_a_container
    cvmfs/container wrapper if --wrap is present in the argv.

    Args:
        script: the tool that was called
        argv: Any valid Genpipes option.
    Returns: Return code of the subprocess or None if not --wrap option in the argv
    """

    if '--wrap' in argv:
        script_dir_current = os.path.dirname(os.path.realpath(__file__))
        parser = argparse.ArgumentParser(conflict_handler='resolve')
        parser.add_argument(
            '--wrap',
            type=str,
            help="Path to the genpipes cvmfs wrapper script. Default: genpipes/resources/container/bin/container_wrapper.sh",
            nargs='?'
            )
        parsed, argv = parser.parse_known_args(argv)
        if parsed.wrap is None:
            genpipes_home = '/'.join(script_dir_current.split('/')[:-2])
            default_wrapper = f'{genpipes_home}/resources/container/bin/container_wrapper.sh'
            parsed.wrap = default_wrapper

        wrap_option = ['--container', 'wrapper', parsed.wrap]

        if ('batch' in argv) and '--no-json' not in argv:
            parser.error('Combining --wrap  and -j "batch" options requires --no-json')

        sys.stderr.write('wrapping\n')
        sys.stderr.write(f"{parsed.wrap} {script} {' '.join(argv)} {' '.join(wrap_option)}\n")
        # Running the pipeline inside the container
        ret_code = subprocess.call([parsed.wrap,  script] + argv + wrap_option)
        sys.exit(ret_code)
    else:
        return None

def strtobool (val):
    """Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))
