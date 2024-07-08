################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################
# Python Standard Modules
import configparser 
import glob
import logging
import os
import re
import subprocess
import io

from ..utils import utils

log = logging.getLogger(__name__)

class Config(configparser.ConfigParser):

    # True only for continuous integration testing
    continuous_integration_testing = 'GENPIPES_CIT' in os.environ
    # continuous_integration_testing = True
    cluster_walltime = "cluster_walltime"
    cit_options = [cluster_walltime]
    cit_prefix = 'cit_'

    # Setting the sanity-check flag to 'False'
    # Will be set to 'True' by the pipeline if it runs with the '--sanity-check' parameter
    sanity = False

    def __init__(self):
        configparser.ConfigParser.__init__(self)

    @property
    def filepath(self):
        return self._filepath

    def parse_files(self, config_files):
        # Make option names case sensitive
        self.optionxform = str
        for config_file in config_files:
            if isinstance(config_file,  io.IOBase):
                self.read_file(config_file)
            else:
                with open(config_file, 'r') as f:
                    self.read_file(f)
        self.check_modules()

    # Check by a system call if all modules defined in config files are available
    def check_modules(self):
        modules = []
        query_module = "show"
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in self.sections():
            for name, value in self.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)
                if re.search("^query_module", name):
                    query_module = value

        if self.sanity:
            log.warning("* Checking modules...")
        else:
            log.info("Check modules...")
        cmd_query_module = ["module", query_module]
        # for module in modules:

        cmd = ' '.join(cmd_query_module + modules)
        os.environ['MODULES_PAGER'] = ''
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        dout, derr = p.communicate()
        if p.returncode != 0 or b':ERROR:' in dout:
            _raise(SanitycheckError("Error in config file(s) with:\n{}".format(dout)))

        log.info("module check finished\n")

    # Retrieve param in config files with optional definition check and type validation
    # By default, parameter is required to be defined in one of the config file
    def global_get(self, section, option, required=True, param_type='string', **kwargs):
        # Store original section for future error message, in case 'DEFAULT' section is used eventually
        original_section = section
        # Keep that if block first, it is only evaluated in testing mode
        if self.continuous_integration_testing and option in self.cit_options:
            # hack because this class becomes a global
            try:
                return super().get(section, '{}{}'.format(self.cit_prefix, option))
            except configparser.Error as e:
                pass

            if option == self.cluster_walltime and self.has_section(section) and self.has_option(section, option):

                from_section = super().get(section, option)
                from_default = super().get('DEFAULT', option)
                if not utils.time_to_datetime(from_section):
                    return from_default
                elif utils.time_to_datetime(from_default) <= utils.time_to_datetime(from_section):
                    return from_default
                else:
                    return from_section

        if not self.has_section(section):
            section = 'DEFAULT'

        if self.has_option(section, option):
            try:
                if param_type == 'int':
                    return self.getint(section, option)
                elif param_type == 'posint':
                    value = self.getint(section, option)
                    if value > 0:
                        return value
                    else:
                        _raise(SanitycheckError("Integer \"" + str(value) + "\" is not > 0!"))
                elif param_type == 'float':
                    return self.getfloat(section, option)
                elif param_type == 'boolean':
                    return self.getboolean(section, option)
                elif param_type == 'filepath':
                    value = os.path.expandvars(super().get(section, option))
                    if not value and not required:
                        return None
                    if os.path.isfile(value):
                        return value
                    else:
                        log.debug(f"{required=}")
                        _raise(SanitycheckError(f"File path \"{value}\" provided in section [{section}] for option {option} does not exist or is not a valid regular file!"))
                elif param_type == 'dirpath':
                    value = os.path.expandvars(super().get(section, option))
                    if not value and not required:
                        return None
                    if os.path.isdir(value):
                        return value
                    else:
                        _raise(SanitycheckError(f"Directory path \"{value}\" provided in section [{section}] for option {option} does not exist or is not a valid directory!"))
                elif param_type == 'prefixpath':
                    value = os.path.expandvars(super().get(section, option))
                    if glob.glob(value + "*"):
                        return value
                    else:
                        _raise(SanitycheckError(f"Prefix path \"{value}\" provided in section [{section}] for option {option} does not match any file!"))
                elif param_type == 'list':
                    # Remove empty strings from list
                    return [x for x in super().get(section, option).split(",") if x]
                elif param_type == 'string':
                    return super().get(section, option)
                else:
                    _raise(SanitycheckError("Unknown parameter type '" + param_type + "'"))
            except Exception as e:
                _raise(SanitycheckError("Error found :\n  " + str(e) + f"\n{section=}, {option=}"))
        elif required:
            _raise(SanitycheckError("Error: REQUIRED parameter \"[" + original_section + "] " + option + "\" is not defined in config file(s)!"))
        # Returning empty string instead of None if nothing catched before as error
        return ""

    @staticmethod
    def config_error(message=''):
        _raise(SanitycheckError(message))

class Error(Exception):
    pass


class SanitycheckError(Error):
    pass


def _raise(error_obj):
    if Config.sanity:
        if isinstance(error_obj, SanitycheckError):
            log.error(error_obj)
    else:
        raise error_obj

# Global config object used throughout the whole pipeline
global_conf = Config()
