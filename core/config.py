#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import ConfigParser
import glob
import logging
import os
import re
import subprocess
import sys

log = logging.getLogger(__name__)

class Config(ConfigParser.SafeConfigParser):

    def __init__(self):
        ConfigParser.SafeConfigParser.__init__(self)

    @property
    def filepath(self):
        return self._filepath

    def parse_files(self, config_files):
        # Make option names case sensitive
        self.optionxform = str
        for config_file in config_files:
            self.readfp(config_file)
        self.check_modules()

    # Check by a system call if all modules defined in config files are available
    def check_modules(self):
        modules = []

        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in self.sections():
            for name, value in self.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        log.info("Check modules...")
        for module in modules:
            # Bash shell must be invoked in order to find "module" cmd
            module_show_output = subprocess.check_output(["bash", "-c", "module show " + module], stderr=subprocess.STDOUT)
            if re.search("Error", module_show_output, re.IGNORECASE):
                raise Exception("Error in config file(s) with " + module + ":\n" + module_show_output)
            else:
                log.info("Module " + module + " OK")
        log.info("Module check finished\n")

    # Retrieve param in config files with optional definition check and type validation
    # By default, parameter is required to be defined in one of the config file
    def param(self, section, option, required=True, type='string'):
        # Store original section for future error message, in case 'DEFAULT' section is used eventually
        original_section = section

        if not self.has_section(section):
            section = 'DEFAULT'

        if self.has_option(section, option):
            try:
                if type == 'int':
                    return self.getint(section, option)
                elif type == 'posint':
                    value = self.getint(section, option)
                    if value > 0:
                        return value
                    else:
                        raise Exception("Integer \"" + str(value) + "\" is not > 0!")
                elif type == 'float':
                    return self.getfloat(section, option)
                elif type == 'boolean':
                    return self.getboolean(section, option)
                elif type == 'filepath':
                    value = os.path.expandvars(self.get(section, option))
                    if os.path.isfile(value):
                        return value
                    else:
                        raise Exception("File path \"" + value + "\" does not exist or is not a valid regular file!")
                elif type == 'dirpath':
                    value = os.path.expandvars(self.get(section, option))
                    if os.path.isdir(value):
                        return value
                    else:
                        raise Exception("Directory path \"" + value + "\" does not exist or is not a valid directory!")
                elif type == 'prefixpath':
                    value = os.path.expandvars(self.get(section, option))
                    if glob.glob(value + "*"):
                        return value
                    else:
                        raise Exception("Prefix path \"" + value + "\" does not match any file!")
                elif type == 'list':
                    # Remove empty strings from list
                    return [x for x in self.get(section, option).split(",") if x]
                elif type == 'string':
                    return self.get(section, option)
                else:
                    raise Exception("Unknown parameter type '" + type + "'")
            except Exception as e:
                raise Exception("Error: parameter \"[" + section + "] " + option + "\" value \"" + self.get(section, option) + "\" is invalid!\n" + e.message)
        elif required:
            raise Exception("Error: parameter \"[" + original_section + "] " + option + "\" is not defined in config file(s)!")
        else:
            return ""

# Global config object used throughout the whole pipeline
config = Config()
