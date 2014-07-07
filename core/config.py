#!/usr/bin/env python

# Python Standard Modules
import ConfigParser
import logging
import os
import re
import subprocess

log = logging.getLogger(__name__)

class Config(ConfigParser.SafeConfigParser):

    def __init__(self):
        ConfigParser.SafeConfigParser.__init__(self)

    @property
    def filepath(self):
        return self._filepath

    def parse_file(self, config_file):
        self._filepath = config_file.name

        # Make option names case sensitive
        self.optionxform = str
        self.readfp(config_file)
        self.check_modules()

    # Check by a system call if all modules defined in config file are available
    def check_modules(self):
        modules = []

        # Retrieve all unique module version values in config file
        # assuming that all module key names contain "moduleVersion"
        for section in self.sections():
            for name, value in self.items(section):
                if re.search("^moduleVersion\.", name) and value not in modules:
                    modules.append(value)

        log.info("Check modules...")
        for module in modules:
            # Bash shell must be invoked in order to find "module" cmd
            module_show_output = subprocess.check_output(["bash", "-c", "module show " + module], stderr=subprocess.STDOUT)
            if re.search("Error", module_show_output, re.IGNORECASE):
                raise Exception("Error in config file with " + module + ":\n" + module_show_output)
            else:
                log.info("Module " + module + " OK")
        log.info("Module check finished\n")

    # Retrieve param in config file with optional definition check and type validation
    # By default, parameter is required to be defined in the config file
    def param(self, section, option, required=True, type='string'):
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
                elif type == 'list':
                    # Remove empty strings from list
                    return [x for x in self.get(section, option).split(",") if x]
                else:
                    return self.get(section, option)
            except:
                print("Error: parameter \"[" + section + "] " + option + "\" value \"" + self.get(section, option) + "\" is invalid!")
                raise
        elif required:
            raise Exception("Error: parameter \"[" + section + "] " + option + "\" is not defined in config file!")
        else:
            return ""

# Global config object used throughout the whole pipeline
config = Config()
