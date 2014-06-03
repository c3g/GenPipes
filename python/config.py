#!/usr/bin/env python

# Python Standard Modules
import os
import re
import subprocess
import ConfigParser

class Config(ConfigParser.SafeConfigParser):

    def __init__(self, config_file):
        ConfigParser.SafeConfigParser.__init__(self)
        # Make option names case sensitive
        self.optionxform = str
        self.readfp(open(config_file))
        self.check_modules()

    # Check by a system call if all modules defined in config file are available
    def check_modules(self):
        modules = []

        # Retrieve all unique module version values in config file
        # assuming that all module key names contain "moduleVersion"
        for section in self.sections():
            for name, value in self.items(section):
                if re.search("moduleVersion", name) and value not in modules:
                    modules.append(value)

        print "Check modules..."
        for module in modules:
            # Bash shell must be invoked in order to find "module" cmd
            module_show_output = subprocess.check_output(["bash", "-c", "module show " + module], stderr=subprocess.STDOUT)
            if re.search("Error", module_show_output, re.IGNORECASE):
                raise Exception("Error in config file with " + module + ":\n" + module_show_output)
            else:
                print "Module " + module + " OK"
        print

    # Retrieve param in config file with optional definition check and type validation
    # By default, parameter is required to be defined in the config file
    def param(self, section, option, required=True, type='string'):
        if self.has_option(section, option):
            try:
                if type == 'int':
                    return self.getint(section, option)
                elif type == 'float':
                    return self.getfloat(section, option)
                elif type == 'boolean':
                    return self.getboolean(section, option)
                elif type == 'float':
                    return self.getfloat(section, option)
                elif type == 'filepath':
                    value = self.get(section, option)
                    if os.path.isfile(os.path.expandvars(value)):
                        return value
                    else:
                        raise Exception("File path \"" + value + "\" does not exist or is not a valid regular file!")
                elif type == 'dirpath':
                    value = self.get(section, option)
                    if os.path.isdir(os.path.expandvars(value)):
                        return value
                    else:
                        raise Exception("Directory path \"" + value + "\" does not exist or is not a valid directory!")
                elif type == 'list':
                    # Remove empty strings from list
                    return [x for x in self.get(section, option).split(",") if x]
                else:
                    return self.get(section, option)
            except:
                print "Error: parameter \"[" + section + "] " + option + "\" value \"" + self.get(section, option) + "\" is invalid!"
                raise
        elif required:
            raise Exception("Error: parameter \"[" + section + "] " + option + "\" is not defined in config file!")
        else:
            return ""

#config = Config("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/dnaSeq.abacus.ini")
#print config.get("trim", "moduleVersion.java")
#print config.param("trim", "toto", False, "list")
