#!/usr/bin/env python

from Config import *

class Command:

    def __init__(self, name, config, input_files, output_files, module_entries = []):
        if re.search("^[a-zA-Z]\w*$", name):
            self._name = name
        else:
            raise Exception("Error: command name \"" + name +
                "\" is invalid (should match [a-zA-Z][a-zA-Z0-9_]*)!")

        self._config = config
        self._input_files = input_files
        self._output_files = output_files

        self._modules = []
        for section, option in module_entries:
            module = self.config.param(section, option)
            if module not in self.modules:
                self.modules.append(module)
            

    def show(self):
        print 'Basic -- name: ' + self._name + ", input_files: " + \
            ", ".join([input_file.name for input_file in self._input_files])

    @property
    def name(self):
        return self._name

    @property
    def config(self):
        return self._config

    @property
    def input_files(self):
        return self._input_files

    @property
    def output_files(self):
        return self._output_files

    @property
    def modules(self):
        return self._modules

    @property
    def line(self):
        return self._line

    @property
    def line_with_modules(self):
        return "module load " + " ".join(self.modules) + "\n" + self.line

    @property
    def is_up2date(self):
        latest_input_time = max([os.stat(os.path.expandvars(input_file)).st_mtime if os.path.isfile(os.path.expandvars(input_file)) else 0 for input_file in self.input_files])
        return latest_input_time

config = Config("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/dnaSeq.abacus.ini")
command = Command(
    "my_cmd",
    config,
    ["input1.bam", "/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/dnaSeq.abacus.ini", "/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/project.nanuq.csv"],
    ["output1.bam"],
    [
        ["trim", "moduleVersion.trimmomatic"],
        ["aln", "moduleVersion.bwa"],
        ["recalibration", "moduleVersion.picard"]
    ])

command.line = "ls -l"
print command.line_with_modules
print command.is_up2date
