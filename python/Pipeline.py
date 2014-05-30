#!/usr/bin/env python

import argparse

from Step import *

class Pipeline:
    def __init__(self, config, step_dict_map, step_range):
        self._config = config
        self._step_map = []
        for step_dict in step_dict_map:
            step_name = step_dict["name"].__name__
            if self.steps_by_name(step_name):
                raise Exception("Error: step name \"" + step_name + "\" is not unique for this pipeline!")
            else:
                self.step_map.append(Step(step_dict["name"], step_dict["loop"]))

        self._steps = self.steps_by_range(step_range)

        self._parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=self.numbered_step_list())
        self._parser.add_argument("-c", "--config", help="config INI-style file", type=file)
        self._parser.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'")
        self._parser.add_argument("-i", "--input", help="input readset file", type=file)
        self._parser.add_argument("-o", "--output", help="output directory (default: current)")

    @property
    def config(self):
        return self._config

    @property
    def step_map(self):
        return self._step_map

    @property
    def steps(self):
        return self._steps

    @property
    def parser(self):
        return self._parser

    def steps_by_name(self, name):
        return [step for step in self.step_map if step.name == name]

    def steps_by_range(self, range):
        return [self.step_map[i - 1] for i in parse_range(range)]

    # Return numbered list of steps
    def numbered_step_list(self):
        return "Steps:\n" + "\n".join([str(idx + 1) + "- " + step.name for idx, step in enumerate(self.step_map)])

    def dependency_jobs(self, current_job):
        dependency_jobs = []
        for step in self.steps:
            for step_job in step.jobs:
                if list(set(current_job.input_files) & set(step_job.output_files)):
                    dependency_jobs.append(step_job)
        return dependency_jobs

    def show(self):
        print 'Pipeline: %s' % self.__class__.__name__
        print '-- steps: '
        for step in self.steps:
            print "Step: " + step.name
            for item in step.loop:
                if item:
                    job = step.create_job(item)
                else:
                    job = step.create_job()
                step.add_job(job)
                print job.id
                for dependency_job in self.dependency_jobs(job):
                    print "  Dependency: " + dependency_job.command
                print "  Input files:" + ",".join(job.input_files)
                print "  " + job.command
                print "  Output files:" + ",".join(job.output_files)
                print

# Return a range list given a string.
# e.g. parse_range('1,3,5-12') returns [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]
def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

# Define pipeline step global loop as 1-element array for practical 1-time iteration
GLOBAL = [None]
