#!/usr/bin/env python

from Step import *

class Pipeline:
    def __init__(self, name, config, steps, step_range):
        self._name = name
        self._config = config
        self._steps = []
        for step in steps:
            if self.get_steps_by_name(step["name"].__name__):
                raise Exception("Error: step name \"" + step["name"].__name__ + "\" is not unique for this pipeline!")
            else:
                self.steps.append(Step(step["name"], step["loop"]))

        self._step_range = step_range

    @property
    def name(self):
        return self._name

    @property
    def config(self):
        return self._config

    @property
    def steps(self):
        return self._steps

    @property
    def step_range(self):
        return self._step_range

    def get_steps_by_name(self, name):
        return [step for step in self.steps if step.name == name]

    def get_steps_by_range(self):
        return [self.steps[i - 1] for i in parse_range(self.step_range)]

    def get_dependency_jobs(self, current_job):
        dependency_jobs = []
        for step in self.get_steps_by_range():
            for step_job in step.jobs:
                if list(set(current_job.input_files) & set(step_job.output_files)):
                    dependency_jobs.append(step_job)
        return dependency_jobs

    def show(self):
        print 'Pipeline: %s' % self.name
        print '-- steps: '
        for step in self.get_steps_by_range():
            print "Step: " + step.name
            for item in step.loop:
                job = step.get_job_by_item(item)
                step.jobs.append(job)
                print step.job_id(job)
                for dependency_job in self.get_dependency_jobs(job):
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
GLOBAL = [1]
