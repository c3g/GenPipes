#!/usr/bin/env python

# Python Standard Modules
import re

class Step:
    def __init__(self, create_jobs):
        # Step name is used in Bash $JOB_ID variable, hence only alphanumeric and "_" characters are allowed
        step_name = create_jobs.__name__
        if re.search("^[a-zA-Z]\w+$", step_name):
            self._name = step_name
        else:
            raise Exception("Error: step name \"" + step_name +
                "\" is invalid (should match [a-zA-Z][a-zA-Z0-9_]+)!")

        self._name = step_name
        self._create_jobs = create_jobs
        self._jobs = []

    @property
    def name(self):
        return self._name

    @property
    def create_jobs(self):
        return self._create_jobs

    @property
    def jobs(self):
        return self._jobs

    def add_job(self, job):
        self.jobs.append(job)
        job.id = self.name + "_" + str(len(self.jobs)) + "_JOB_ID"
