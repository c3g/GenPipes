#!/usr/bin/env python

import re

class Step:
    def __init__(self, get_job_by_item, loop):
        # Step name is used in Bash $JOB_ID variable, hence only alphanumeric and "_" characters are allowed
        if re.search("^[a-zA-Z]\w+$", get_job_by_item.__name__):
            self._name = get_job_by_item.__name__
        else:
            raise Exception("Error: step name \"" + get_job_by_item.__name__ +
                "\" is invalid (should match [a-zA-Z][a-zA-Z0-9_]+)!")

        self._name = get_job_by_item.__name__
        self._loop = loop
        self._get_job_by_item = get_job_by_item
        self._jobs = []

    @property
    def name(self):
        return self._name

    @property
    def loop(self):
        return self._loop

    @property
    def get_job_by_item(self):
        return self._get_job_by_item

    @property
    def jobs(self):
        return self._jobs

    def job_id(self, job):
        return self.name + "_" + str(self.jobs.index(job) + 1) + "_JOB_ID"
