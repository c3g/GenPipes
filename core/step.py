################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import re

class Step(object):
    def __init__(self, create_jobs, analyse_type=None):
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
        self._analyse_type = analyse_type

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def create_jobs(self):
        return self._create_jobs

    @create_jobs.setter
    def create_jobs(self, value):
        self._create_jobs = value

    @property
    def jobs(self):
        return self._jobs

    @jobs.setter
    def jobs(self, value):
        self._jobs = value

    @property
    def analyse_type(self):
        return self._analyse_type

    @analyse_type.setter
    def analyse_type(self, value):
        self._analyse_type = value

    def add_job(self, job):
        self.jobs.append(job)
        job.id = self.name + "_" + str(len(self.jobs)) + "_JOB_ID"
        
    def add_job(self, job):
        self.jobs.append(job)
        job.id = self.name + "_" + str(len(self.jobs)) + "_JOB_ID"