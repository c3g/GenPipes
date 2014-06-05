#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules

def create_scheduler(type):
    if type == "torque":
        return TorqueScheduler()
    elif type == "batch":
        return BatchScheduler()
    elif type == "daemon":
        return DaemonScheduler()
    else:
        raise Exception("Error: scheduler type \"" + type + "\" is invalid!")

class Scheduler:
    def submit(self, jobs):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

class TorqueScheduler(Scheduler):
    def submit(self, jobs):
       for job in jobs:
            print("echo(\"\n" + job.command_with_modules + " \\\n | qsub depend=afterok:" +
                ":".join(["$" + dependency_job.id for dependency_job in job.dependency_jobs])) + "\")"

class BatchScheduler(Scheduler):
    def submit(self, jobs):
       for job in jobs:
            print(job.command_with_modules)

class DaemonScheduler(Scheduler):
    def submit(self, jobs):
       for job in jobs:
            print("daemon(" + job.command_with_modules + ")")
