#!/usr/bin/env python3
import argparse
import csv
import datetime
import logging
import os
import re
import sys

logger = logging.getLogger(__name__)

class JobStat(object):

    JOBID = 'JobId'
    JOBSTATE = 'JobState'
    RUNTIME = 'RunTime'
    SUBMITTIME = 'SubmitTime'
    STARTTIME = 'StartTime'
    NUMCPUS = 'NumCPUs'
    MEM = 'mem'
    REQUEUE = 'Requeue'
    RESTARTS = 'Restarts'


    def __init__(self, path, number=None, depenendency=None, name=None):
        self._path = path
        self._number = number
        self._dependency = depenendency
        self._name = name
        self._prologue = {}
        self._epiloque = {}
        self._delta = None
        self._id = None
        self._start = None
        self._stop = None

        self.fill_from_file(self._path)

    def fill_from_file(self, path):
        # parse stdout and get all the prologue and epilogues values
        slurm_stdtout = open(path).read()
        bidon = re.findall(r"(\w+)=(\S+)\s", slurm_stdtout)
        all_value = {}
        for k, v in bidon:
            all_value.setdefault(k, []).append(v)

        try:
            pro, epi = [i for i, x in enumerate(all_value[self.JOBID])
                        if self.number == int(x)]
        except ValueError:
            logger.error('{} log is not in {}, there is a mismatch '
                         'between the job number and the output'
                         .format(self.number, self._path))

            # Get the closest number
            diff = [abs(int(v) - self.number) for v in all_value[self.JOBID]]
            self._id = int(all_value[self.JOBID][diff.index(min(diff))])

            pro, epi = [i for i, x in enumerate(all_value[self.JOBID])
                        if self._id == int(x)]

        for k, v in all_value.items():
             try :
                self._prologue[k] = v[pro]
                self._epiloque[k] = v[epi]
             except IndexError:
                logger.warning('{} = {} in not a slurm prologue/epilogue value or is ambiguous'.format(k, v))

        tres = re.findall(r"TRES=cpu=(\d+),mem=(\d+\w)", slurm_stdtout)
        self._prologue[self.NUMCPUS] = tres[pro][0]
        self._epiloque[self.NUMCPUS] = tres[epi][0]
        self._prologue[self.MEM] = tres[pro][1]
        self._epiloque[self.MEM] = tres[epi][1]



    @property
    def runtime(self):

        if self._delta is None:
            stop = [int(s) for s in self._epiloque[self.RUNTIME].split(':')]
            start = [int(s) for s in self._prologue[self.RUNTIME].split(':')]
            self._stop = datetime.timedelta(hours=stop[0], minutes=stop[1], seconds=stop[2])
            self._start = datetime.timedelta(hours=start[0], minutes=start[1], seconds=start[2])
            self._delta = self._stop - self._start

        return self._delta


    @property
    def prologue_time(self):
        return self.start_time + self._start

    
    @property
    def epilogue_time(self):
        return self.start_time + self._stop


    @property
    def start_time(self):
        return datetime.datetime.strptime(self._prologue[self.STARTTIME], '%Y-%m-%dT%H:%M:%S')


    @property
    def n_cpu(self):
        return int(self._prologue[self.NUMCPUS])

    @property
    def mem(self):
        return self._prologue[self.MEM]

    @property
    def name(self):
        return self._name

    @property
    def number(self):
        return int(self._number)

def get_report(job_list_tsv=None):


    job_output_path = os.path.dirname(job_list_tsv)

    with open(job_list_tsv) as tsvin:
        jobs = csv.reader(tsvin, delimiter='\t')

        stats = []
        i = 0
        for job in jobs:

            stats.append(JobStat(path=os.path.join(job_output_path, job[3]),
                                 name=job[1], number=job[0], depenendency=job[2]))
            # j = stats[-1]
            # print ('{} {} {} {} {} {}'.format(j.number, j._id , j.name, j.n_cpu, j.mem, j.runtime))
            # i = i + 1

    return stats


def print_report(stats):
    min_start = []
    max_stop = []
    total_machine = datetime.timedelta(0)
    total_machine_core = datetime.timedelta(0)
    for j in stats:
        print ('{} {} {} {} {} {}'.format(j.number, j._id, j.name, j.n_cpu, j.mem, j.runtime))

        min_start.append(j.prologue_time)
        max_stop.append(j.epilogue_time)
        total_machine += j.runtime
        total_machine_core += j.runtime * j.n_cpu

    min_start = min(min_start)
    max_stop = max(max_stop)

    # total on node
    total_humain = max_stop - min_start

    print ("Humain time: {}\n Cumulative machine time: {}\n"
           " cumulative core time: {}".format(total_humain, total_machine,
                                              total_machine_core))

    # keys = stats[0].keys()
    # with open('output.csv', 'wb') as output_file:
    #     dict_writer = csv.DictWriter(output_file, keys)
    #     dict_writer.writeheader()
    #     dict_writer.writerows(toCSV)

if __name__ == '__main__':



        logging.basicConfig(level=logging.INFO)

        parser = argparse.ArgumentParser()
        parser.add_argument('path', help="Path to a Genpipes job list")

        args = parser.parse_args()

        path = args.path

        stats = get_report(path)

        print_report(stats)
