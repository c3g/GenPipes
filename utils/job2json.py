#!/usr/bin/env python

### Edouard Henrion (2017/07/18)
### job2json

import os
import sys
import getopt
import json

def getarg(argument):
    step_name = ""
    job_name = ""
    job_id = ""
    command = ""
    input_files = ""
    output_files = ""
    job_dependencies = ""
    job_log = ""
    job_done = ""
    json_files = ""
    success = True

    optli,arg = getopt.getopt(argument[1:], "s:n:i:c:f:o:d:l:b:j:g:h", ['step', 'job_name', 'job_id', 'command', 'infiles', 'outfiles', 'dependencies', 'job_log', 'job_done', 'json_outfiles', 'success', 'help'])

    if len(optli) == 0 :
        usage()
        sys.exit("Error : No argument given")

    for option, value in optli:
        if option in ("-s", "--step"):
            if str(value) == "" :
                sys.exit("Error - step (-s, --step) not provided...\n" + str(value))
            else :
                step_name = str(value)
        if option in ("-n", "--job_name"):
            if str(value) == "" :
                sys.exit("Error - job_name (-n, --job_name) not provided...\n")
            else :
                job_name = str(value)
        if option in ("-i", "--job_id"):
            if str(value) == "" :
                sys.exit("Error - job_id (-i, --job_id) not provided...\n")
            else :
                job_id = str(value)
        if option in ("-c", "--command"):
            if str(value) == "" :
                sys.exit("Error - command (-c, --command) not provided...\n")
            else :
                command = str(value)
        if option in ("-f", "--infiles"):
            if str(value) == "" :
                sys.exit("Error - infiles (-f, --infiles) not provided...\n")
            else :
                input_files = str(value)
        if option in ("-o", "--outfiles"):
            if str(value) == "" :
                sys.exit("Error - outfiles (-o, --outfiles) not provided...\n")
            else :
                output_files = str(value)
        if option in ("-d", "--dependencies"):
            job_dependencies = str(value)
        if option in ("-l", "--job_log"):
            if str(value) == "" :
                sys.exit("Error - job_log (-l, --job_log) not provided...\n")
            else :
                job_done = str(value)
        if option in ("-b", "--job_done"):
            if str(value) == "" :
                sys.exit("Error - job_done (-b, --job_done) not provided...\n")
            else :
                job_done = str(value)
        if option in ("-j", "--json_outfiles"):
            if str(value) == "" :
                sys.exit("Error - json_outfiles (-j, --json_outfiles) not provided...\n")
            else :
                json_files = str(value)
        if option in ("-g", "--success"):
            success = str(value)
        if option in ("-h", "--help"):
            usage()
            sys.exit()

    return step_name, job_name, job_id, command, input_files, output_files, job_dependencies, job_log, job_done, json_files, success

def usage():
    print "USAGE : job2json.py [option] "
    print "       -s    --step          : name of the step of the current job"
    print "       -n    --job_name      : name of the current job"
    print "       -i    --job_id        : ID of the current job"
    print "       -c    --command       : command line executed for the current job"
    print "       -f    --infiles       : comma-separated list of input files for the curent job"
    print "       -o    --outfiles      : comma-separated list of output files for the curent job"
    print "       -d    --dependencies  : comma-separated list of IDs of job dependencies for the current job"
    print "       -l    --job_log       : name of the log file for the current job"
    print "       -b    --job_done      : name of the done file for the current job"
    print "       -j    --json_outfiles : comma-separated list of names of json files which need to be appended affected by the current job"
    print "       -g    --success       : boolean value to indicate if the job has failed (False/0) or succeeded (True/1) - Default : True"
    print "       -h    --help          : this help \n"

def main():
    print "\n-------------------------------------------------------------------------------------------"
    print "job2json.py will append a JSON section describing a pipeline job that has just finished"
    print "to a JSON file which was pre-generated when the pipeline was launched."
    print "This script is usually launched automatically after each pipeline job has run successfully."
    print "This program was written by Edouard HENRION"
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "-------------------------------------------------------------------------------------------\n"

    print "command line used :\n" + " ".join(sys.argv[:])

    step_name, job_name, job_id, command, input_files, output_files, job_dependencies, job_log, job_done, json_files, success = getarg(sys.argv)

    for jfile in json_files.split(","):
        with open(jfile, 'r') as json_file:
            current_json = json.load(json_file)

        step_found = False
        job_found = False
        # check if the step already exists in the json object
        for jstep in current_json['sample']['pipeline']['step']:
            # if step exists already, then add the current job to it 
            if jstep['name'] == step_name:
                step_found = True
                # check if the job already exists in the json object
                for jjob in jstep['job']:
                    # if job exists already, then replace it with the current job
                    if jjob['name'] == job_name:
                        job_found = True
                        jjob['name'] = job_name
                        jjob['id'] = job_id
                        jjob['command'] = command
                        jjob['input_file'] = [input_files]
                        jjob['output_file'] = [output_files]
                        jjob['dependency'] = [job_dependencies]
                        jjob['log_file'] = job_log
                        if success
                            jjob['completion'] = "job successfully completed"
                            jjob['done_file'] = job_done
                        else:
                            jjob['completion'] = "job failed..."
                # if job does not exists already, add it to the current step
                if not job_found :
                    if success:
                        jstep['job'].append(
                            {
                                "name": job_name,
                                "id": job_id,
                                "command": command,
                                "input_file": [input_files],
                                "output_file": [output_files],
                                "dependency": [job_dependencies],
                                "log_file": job_log,
                                "done_file": job_done,
                                "completion": "job successfully completed"
                            }
                        )
                    else:
                        jstep['job'].append(
                            {
                                "name": job_name,
                                "id": job_id,
                                "command": command,
                                "input_file": [input_files],
                                "output_file": [output_files],
                                "dependency": [job_dependencies],
                                "log_file": job_log,
                                "completion": "job failed..."
                            }
        # if step does not exists already, create it as well as the current job
        if not step_found:
            if success:
                current_json['sample']['pipeline']['step'].append(
                    {
                        'name': step_name,
                        'job': [{
                            "name": job_name,
                            "id": job_id,
                            "command": command,
                            "input_file": [input_files],
                            "output_file": [output_files],
                            "dependency": [job_dependencies],
                            "log_file": job_log,
                            "done_file": job_done,
                            "completion": "job successfully completed"
                        }]
                    }
                )
            else:
                current_json['sample']['pipeline']['step'].append(
                    {
                        'name': step_name,
                        'job': [{
                            "name": job_name,
                            "id": job_id,
                            "command": command,
                            "input_file": [input_files],
                            "output_file": [output_files],
                            "dependency": [job_dependencies],
                            "log_file": job_log,
                            "completion": "job failed..."
                        }]
                    }
                )

        # Print to file
        with open(jfile, 'w') as out_json:
            json.dump(current_json, out_json, indent=4, sort_keys=True)

main()
