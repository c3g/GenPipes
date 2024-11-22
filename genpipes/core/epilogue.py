#!/usr/bin/env python3

import csv
import os
import subprocess
import logging
from io import StringIO
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO, format='EPILOGUE - %(message)s')

def get_slurm_job_info(job_id):
    try:
        logging.info(f"Fetching job info for job ID: {job_id}")
        result = subprocess.run(
            ["sacct", "-j", f"{job_id}", "--parsable", "--format=JobID,JobName,User,NodeList,Priority,Submit,Eligible,Timelimit,ReqCPUS,ReqMem,State,Start,End,Elapsed,AveCPU,AveRSS,MaxRSS,AveDiskRead,MaxDiskRead,AveDiskWrite,MaxDiskWrite"],
            capture_output=True, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error retrieving job info: {e}")
        return None

def parse_slurm_job_info(job_info, job_id):
    print(job_info)
    job_details = {}
    reader = csv.DictReader(StringIO(job_info), delimiter='|')
    for row in reader:
        if row['JobID'] == job_id:
            # Extracting from the first line (line starting with JobID)
            job_details['JobID'] = row['JobID']
            job_details['JobName'] = row['JobName']
            job_details['User'] = row['User']
            job_details['NodeList'] = row['NodeList']
            job_details['Priority'] = row['Priority']
            job_details['Submit'] = row['Submit']
            job_details['Eligible'] = row['Eligible']
            job_details['Timelimit'] = row['Timelimit']
            job_details['ReqCPUS'] = row['ReqCPUS']
            job_details['ReqMem'] = row['ReqMem']

        # Extracting from the third line (line starting with ${JOBID}.0)
        elif row['JobID'] == f"{job_id}.0":
            job_details['State'] = row['State']
            job_details['Start'] = row['Start']
            job_details['End'] = row['End']
            job_details['Elapsed'] = row['Elapsed']
            job_details['AveCPU'] = row['AveCPU']
            job_details['AveMem'] = row['AveRSS']
            job_details['MaxMem'] = row['MaxRSS']
            job_details['AveDiskRead'] = row['AveDiskRead']
            job_details['MaxDiskRead'] = row['MaxDiskRead']
            job_details['AveDiskWrite'] = row['AveDiskWrite']
            job_details['MaxDiskWrite'] = row['MaxDiskWrite']
            break
    return job_details

def get_pbs_job_info(job_id):
    try:
        logging.info(f"Fetching PBS job info for job ID: {job_id}")
        result = subprocess.run(
            ["qstat", "-f", f"{job_id}"],
            capture_output=True, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error retrieving PBS job info: {e}")
        return None

def parse_pbs_job_info(job_info):
    job_details = {}
    lines = job_info.strip().split('\n')
    for line in lines:
        if '=' in line:
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip()
            if key == 'Job Id':
                job_details['JobID'] = value
            elif key == 'Job_Name':
                job_details['JobName'] = value
            elif key == 'Job_Owner':
                job_details['User'] = value.split('@')[0]
            elif key == 'exec_host':
                job_details['NodeList'] = value
            elif key == 'Priority':
                job_details['Priority'] = value
            elif key == 'qtime':
                job_details['Submit'] = value
            elif key == 'etime':
                job_details['Eligible'] = value
            elif key == 'Resource_List.walltime':
                job_details['Timelimit'] = value
            elif key == 'Resource_List.ncpus':
                job_details['ReqCPUS'] = value
            elif key == 'Resource_List.mem':
                job_details['ReqMem'] = value
            elif key == 'job_state':
                job_details['State'] = value
            elif key == 'start_time':
                job_details['Start'] = value
            elif key == 'comp_time':
                job_details['End'] = value
            elif key == 'resources_used.walltime':
                job_details['Elapsed'] = value
            elif key == 'resources_used.cput':
                job_details['AveCPU'] = value
            elif key == 'resources_used.mem':
                job_details['AveMem'] = value
            elif key == 'resources_used.vmem':
                job_details['MaxMem'] = value
            elif key == 'resources_used.read_bytes':
                job_details['AveDiskRead'] = value
            elif key == 'resources_used.write_bytes':
                job_details['AveDiskWrite'] = value
    return job_details

def time_str_to_seconds(time_str):
    h, m, s = map(int, time_str.split(':'))
    return h * 3600 + m * 60 + s

def calculate_time_difference(start_time, end_time):
    start_dt = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S")
    end_dt = datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S")
    delta = end_dt - start_dt
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{days:02}:{hours:02}:{minutes:02}:{seconds:02}"

def calculate_time_efficency(elapsed, timelimit):
    elapsed_seconds = sum(int(x) * 60 ** i for i, x in enumerate(reversed(elapsed.split(':'))))
    timelimit_seconds = sum(int(x) * 60 ** i for i, x in enumerate(reversed(timelimit.split(':'))))
    return 100 * elapsed_seconds / timelimit_seconds if timelimit_seconds > 0 else 0

def calculate_percentage(used, requested):
    return (used / requested * 100) if requested > 0 else 0

def convert_memory_to_gb(memory_str):
    if memory_str.endswith('K'):
        return float(memory_str[:-1]) / (1024 ** 2)
    elif memory_str.endswith('M'):
        return float(memory_str[:-1]) / 1024
    elif memory_str.endswith('G'):
        return float(memory_str[:-1])
    elif memory_str.endswith('T'):
        return float(memory_str[:-1]) * 1024
    return float(memory_str)

def main():
    job_id = os.getenv('SLURM_JOB_ID') or os.getenv('PBS_JOBID')
    if not job_id:
        logging.error("Unknown scheduler")
        return

    if 'SLURM_JOB_ID' in os.environ:
        job_info = get_slurm_job_info(job_id)
        if job_info:
            job_details = parse_slurm_job_info(job_info, job_id)
        else:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    elif 'PBS_JOBID' in os.environ:
        job_info = get_pbs_job_info(job_id)
        if job_info:
            job_details = parse_pbs_job_info(job_info)
        else:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return

    # Calculate time spent in queue format DD:HH:MM:SS
    time_in_queue = calculate_time_difference(job_details['Submit'], job_details['Start'])
    # Calculate time efficency between walltime and time used
    time_efficency = calculate_time_efficency(job_details['Elapsed'], job_details['Timelimit'])
    # Convert memory to GB
    req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])
    ave_mem_gb = convert_memory_to_gb(job_details['AveMem'])
    max_mem_gb = convert_memory_to_gb(job_details['MaxMem'])
    # Calculate percentages for CPU and memory usage
    req_cpus = time_str_to_seconds(job_details.get('ReqCPUS', '00:00:00'))
    ave_cpu = time_str_to_seconds(job_details.get('AveCPU', '00:00:00'))
    cpu_usage_percentage = calculate_percentage(ave_cpu, req_cpus)
    mem_usage_percentage = calculate_percentage(ave_mem_gb, req_mem_gb)

    print("-" * 90)
    logging.info(f"GenPipes Epilogue for job {job_id}")
    logging.info(f"Job name: {job_details.get('JobName', 'Unknown')}")
    logging.info(f"User: {job_details.get('User', 'Unknown')}")
    logging.info(f"Node(s): {job_details.get('NodeList', 'Unknown')}")
    logging.info(f"Priority: {job_details.get('Priority', 'Unknown')}")
    logging.info(f"Status: {job_details.get('State', 'Unknown')}")
    logging.info(f"Submit Time: {job_details.get('Submit', 'Unknown')}")
    logging.info(f"Eligible Time: {job_details.get('Eligible', 'Unknown')}")
    logging.info(f"Start Time: {job_details.get('Start', 'Unknown')}")
    logging.info(f"Time Spent in Queue (DD:HH:MM:SS): {time_in_queue}")
    logging.info(f"End Time: {job_details.get('End', 'Unknown')}")
    logging.info(f"Runtime: {job_details.get('Elapsed', 'Unknown')}")
    logging.info(f"Time Limit: {job_details.get('Timelimit', 'Unknown')}")
    logging.info(f"Time efficiency (Percentage of Walltime to Time used): {time_efficency:.1f}")
    logging.info(f"Number of CPU(s) requested: {job_details.get('ReqCPUS', 'Unknown')}")
    logging.info(f"Average CPU usage: {job_details.get('AveCPU', 'Unknown')}")
    logging.info(f"CPU efficiency (Percentage of CPU requested to CPU used): {cpu_usage_percentage:.1f}")
    logging.info(f"Memory Requested: {req_mem_gb:.2f} GB")
    logging.info(f"Average memory usage: {ave_mem_gb:.2f} GB")
    logging.info(f"Max memory usage: {max_mem_gb:.2f} GB")
    logging.info(f"Memory efficiency (Percentage of Memory requested to Memory used in average): {mem_usage_percentage:.1f}")
    logging.info(f"Average Disk Read: {job_details.get('AveDiskRead', 'Unknown')}")
    logging.info(f"Max Disk Read: {job_details.get('MaxDiskRead', 'Unknown')}")
    logging.info(f"Average Disk Write: {job_details.get('AveDiskWrite', 'Unknown')}")
    logging.info(f"Max Disk Write: {job_details.get('MaxDiskWrite', 'Unknown')}")
    print("-" * 90)

if __name__ == "__main__":
    main()
