#!/usr/bin/env python3

import csv
import os
import subprocess
import logging
import time
from io import StringIO

# Configure logging
logging.basicConfig(level=logging.INFO, format='PROLOGUE - %(message)s')

def get_slurm_job_info(job_id, retries=10, delay=5):
    """
    Retrieve job information from SLURM using sacct command with retries.
    """
    for _ in range(retries):
        try:
            result = subprocess.run(
                ["sacct", "-j", f"{job_id}", "--parsable", "--format=JobID,JobName,User,NodeList,Priority,Submit,Timelimit,ReqCPUS,ReqMem"],
                capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                job_info = result.stdout
                job_details = parse_slurm_job_info(job_info, job_id)
                if job_details:
                    if job_details['AveRSS']:
                        return job_details
                time.sleep(delay)
            else:
                time.sleep(delay)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error retrieving job info: {e}")
            return None
    logging.error(f"{job_info}")
    logging.error(f"Failed to retrieve complete job info for job {job_id} after {retries} attempts")
    return None

def parse_slurm_job_info(job_info, job_id):
    """
    Parse the job information retrieved from SLURM.
    """
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
            job_details['Timelimit'] = row['Timelimit']
            job_details['ReqCPUS'] = row['ReqCPUS']
            job_details['ReqMem'] = row['ReqMem']
    # Check if all necessary fields are populated
    required_fields = ['JobID', 'JobName', 'User', 'NodeList', 'Priority', 'Submit', 'Timelimit', 'ReqCPUS', 'ReqMem',]
    missing_fields = [field for field in required_fields if field not in job_details]
    if missing_fields:
        logging.warning(f"Missing fields: {', '.join(missing_fields)}")
        return None
    return job_details

def get_pbs_job_info(job_id):
    """
    Retrieve job information from PBS using qstat command
    """
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
    """
    Parse the job information retrieved from PBS.
    """
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
            elif key == 'Resource_List.walltime':
                job_details['Timelimit'] = value
            elif key == 'Resource_List.ncpus':
                job_details['ReqCPUS'] = value
            elif key == 'Resource_List.mem':
                job_details['ReqMem'] = value
    return job_details

def convert_memory_to_gb(memory_str):
    """
    Convert memory string to GB.
    """
    if memory_str.endswith('K'):
        return float(memory_str[:-1]) / (1024 ** 2)
    if memory_str.endswith('M'):
        return float(memory_str[:-1]) / 1024
    if memory_str.endswith('G'):
        return float(memory_str[:-1])
    if memory_str.endswith('T'):
        return float(memory_str[:-1]) * 1024
    return float(memory_str)

def main():
    """
    Main function to run the epilogue script.
    """
    job_id = os.getenv('SLURM_JOB_ID') or os.getenv('PBS_JOBID')
    if not job_id:
        logging.error("Unknown scheduler")
        return

    if 'SLURM_JOB_ID' in os.environ:
        job_details = get_slurm_job_info(job_id)
        if not job_details:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    elif 'PBS_JOBID' in os.environ:
        job_info = get_pbs_job_info(job_id)
        if job_info:
            job_details = parse_pbs_job_info(job_info)
        else:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return

    # Convert memory to GB
    req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])

    logging.info(f"GenPipes Prologue for job {job_id}")
    logging.info(f"Job name:                                                         {job_details.get('JobName', 'Unknown')}")
    logging.info(f"User:                                                             {job_details.get('User', 'Unknown')}")
    logging.info(f"Node(s):                                                          {job_details.get('NodeList', 'Unknown')}")
    logging.info(f"Priority:                                                         {job_details.get('Priority', 'Unknown')}")
    logging.info(f"Submit time:                                                      {job_details.get('Submit', 'Unknown')}")
    logging.info(f"Time limit:                                                       {job_details.get('Timelimit', 'Unknown')}")
    logging.info(f"Number of CPU(s) requested:                                       {job_details.get('ReqCPUS', 'Unknown')}")
    logging.info(f"Memory Requested:                                                 {req_mem_gb:.2f} GB")

if __name__ == "__main__":
    main()
