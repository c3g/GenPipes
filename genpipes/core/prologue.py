#!/cvmfs/soft.mugqic/CentOS6/software/python/Python-3.12.2/bin/python

import csv
import os
import subprocess
import logging
import time
import sys
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
                    if job_details['Timelimit']:
                        return job_details
                time.sleep(delay)
            else:
                time.sleep(delay)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error retrieving job info: {e}")
            return None
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

def parse_pbs_job_info():
    """
    Parse the job information retrieved from PBS.
    """
    job_details = {}

    requested_resource_limits = sys.argv[5].split(",")
    job_details['JobID'] = sys.argv[1]
    job_details['JobName'] = sys.argv[4]
    job_details['User'] = sys.argv[2]
    job_details['NodeList'] = "Unknown"
    job_details['Priority'] = "Unknown"
    job_details['Submit'] = "Unknown"
    job_details['Timelimit'] = requested_resource_limits[3].split("=")[1]
    job_details['ReqCPUS'] = requested_resource_limits[1].split("=")[2]
    job_details['ReqMem'] = "Unknown"

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
    print("Running GenPipes Prologue")
    job_id = os.getenv('SLURM_JOB_ID')

    if 'SLURM_JOB_ID' in os.environ:
        job_details = get_slurm_job_info(job_id)
        if not job_details:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    # Built-in support for PBS
    elif len(sys.argv) > 5 and sys.argv[5].startswith('neednodes'):
        job_details = parse_pbs_job_info()
        if not job_details:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    else:
        logging.error("Unsupported scheduler")
        return

    # Convert memory to GB
    req_mem_gb = None
    if job_details['ReqMem']:
        req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])

    logging.info(f"GenPipes Prologue for job {job_id}")
    logging.info(f"Job name:                                                         {job_details.get('JobName', 'Unknown')}")
    logging.info(f"User:                                                             {job_details.get('User', 'Unknown')}")
    logging.info(f"Node(s):                                                          {job_details.get('NodeList', 'Unknown')}")
    logging.info(f"Priority:                                                         {job_details.get('Priority', 'Unknown')}")
    logging.info(f"Submit time:                                                      {job_details.get('Submit', 'Unknown')}")
    logging.info(f"Time limit:                                                       {job_details.get('Timelimit', 'Unknown')}")
    logging.info(f"Number of CPU(s) requested:                                       {job_details.get('ReqCPUS', 'Unknown')}")
    if req_mem_gb:
        logging.info(f"Memory Requested:                                                 {req_mem_gb:.2f} GB")
    else:
        logging.info(f"Memory Requested:                                                 Unknown")

if __name__ == "__main__":
    main()
