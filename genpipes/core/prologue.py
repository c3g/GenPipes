#!/cvmfs/soft.mugqic/CentOS6/software/python/Python-3.12.2/bin/python

import csv
import os
import subprocess
import logging
import re
import time
from io import StringIO
from datetime import datetime

# Create a custom logger
logger = logging.getLogger('logging')
logger.setLevel(logging.INFO)

# Create a console handler and set the format
console_handler = logging.StreamHandler()
formatter = logging.Formatter('PROLOGUE - %(message)s')
console_handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(console_handler)

def log_separator():
    """
    Use the handler's stream to print the separator without the prefix.
    """
    console_handler.stream.write('\n' + '-' * 90 + '\n')
    console_handler.flush()

def get_slurm_job_info(job_id, retries=20, delay=10):
    """
    Retrieve job information from SLURM using sacct command with retries.
    Args:
        job_id (str): Job ID.
        retries (int): Number of retries.
        delay (int): Delay in seconds between retries.
    Returns:
        dict: Job details if found, otherwise None.
    """
    for _ in range(retries):
        try:
            result = subprocess.run(
                ["sacct", "-j", job_id, "--parsable", "--format=JobID,JobName,User,NodeList,Priority,Submit,Timelimit,ReqCPUS,ReqMem"],
                capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                job_info = result.stdout
                job_details = parse_slurm_job_info(job_info, job_id)
                if job_details:
                    if job_details['Timelimit'] and job_details['NodeList'] != "None assigned":
                        return job_details
                time.sleep(delay)
            else:
                time.sleep(delay)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error retrieving job info: {e}")
            return None
    logger.error(f"Failed to retrieve complete job info for job {job_id} after {retries} attempts")
    return None

def parse_slurm_job_info(job_info, job_id):
    """
    Parse the job information retrieved from SLURM.
    Args:
        job_info (str): Job information.
        job_id (str): Job ID.
    Returns:
        dict: Job details if found, otherwise None.
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
        # logger.warning(f"Missing fields: {', '.join(missing_fields)}")
        return None
    return job_details

def get_pbs_job_info(job_id):
    """
    Retrieve job information from PBS using qstat command.
    Args:
        job_id (str): Job ID.
    Returns:
        dict: Job details if found, otherwise None.
    """
    result = subprocess.run(
        ["qstat", "-f", job_id],
        capture_output=True, text=True, check=True
    )
    job_info = result.stdout
    job_details = parse_pbs_job_info(job_id, job_info)
    if not result.stdout.strip():
        logger.warning("Error retrieving job info. The prologue will be incomplete.")
    return job_details

def parse_datetime(job_details, field_name):
    """
    Parse datetime field from job details.
    Args:
        job_details (str): Job details string.
        field_name (str): Field name to search for.
    Returns:
        str: Datetime string in ISO format if found, otherwise "Unknown".
    """
    pattern = rf"{field_name} = (.+)"
    match = re.search(pattern, job_details)
    if match:
        return datetime.strptime(match.group(1), "%a %b %d %H:%M:%S %Y").strftime("%Y-%m-%dT%H:%M:%S")
    return "Unknown"

def parse_pbs_job_info(job_id, job_info):
    """
    Parse the job information retrieved from PBS.
    Args:
        job_id (str): Job ID.
        job_info (str): Job information.
    Returns:
        dict: Job details if found, otherwise None.
    """
    job_details = {}

    job_details['JobID'] = job_id
    job_name_match = re.search(r"Job_Name\s*=\s*(.+)", job_info)
    job_details['JobName'] = job_name_match.group(1) if job_name_match else "Unknown"
    job_owner_match = re.search(r"Job_Owner\s*=\s*(.+)@", job_info)
    job_details['User'] = job_owner_match.group(1) if job_owner_match else "Unknown"
    exec_host_match = re.search(r"exec_host\s*=\s*(.+)/", job_info)
    job_details['NodeList'] = exec_host_match.group(1) if exec_host_match else "Unknown"
    priority_match = re.search(r"Priority\s*=\s*(.+)", job_info)
    job_details['Priority'] = priority_match.group(1) if priority_match else "Unknown"
    job_details['Submit'] = parse_datetime(job_info, "qtime")
    walltime_match = re.search(r'Resource_List\.walltime\s*=\s*(.+)', job_info)
    job_details['Timelimit'] = walltime_match.group(1) if walltime_match else "Unknown"
    ppn_match = re.search(r'Resource_List\.nodes\s*=\s*\d+:ppn=(\d+)', job_info)
    job_details['ReqCPUS'] = ppn_match.group(1) if ppn_match else "Unknown"

    queue_match = re.search(r"queue\s*=\s*(.+)", job_info)
    if queue_match.group(1) == "lm":
        job_details['ReqMem'] = f"{int(job_details['ReqCPUS']) * 15}G"
    else:
        job_details['ReqMem'] = f"{int(job_details['ReqCPUS']) * 5}G"

    return job_details

def convert_memory_to_gb(memory_str):
    """
    Convert memory string to GB.
    Args:
        memory_str (str): Memory string.
    Returns:
        float: Memory in GB.
    """
    if memory_str.endswith('K'):
        return float(memory_str[:-1]) / (1024 ** 2)
    if memory_str.endswith('kb'):
        return float(memory_str[:-2]) / (1024 ** 2)
    if memory_str.endswith('M'):
        return float(memory_str[:-1]) / 1024
    if memory_str.endswith('mb'):
        return float(memory_str[:-2]) / 1024
    if memory_str.endswith('G'):
        return float(memory_str[:-1])
    if memory_str.endswith('gb'):
        return float(memory_str[:-2])
    if memory_str.endswith('T'):
        return float(memory_str[:-1]) * 1024
    if memory_str.endswith('tb'):
        return float(memory_str[:-2]) * 1024
    return float(memory_str)

def custom_get(dictionary, key, default="Unknown"):
    """
    Custom get method for dictionaries.
    Args:
        dictionary (dict): Dictionary.
        key (str): Key to search for.
        default (str): Default value if key is not found.
    Returns:
        str: Value if found, otherwise default value.
    """
    value = dictionary.get(key, default)
    return value if value is not None else default

def main():
    """
    Main function to run the prologue script.
    """
    job_id = os.getenv('SLURM_JOB_ID') or os.getenv('PBS_JOBID')

    if 'SLURM_JOB_ID' in os.environ:
        job_details = get_slurm_job_info(job_id)
        if not job_details:
            logger.error(f"Failed to retrieve job info for job {job_id}")
            return
    elif 'PBS_JOBID' in os.environ:
        job_details = get_pbs_job_info(job_id)
        if not job_details:
            logger.error(f"Failed to retrieve job info for job {job_id}")
            return
    else:
        logger.error("Unsupported scheduler")
        return

    # Convert memory to GB
    req_mem_gb = None
    if job_details['ReqMem'] != "Unknown":
        req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])

    logger.info(f"GenPipes Prologue for job {job_id}")
    logger.info(f"Job Name:                                                         {custom_get(job_details, 'JobName')}")
    logger.info(f"User:                                                             {custom_get(job_details, 'User')}")
    logger.info(f"Node(s):                                                          {custom_get(job_details, 'NodeList')}")
    logger.info(f"Priority:                                                         {custom_get(job_details, 'Priority')}")
    logger.info(f"Submit Time:                                                      {custom_get(job_details, 'Submit')}")
    logger.info(f"Time Limit:                                                       {custom_get(job_details, 'Timelimit')}")
    logger.info(f"Number of CPU(s) Requested:                                       {custom_get(job_details, 'ReqCPUS')}")
    logger.info(f"Memory Requested:                                                 {f'{req_mem_gb:.2f} GB' if req_mem_gb is not None else 'Unknown'}")
    log_separator()

if __name__ == "__main__":
    main()
