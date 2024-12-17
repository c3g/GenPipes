#!/cvmfs/soft.mugqic/CentOS6/software/python/Python-3.12.2/bin/python

import csv
import os
import subprocess
import logging
import re
import time
import signal
import sys
from io import StringIO
from datetime import datetime

# Create a custom logger
logger = logging.getLogger('logging')
logger.setLevel(logging.INFO)

# Create a console handler and set the format
console_handler = logging.StreamHandler()
formatter = logging.Formatter('EPILOGUE - %(message)s')
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
                ["sacct", "-j", job_id, "--parsable", "--format=JobID,JobName,User,NodeList,Priority,Submit,Eligible,Timelimit,ReqCPUS,ReqMem,State,Start,End,Elapsed,TotalCPU,AveRSS,MaxRSS,AveDiskRead,MaxDiskRead,AveDiskWrite,MaxDiskWrite"],
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
    # logger.info(f"Job info: {job_info}")
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

        # Extracting from the third line (line starting with ${JobID}.0)
        elif row['JobID'] == f"{job_id}.0":
            if 'CANCELLED by ' in row['State']:
                job_details['State'] = 'OUT_OF_MEMORY'
            else:
                job_details['State'] = row['State']
            job_details['Start'] = row['Start']
            job_details['End'] = row['End']
            job_details['Elapsed'] = row['Elapsed']
            job_details['TotalCPU'] = row['TotalCPU']
            job_details['AveRSS'] = row['AveRSS']
            job_details['MaxRSS'] = row['MaxRSS']
            job_details['AveDiskRead'] = row['AveDiskRead']
            job_details['MaxDiskRead'] = row['MaxDiskRead']
            job_details['AveDiskWrite'] = row['AveDiskWrite']
            job_details['MaxDiskWrite'] = row['MaxDiskWrite']
            break

    # Check if all necessary fields are populated
    required_fields = ['JobID', 'JobName', 'User', 'NodeList', 'Priority', 'Submit', 'Eligible', 'Timelimit', 'ReqCPUS', 'ReqMem', 'State', 'Start', 'End', 'Elapsed', 'TotalCPU', 'AveRSS', 'MaxRSS', 'AveDiskRead', 'MaxDiskRead', 'AveDiskWrite', 'MaxDiskWrite']
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
        ["/opt/torque/x86_64/bin/qstat", "-f", job_id],
        capture_output=True, text=True, check=True
    )
    job_info = result.stdout
    job_details = parse_pbs_job_info(job_id, job_info)
    if not result.stdout.strip():
        logger.warning("Error retrieving job info. The epilogue will be incomplete.")
    return job_details

def parse_datetime(job_details, field_name):
    """
    Parse datetime field from job details.
    Args:
        job_details (str): Job details string.
        field_name (str): Field name to search for.
    Returns:
        str: Datetime string in ISO format if found, otherwise None.
    """
    pattern = rf"{field_name} = (.+)"
    match = re.search(pattern, job_details)
    if match:
        return datetime.strptime(match.group(1), "%a %b %d %H:%M:%S %Y").strftime("%Y-%m-%dT%H:%M:%S")
    return None

def pbs_exit_code_to_string(exit_code):
    """
    Convert PBS exit code to string.
    Args:
        exit_code (int): Exit code.
    Returns:
        str: Exit code string if recognized, otherwise None.
    """
    if exit_code == 0:
        return "COMPLETED"
    if 1 <= exit_code <= 127:
        return f"ERROR ({os.strerror(exit_code)})"
    if exit_code >= 128:
        signal_number = exit_code - 256
        if signal_number in signal.Signals.__members__.values():
            signal_name = signal.Signals(signal_number).name
            return f"TERMINATED ({signal_name})"
        return f"TERMINATED (Unknown signal {signal_number})"
    return None

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

    used_resource = sys.argv[7]

    job_details['JobID'] = job_id
    job_name_match = re.search(r"Job_Name\s*=\s*(.+)", job_info)
    job_details['JobName'] = job_name_match.group(1) if job_name_match else None
    job_owner_match = re.search(r"Job_Owner\s*=\s*(.+)@", job_info)
    job_details['User'] = job_owner_match.group(1) if job_owner_match else None
    exec_host_match = re.search(r"exec_host\s*=\s*(.+)/", job_info)
    job_details['NodeList'] = exec_host_match.group(1) if exec_host_match else None
    priority_match = re.search(r"Priority\s*=\s*(.+)", job_info)
    job_details['Priority'] = priority_match.group(1) if priority_match else None
    job_details['Submit'] = parse_datetime(job_info, "qtime")
    job_details['Eligible'] = parse_datetime(job_info, "etime")
    walltime_match = re.search(r'Resource_List\.walltime\s*=\s*(.+)', job_info)
    job_details['Timelimit'] = walltime_match.group(1) if walltime_match else None
    ppn_match = re.search(r'Resource_List\.nodes\s*=\s*\d+:ppn=(\d+)', job_info)
    job_details['ReqCPUS'] = ppn_match.group(1) if ppn_match else None
    job_details['State'] = pbs_exit_code_to_string(int(sys.argv[10]))
    job_details['Start'] = parse_datetime(job_info, "start_time")
    job_details['End'] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    walltime_match_used = re.search(r'walltime=([\d:]+)', used_resource)
    job_details['Elapsed'] = walltime_match_used.group(1) if walltime_match_used else None
    cput_match = re.search(r'cput=([\d:]+)', used_resource)
    job_details['TotalCPU'] = cput_match.group(1) if cput_match else None
    mem_match = re.search(r',mem=([\d]+\w\w)', used_resource)
    job_details['AveRSS'] = None
    job_details['MaxRSS'] = mem_match.group(1) if mem_match else None
    job_details['AveDiskRead'] = None
    job_details['MaxDiskRead'] = None
    job_details['AveDiskWrite'] = None
    job_details['MaxDiskWrite'] = None

    queue_match = re.search(r"queue\s*=\s*(.+)", job_info)
    if queue_match.group(1) == "lm":
        job_details['ReqMem'] = f"{int(job_details['ReqCPUS']) * 15}G"
    else:
        job_details['ReqMem'] = f"{int(job_details['ReqCPUS']) * 5}G"

    return job_details

def time_str_to_seconds(time_str):
    """
    Convert time string in HH:MM:SS format to total seconds.
    Args:
        time_str (str): Time string.
    Returns:
        int: Total seconds.
    """
    if not time_str:
        return 0
    parts = time_str.split(':')
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h = '0'
        m, s = parts
    elif len(parts) == 1:
        h = '0'
        m = '0'
        s = parts[0]
    else:
        raise ValueError(f"Unexpected time format: {time_str}")

    # Convert to float to handle decimal seconds
    h = int(h)
    m = int(m)
    s = float(s)

    total_seconds = h * 3600 + m * 60 + s
    return round(total_seconds)

def calculate_time_difference(start_time, end_time):
    """
    Calculate the time difference between two datetime strings.
    Args:
        start_time (str): Start time in ISO format.
        end_time (str): End time in ISO format.
    Returns:
        str: Time difference in D-HH:MM:SS format.
    """
    start_dt = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S")
    end_dt = datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S")
    delta = end_dt - start_dt
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    if days > 0:
        return f"{days}-{hours:02}:{minutes:02}:{seconds:02}"
    return f"{hours:02}:{minutes:02}:{seconds:02}"

def calculate_time_efficiency(elapsed, timelimit):
    """
    Calculate the time efficiency of a job.
    Args:
        elapsed (str): Elapsed time.
        timelimit (str): Time limit.
    Returns:
        float: Time efficiency percentage.
    """
    def parse_time(time_str):
        if '-' in time_str:
            # Format is D-HH:MM:SS
            days, time_part = time_str.split('-')
            days = int(days)
            hours, minutes, seconds = map(int, time_part.split(':'))
        else:
            # Format is HH:MM:SS
            days = 0
            hours, minutes, seconds = map(int, time_str.split(':'))

        return days * 86400 + hours * 3600 + minutes * 60 + seconds

    try:
        elapsed_seconds = parse_time(elapsed)
        timelimit_seconds = parse_time(timelimit)
    except ValueError as e:
        logger.error(f"Error parsing time: {e}")
        return 0

    return 100 * elapsed_seconds / timelimit_seconds if timelimit_seconds > 0 else 0

def calculate_percentage(used, requested):
    """
    Calculate the percentage of used resources to requested resources.
    Args:
        used (float): Used resources.
        requested (float): Requested resources.
    Returns:
        float: Percentage.
    """
    if used is None or requested is None:
        return None
    return (used / requested * 100) if requested > 0 else 0

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
    Main function to run the epilogue script.
    """
    log_separator()
    if 'SLURM_JOB_ID' in os.environ:
        job_id = os.getenv('SLURM_JOB_ID')
        job_details = get_slurm_job_info(job_id)
        if not job_details:
            logger.error(f"Failed to retrieve job info for job {job_id}")
            return
    # Built-in support for PBS
    elif len(sys.argv) > 8:
        job_id = sys.argv[1]
        job_details = get_pbs_job_info(job_id)
        if not job_details:
            logger.error(f"Failed to retrieve job info for job {job_id}")
            return
    else:
        logger.error("Unsupported scheduler")
        return

    # Calculate time spent in queue format HH:MM:SS
    time_in_queue = None
    if job_details['Submit'] is not None:
        time_in_queue = calculate_time_difference(job_details['Eligible'], job_details['Start'])
    # Calculate time efficency between walltime and time used
    time_efficency = calculate_time_efficiency(job_details['Elapsed'], job_details['Timelimit'])
    # Convert memory to GB
    req_mem_gb = None
    if job_details['ReqMem'] is not None:
        req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])
    ave_mem_gb = None
    if job_details['AveRSS'] is not None:
        ave_mem_gb = convert_memory_to_gb(job_details['AveRSS'])
    max_mem_gb = None
    if job_details['MaxRSS'] is not None:
        max_mem_gb = convert_memory_to_gb(job_details['MaxRSS'])
    # Calculate percentages for CPU and memory usage
    total_cpu = None
    if job_details['TotalCPU'] is not None:
        total_cpu = time_str_to_seconds(job_details.get('TotalCPU'))
    elapsed = None
    if job_details['Elapsed'] is not None:
        elapsed = time_str_to_seconds(job_details.get('Elapsed'))
    cpu_usage_percentage = calculate_percentage(total_cpu, elapsed)
    mem_usage_percentage = calculate_percentage(max_mem_gb, req_mem_gb)

    logger.info(f"GenPipes Epilogue for job {job_id}")
    logger.info(f"Job Name:                                                         {custom_get(job_details, 'JobName')}")
    logger.info(f"User:                                                             {custom_get(job_details, 'User')}")
    logger.info(f"Node(s):                                                          {custom_get(job_details, 'NodeList')}")
    logger.info(f"Priority:                                                         {custom_get(job_details, 'Priority')}")
    logger.info(f"Status:                                                           {custom_get(job_details, 'State')}")
    logger.info(f"Submit Time:                                                      {custom_get(job_details, 'Submit')}")
    logger.info(f"Eligible Time:                                                    {custom_get(job_details, 'Eligible')}")
    logger.info(f"Start Time:                                                       {custom_get(job_details, 'Start')}")
    logger.info(f"Time Spent in Queue:                                              {time_in_queue if time_in_queue is not None else 'Unknown'}")
    logger.info(f"End Time:                                                         {custom_get(job_details, 'End')}")
    logger.info(f"Total Wall-clock Time:                                            {custom_get(job_details, 'Elapsed')}")
    logger.info(f"Time Limit:                                                       {custom_get(job_details, 'Timelimit')}")
    logger.info(f"Time Efficiency (% of Wall-clock Time to Time Limit):             {time_efficency:.1f}%")
    logger.info(f"Number of CPU(s) Requested:                                       {custom_get(job_details, 'ReqCPUS')}")
    logger.info(f"Total CPU Time:                                                   {custom_get(job_details, 'TotalCPU')}")
    logger.info(f"CPU Efficiency (% of CPU Time to Wall-clock Time):                {f'{cpu_usage_percentage:.1f}%' if cpu_usage_percentage is not None else 'Unknown'}")
    logger.info(f"Memory Requested:                                                 {f'{req_mem_gb:.2f} GB' if req_mem_gb is not None else 'Unknown'}")
    logger.info(f"Average Memory Usage:                                             {f'{ave_mem_gb:.2f} GB' if ave_mem_gb is not None else 'Unknown'}")
    logger.info(f"Maximum Memory Usage:                                             {f'{max_mem_gb:.2f} GB' if max_mem_gb is not None else 'Unknown'}")
    logger.info(f"Memory Efficiency (% of Memory Requested to Maximum Memory Used): {f'{mem_usage_percentage:.1f}%' if mem_usage_percentage is not None else 'Unknown'}")
    logger.info(f"Average Disk Read:                                                {custom_get(job_details, 'AveDiskRead')}")
    logger.info(f"Maximum Disk Read:                                                {custom_get(job_details, 'MaxDiskRead')}")
    logger.info(f"Average Disk Write:                                               {custom_get(job_details, 'AveDiskWrite')}")
    logger.info(f"Maximum Disk Write:                                               {custom_get(job_details, 'MaxDiskWrite')}")

if __name__ == "__main__":
    main()
