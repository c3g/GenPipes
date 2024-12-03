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

# Configure logging
logging.basicConfig(level=logging.INFO, format='EPILOGUE - %(message)s')

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
                ["sacct", "-j", f"{job_id}", "--parsable", "--format=JobID,JobName,User,NodeList,Priority,Submit,Eligible,Timelimit,ReqCPUS,ReqMem,State,Start,End,Elapsed,TotalCPU,AveRSS,MaxRSS,AveDiskRead,MaxDiskRead,AveDiskWrite,MaxDiskWrite"],
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
    logging.error(f"Failed to retrieve complete job info for job {job_id} after {retries} attempts")
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
    # logging.info(f"Job info: {job_info}")
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
    job_details['TotalMem'] = 'Unknown'

    # Check if all necessary fields are populated
    required_fields = ['JobID', 'JobName', 'User', 'NodeList', 'Priority', 'Submit', 'Eligible', 'Timelimit', 'ReqCPUS', 'ReqMem', 'State', 'Start', 'End', 'Elapsed', 'TotalCPU', 'AveRSS', 'MaxRSS', 'AveDiskRead', 'MaxDiskRead', 'AveDiskWrite', 'MaxDiskWrite']
    missing_fields = [field for field in required_fields if field not in job_details]
    if missing_fields:
        logging.warning(f"Missing fields: {', '.join(missing_fields)}")
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
        ["qstat", "-f", f"{job_id}"],
        capture_output=True, text=True, check=True
    )
    job_info = result.stdout
    job_details = parse_pbs_job_info(job_id, job_info)
    if not result.stdout.strip():
        logging.warning("Error retrieving job info. The prologue will be incomplete.")
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

def pbs_exit_code_to_string(exit_code):
    """
    Convert PBS exit code to string.
    Args:
        exit_code (int): Exit code.
    Returns:
        str: Exit code string if recognized, otherwise "Unknown".
    """
    if exit_code == 0:
        return "COMPLETED"
    if 1 <= exit_code <= 127:
        return f"ERROR ({os.strerror(exit_code)})"
    if exit_code >= 128:
        signal_number = exit_code - 128
        return f"TERMINATED ({signal.Signals(signal_number).name})"
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
    logging.info(f"Job info: {job_info}")

    job_details['JobID'] = job_id
    job_name_match = re.search(r"Job_Name\s*=\s*(.+)", job_info)
    job_details['JobName'] = job_name_match.group(1) if job_name_match else "Unknown"
    job_owner_match = re.search(r"Job_Owner\s*=\s*(.+)@", job_info)
    job_details['User'] = job_owner_match.group(1) if job_owner_match else "Unknown"
    exec_host_match = re.search(r"exec_host\s*=\s*(.+)", job_info)
    job_details['NodeList'] = exec_host_match.group(1) if exec_host_match else "Unknown"
    priority_match = re.search(r"Priority\s*=\s*(.+)", job_info)
    job_details['Priority'] = priority_match.group(1) if priority_match else "Unknown"
    job_details['Submit'] = parse_datetime(job_info, "qtime")
    job_details['Eligible'] = parse_datetime(job_info, "etime")
    walltime_match = re.search(r'Resource_List\.walltime\s*=\s*(.+)', job_info)
    job_details['Timelimit'] = walltime_match.group(1) if walltime_match else "Unknown"
    ppn_match = re.search(r'Resource_List\.nodes\s*=\s*\d+:ppn=(\d+)', job_info)
    job_details['ReqCPUS'] = ppn_match.group(1) if ppn_match else "Unknown"
    # job_details['State'] = pbs_exit_code_to_string(int(sys.argv[10]))
    job_details['State'] = "Unknown"
    job_details['Start'] = parse_datetime(job_info, "start_time")
    job_details['End'] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    # walltime_match_used = re.search(r'walltime=([\d:]+)', used_resource)
    # job_details['Elapsed'] = walltime_match_used.group(1) if walltime_match_used else "Unknown"
    job_details['Elapsed'] = parse_datetime(job_info, "resources_used.walltime")
    # cput_match = re.search(r'cput=([\d:]+)', used_resource)
    cput_match = re.search(r'resources_used\.cput=([\d:]+)', job_info)
    job_details['TotalCPU'] = cput_match.group(1) if cput_match else "Unknown"
    # mem_match = re.search(r',mem=([\d]+\w\w)', used_resource)
    mem_match = re.search(r'resources_used\.mem=([\d]+\w\w)', job_info)
    job_details['TotalMem'] = mem_match.group(1) if mem_match else "Unknown"
    job_details['AveRSS'] = "Unknown"
    job_details['MaxRSS'] = "Unknown"
    job_details['AveDiskRead'] = "Unknown"
    job_details['MaxDiskRead'] = "Unknown"
    job_details['AveDiskWrite'] = "Unknown"
    job_details['MaxDiskWrite'] = "Unknown"

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
        str: Time difference in DD:HH:MM:SS format.
    """
    start_dt = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S")
    end_dt = datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S")
    delta = end_dt - start_dt
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{days:02}:{hours:02}:{minutes:02}:{seconds:02}"

def calculate_time_efficency(elapsed, timelimit):
    """
    Calculate the time efficiency of a job.
    Args:
        elapsed (str): Elapsed time.
        timelimit (str): Time limit.
    Returns:
        float: Time efficiency percentage.
    """
    elapsed_seconds = sum(int(x) * 60 ** i for i, x in enumerate(reversed(elapsed.split(':'))))
    timelimit_seconds = sum(int(x) * 60 ** i for i, x in enumerate(reversed(timelimit.split(':'))))
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

def main():
    """
    Main function to run the epilogue script.
    """
    job_id = os.getenv('SLURM_JOB_ID') or os.getenv('PBS_JOBID')

    if 'SLURM_JOB_ID' in os.environ:
        job_details = get_slurm_job_info(job_id)
        if not job_details:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    elif 'PBS_JOBID' in os.environ:
        job_details = get_pbs_job_info(job_id)
        if not job_details:
            logging.error(f"Failed to retrieve job info for job {job_id}")
            return
    else:
        logging.error("Unsupported scheduler")
        return

    # Calculate time spent in queue format DD:HH:MM:SS
    time_in_queue = None
    if job_details['Submit'] != "Unknown":
        time_in_queue = calculate_time_difference(job_details['Submit'], job_details['Start'])
    # Calculate time efficency between walltime and time used
    time_efficency = calculate_time_efficency(job_details['Elapsed'], job_details['Timelimit'])
    # Convert memory to GB
    req_mem_gb = None
    if job_details['ReqMem'] != "Unknown":
        req_mem_gb = convert_memory_to_gb(job_details['ReqMem'])
    ave_mem_gb = None
    if job_details['AveRSS'] != "Unknown":
        ave_mem_gb = convert_memory_to_gb(job_details['AveRSS'])
    max_mem_gb = None
    if job_details['MaxRSS'] != "Unknown":
        max_mem_gb = convert_memory_to_gb(job_details['MaxRSS'])
    total_mem_gb = None
    if job_details['TotalMem'] != "Unknown":
        total_mem_gb = convert_memory_to_gb(job_details['TotalMem'])
    # Calculate percentages for CPU and memory usage
    elapsed = time_str_to_seconds(job_details.get('Elapsed', '00:00:00'))
    total_cpu = time_str_to_seconds(job_details.get('TotalCPU', '00:00:00'))
    cpu_usage_percentage = calculate_percentage(total_cpu, elapsed)
    mem_usage_percentage = None
    if req_mem_gb and ave_mem_gb:
        mem_usage_percentage = calculate_percentage(ave_mem_gb, req_mem_gb)

    logging.info(f"GenPipes Epilogue for job {job_id}")
    logging.info(f"Job Name:                                                         {job_details.get('JobName', 'Unknown')}")
    logging.info(f"User:                                                             {job_details.get('User', 'Unknown')}")
    logging.info(f"Node(s):                                                          {job_details.get('NodeList', 'Unknown')}")
    logging.info(f"Priority:                                                         {job_details.get('Priority', 'Unknown')}")
    logging.info(f"Status:                                                           {job_details.get('State', 'Unknown')}")
    logging.info(f"Submit Time:                                                      {job_details.get('Submit', 'Unknown')}")
    logging.info(f"Eligible Time:                                                    {job_details.get('Eligible', 'Unknown')}")
    logging.info(f"Start Time:                                                       {job_details.get('Start', 'Unknown')}")
    if time_in_queue:
        logging.info(f"Time Spent in Queue (DD:HH:MM:SS):                                {time_in_queue}")
    else:
        logging.info(f"Time Spent in Queue:                                              Unknown")
    logging.info(f"End Time:                                                         {job_details.get('End', 'Unknown')}")
    logging.info(f"Total Wall-clock Time:                                            {job_details.get('Elapsed', 'Unknown')}")
    logging.info(f"Time Limit:                                                       {job_details.get('Timelimit', 'Unknown')}")
    logging.info(f"Time Efficiency (% of Wall-clock Time to Time Limit):             {time_efficency:.1f}%")
    logging.info(f"Number of CPU(s) Requested:                                       {job_details.get('ReqCPUS', 'Unknown')}")
    logging.info(f"Total CPU Time:                                                   {job_details.get('TotalCPU', 'Unknown')}")
    logging.info(f"CPU Efficiency (% of CPU Time to Wall-clock Time):                {cpu_usage_percentage:.1f}%")
    logging.info(f"Memory Requested:                                                 {req_mem_gb:.2f} GB")
    if total_mem_gb is not None:
        logging.info(f"Total Memory:                                                     {total_mem_gb:.2f} GB")
    else:
        logging.info(f"Total Memory:                                                     Unknown")
    if ave_mem_gb is not None:
        logging.info(f"Average Memory Usage:                                             {ave_mem_gb:.2f} GB")
    else:
        logging.info(f"Average Memory Usage:                                             Unknown")
    if max_mem_gb is not None:
        logging.info(f"Maximum Memory Usage:                                             {max_mem_gb:.2f} GB")
    else:
        logging.info(f"Maximum Memory Usage:                                             Unknown")
    if mem_usage_percentage:
        logging.info(f"Memory Efficiency (% of Memory Requested to Average Memory Used): {mem_usage_percentage:.1f}%")
    else:
        logging.info(f"Memory Efficiency (% of Memory Requested to Average Memory Used): Unknown")
    logging.info(f"Average Disk Read:                                                {job_details.get('AveDiskRead', 'Unknown')}")
    logging.info(f"Maximum Disk Read:                                                {job_details.get('MaxDiskRead', 'Unknown')}")
    logging.info(f"Average Disk Write:                                               {job_details.get('AveDiskWrite', 'Unknown')}")
    logging.info(f"Maximum Disk Write:                                               {job_details.get('MaxDiskWrite', 'Unknown')}")

if __name__ == "__main__":
    main()
