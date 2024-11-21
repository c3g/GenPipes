#!/usr/bin/env python3
import os
import subprocess
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def get_slurm_job_info(job_id):
    try:
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=JobID,JobName,Submit,Start,State,AllocCPUS,ReqMem,Elapsed"],
            capture_output=True, text=True, check=True
        )
        logging.debug(f"sacct output: {result.stdout}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error retrieving job info: {e}")
        return None

def parse_slurm_output(output):
    if not output:
        logging.warning("Empty SLURM output")
        return {}
    lines = output.strip().split('\n')
    if len(lines) < 2:
        logging.warning("Unexpected SLURM output format")
        return {}
    headers = lines[0].split()
    values = lines[1].split()
    job_info = dict(zip(headers, values))
    return job_info

def main():
    logging.debug("Starting prologue script")
    # Determine scheduler and set environment variables accordingly
    if 'PBS_JOBID' in os.environ:
        logging.debug("Detected PBS scheduler")
        # PBS Environment Variables
        job_id = os.getenv('PBS_JOBID', 'Unknown')
        job_state = os.getenv('PBS_JOBSTATE', 'Unknown')
        submit_time = os.getenv('PBS_SUBMIT_TIME', 'Unknown')
        start_time = os.getenv('PBS_START_TIME', 'Unknown')
        num_cpus = os.getenv('PBS_NUM_PPN', 'Unknown')
        mem = os.getenv('PBS_MEM', 'Unknown')
        requeue = os.getenv('PBS_REQUEUE', 'Unknown')
        runtime = os.getenv('PBS_RUNTIME', 'Unknown')
    elif 'SLURM_JOB_ID' in os.environ:
        logging.debug("Detected SLURM scheduler")
        # SLURM Environment Variables
        job_id = os.getenv('SLURM_JOB_ID', 'Unknown')
        job_info_output = get_slurm_job_info(job_id)
        if job_info_output:
            job_info = parse_slurm_output(job_info_output)
            job_state = job_info.get('State', 'Unknown')
            submit_time = job_info.get('Submit', 'Unknown')
            start_time = job_info.get('Start', 'Unknown')
            num_cpus = job_info.get('AllocCPUS', 'Unknown')
            mem = job_info.get('ReqMem', 'Unknown')
            runtime = job_info.get('Elapsed', 'Unknown')
        else:
            logging.warning("Failed to retrieve SLURM job info")
            return
    else:
        logging.warning("Unknown scheduler")
        return

    logging.info(f"Prologue: Setting up environment for job {job_id}")
    logging.info(f"Job State: {job_state}")
    logging.info(f"Submit Time: {submit_time}")
    logging.info(f"Start Time: {start_time}")
    logging.info(f"Number of CPUs: {num_cpus}")
    logging.info(f"Memory: {mem}")
    logging.info(f"Runtime: {runtime}")
    print()

if __name__ == "__main__":
    main()
