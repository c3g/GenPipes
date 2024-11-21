#!/usr/bin/env python3
import os
import subprocess
import time
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_slurm_job_info(job_id):
    try:
        logging.info(f"Fetching job info for job ID: {job_id}")
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=JobID,JobName,Submit,Start,State,AllocCPUS,ReqMem,Elapsed", "--noheader"],
            capture_output=True, text=True, check=True
        )
        logging.info("Job info fetched successfully")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error retrieving job info: {e}")
        return None

def parse_job_info(job_info):
    lines = job_info.strip().split('\n')
    job_details = {}
    if lines:
        # Assuming the first line contains the job details
        details = lines[0].split()
        job_details['JobID'] = details[0]
        job_details['JobName'] = details[1]
        job_details['Submit'] = details[2]
        job_details['Start'] = details[3]
        job_details['State'] = details[4]
        job_details['AllocCPUS'] = details[5]
        job_details['ReqMem'] = details[6]
        job_details['Elapsed'] = details[7]
    return job_details

def main():
    if 'SLURM_JOB_ID' in os.environ:
        job_id = os.getenv('SLURM_JOB_ID', 'Unknown')
        logging.info(f"SLURM_JOB_ID: {job_id}")

        # Adding a delay to ensure job status is updated
        time.sleep(5)

        job_info = get_slurm_job_info(job_id)
        if job_info:
            job_details = parse_job_info(job_info)
            logging.info("-" * 79)
            logging.info(f"Epilogue: Cleaning up environment for job {job_id}")
            logging.info(f"Job State: {job_details.get('State', 'Unknown')}")
            logging.info(f"Submit Time: {job_details.get('Submit', 'Unknown')}")
            logging.info(f"Start Time: {job_details.get('Start', 'Unknown')}")
            logging.info(f"Number of CPUs: {job_details.get('AllocCPUS', 'Unknown')}")
            logging.info(f"Memory: {job_details.get('ReqMem', 'Unknown')}")
            logging.info(f"Runtime: {job_details.get('Elapsed', 'Unknown')}")
        else:
            logging.error(f"Failed to retrieve job info for job {job_id}")
        return
    else:
        logging.error("Unknown scheduler")
        return

if __name__ == "__main__":
    main()
