#!/usr/bin/env python3
import os
import subprocess
import time

def get_slurm_job_info(job_id):
    try:
        print(f"Fetching job info for job ID: {job_id}")
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=JobID,JobName,Submit,Start,State,AllocCPUS,ReqMem,Elapsed"],
            capture_output=True, text=True, check=True
        )
        print("Job info fetched successfully")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error retrieving job info: {e}")
        return None

def main():
    if 'SLURM_JOB_ID' in os.environ:
        job_id = os.getenv('SLURM_JOB_ID', 'Unknown')
        print(f"SLURM_JOB_ID: {job_id}")

        # Adding a delay to ensure job status is updated
        time.sleep(5)

        job_info = get_slurm_job_info(job_id)
        if job_info:
            print("-" * 79)
            print(f"Epilogue: Cleaning up environment for job {job_id}")
            print(job_info)
        else:
            print(f"Failed to retrieve job info for job {job_id}")
        return
    else:
        print("Unknown scheduler")
        return

if __name__ == "__main__":
    main()
