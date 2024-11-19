#!/usr/bin/env python3
import os
import subprocess

def get_slurm_job_info(job_id):
    try:
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=JobID,JobName,Submit,Start,State,AllocCPUS,ReqMem,Elapsed,RestartCount"],
            capture_output=True, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error retrieving job info: {e}")
        return None

def main():
    # Determine scheduler and set environment variables accordingly
    if 'PBS_JOBID' in os.environ:
        # PBS Environment Variables
        job_id = os.getenv('PBS_JOBID', 'Unknown')
        job_state = os.getenv('PBS_JOBSTATE', 'Unknown')
        submit_time = os.getenv('PBS_SUBMIT_TIME', 'Unknown')
        start_time = os.getenv('PBS_START_TIME', 'Unknown')
        num_cpus = os.getenv('PBS_NUM_PPN', 'Unknown')
        mem = os.getenv('PBS_MEM', 'Unknown')
        requeue = os.getenv('PBS_REQUEUE', 'Unknown')
        restarts = os.getenv('PBS_RESTARTS', 'Unknown')
        runtime = os.getenv('PBS_RUNTIME', 'Unknown')
    elif 'SLURM_JOB_ID' in os.environ:
        # SLURM Environment Variables
        job_id = os.getenv('SLURM_JOB_ID', 'Unknown')
        job_info = get_slurm_job_info(job_id)
        if job_info:
            print("-" * 79)
            print(f"Epilogue: Cleaning up environment for job {job_id}")
            print(job_info)
        return
    else:
        print("Unknown scheduler")
        return

    print("-" * 79)
    print(f"Epilogue: Cleaning up environment for job {job_id}")
    print(f"Job State: {job_state}")
    print(f"Submit Time: {submit_time}")
    print(f"Start Time: {start_time}")
    print(f"Number of CPUs: {num_cpus}")
    print(f"Memory: {mem}")
    print(f"Requeue: {requeue}")
    print(f"Restarts: {restarts}")
    print(f"Runtime: {runtime}")

if __name__ == "__main__":
    main()
