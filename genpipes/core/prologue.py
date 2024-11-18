#!/usr/bin/env python3
import os

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
        job_state = os.getenv('SLURM_JOB_STATE', 'Unknown')
        submit_time = os.getenv('SLURM_SUBMIT_TIME', 'Unknown')
        start_time = os.getenv('SLURM_START_TIME', 'Unknown')
        num_cpus = os.getenv('SLURM_CPUS_ON_NODE', 'Unknown')
        mem = os.getenv('SLURM_MEM_PER_NODE', 'Unknown')
        requeue = os.getenv('SLURM_REQUEUE', 'Unknown')
        restarts = os.getenv('SLURM_RESTART_COUNT', 'Unknown')
        runtime = os.getenv('SLURM_RUNTIME', 'Unknown')
    else:
        print("Unknown scheduler")
        return

    print(f"Prologue: Setting up environment for job {job_id}")
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
