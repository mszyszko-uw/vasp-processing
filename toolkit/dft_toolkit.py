import sys
from simple_slurm import Slurm
import subprocess
import os
import copy
from datetime import datetime
from slurm_waiting_prediction import SLURM_make_times_report, get_available_resources
from parsers import parse_part_values, read_cli_params, print_avalable_steps
from read_steps import read_steps, create_job_steps

from set_parallel import set_parallelization
from run_local import run_local

from job_class import slurm_job

def main():
    params = read_cli_params()
    print(params.action[0])
    job = slurm_job()
    steps = []
    steps_jobs = []

    if params.config:
        job.read_config(params.config)
    if params.steps:
        steps=parse_part_values(params.steps)
        steps_def = read_steps()
        steps_jobs = create_job_steps(steps, steps_def, job)
    if params.path:
        if params.steps:
            for step in steps:
                steps_jobs[step].set_path(params.path)
        else:
            job.set_path(params.path)
    if params.optimize_cpus:
	#RUN dry_run
        if params.steps:
            for step in steps:
                run_local(steps_def, step, vasp_command="mpiexec vasp_std" )
                cpus= set_parallelization("./", params.optimize_cpus)
                steps_jobs[step].add_option('ntasks', cpus)
        else:
            cpus= set_parallelization(params.optimize_cpus)
            job.add_option('ntasks', cpus)
    if params.array:
        if params.steps:
            for step in steps:
                steps_jobs[step].add_array_from_path()
#                steps_jobs[step].submit()
        else:
            job.add_array_from_path()
#            job.submit()

    if params.dependency_step:
        if params.steps:
            for step in steps:
                steps_jobs[step].dependency = params.dependency_step
        else:
            job.dependency = params.dependency_step

    action = params.action[0]
    if action == "freenodes":
        get_available_resources(printing=True)
    if action == "wating":
        SLURM_make_times_report()
    if action == "print":
        if params.steps:
            for step in steps:
                steps_jobs[step].print()
        else:
            job.print()
    if action == "create":
        if params.steps:
            for step in steps:
                sub_name="./JOBS/STEP0"+str(step)
                steps_jobs[step].write(sub_name)
        else:
            job.write()
    if action == 'print_avalable_steps':
        print_avalable_steps()

    if action == 'submit':
        if params.steps:
#            job.submit()
            job_id = 0
            for step in steps:
                if job_id:
                    print(job_id)
                    job_id = steps_jobs[step].submit_dependence(job_id)
                else:
                    job_id = steps_jobs[step].submit()
        else:
            job.submit()
    if action == 'checkqueue':
        subprocess.run(["squeue", "--me",])
    if params.id:
        if action == 'canceljob':
            subprocess.run(["scancel", params.id])
        if action == 'jobinfo':
            subprocess.run(["scontrol", "show", "job", params.id])
#    job.predict_time()


if __name__ == "__main__":
    main()
