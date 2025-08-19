import yaml
import copy
from simple_slurm import Slurm
from job_class import slurm_job

def read_steps(steps_path: str = './steps.yaml'):
    with open(steps_path, "r") as f:
        steps = yaml.safe_load(f)
    return steps

def create_job_steps(steps, steps_def, job_org: slurm_job):
    step_jobs = dict()
    for step in steps:
        if step not in steps_def:
            raise KeyError("Undefined step")
        step_jobs[step] = copy.deepcopy(job_org)
        try:
            step_jobs[step].options.update(steps_def[step]['slurm'])
            step_jobs[step].slurm = Slurm(**step_jobs[step].options)
        except:
            pass
        step_jobs[step].run_commands = steps_def[step]['cmd']
        step_jobs[step].add_cmds()
    return step_jobs

if __name__ == "__main__":
    print(read_steps('/scratch/kwlg019/TESTs/slurm-pred/steps.yaml')['test_1']['cmd'])
