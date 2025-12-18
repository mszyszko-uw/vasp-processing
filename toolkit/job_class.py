import sys
import argparse
from simple_slurm import Slurm
import datetime
import json
import yaml
import subprocess
import os
import time
import copy
from datetime import datetime
from slurm_waiting_prediction import SLURM_make_times_report, get_available_resources
from parsers import parse_part_values, read_cli_params

class slurm_job:
    def __init__(self,
                path: str = "./",
                job_name: str = "my_job",
                output: str = "output-%j.out",
                error: str = "output-%j.err",
                nodes: int = 1,
                ntasks: int = 1,
                cpus_per_task: int = 1,
                mem: str = "1G",
                time: str = "00:30:00",
                partition: str = "general",
                account: str = None,
                additional_options: dict = None,
                commands: list = None,
                module_name: str = "VASP",
                run_commands: str = "vasp_parallel",
                dependency: str = '',
                input_file: str = ''):
        self.options = {
                "job_name": job_name,
                "output": output,
#                "ntasks_per_node": 72,
#                "nodes": nodes,
#                "ntasks": ntasks,
#                "cpus_per_task": cpus_per_task,
#                "mem": mem,
#                "time": time,
                "partition": partition,}
        if error:
            self.options["error"] = error
        if account:
            self.options["account"] = account
        if additional_options:
            self.options.update(additional_options)
        self.slurm = Slurm(**self.options)
        self.path = path
        self.module_name = module_name
        self.env_cmd = ''
        self.run_commands = run_commands
        self.env_array = []
        self.input_file = input_file
        self.cmds = ["module load " + self.module_name, self.run_commands]
        self.add_cmds()
        self.dependency = dependency

    def add_cmds(self):
        self.cmds = self.env_array + ["module load " + self.module_name, self.env_cmd, f'INPUT="{self.input_file}"', self.run_commands]
        for cmd in self.cmds:
            self.slurm.add_cmd(cmd)
        self.slurm.set_shell("/bin/bash -l")
         
    def print(self):
        print(self.slurm.__str__())

    def write(self, sub_name: str = './'):
        try:
            os.makedirs(sub_name, exist_ok=True)
        except:
            pass
        with open(sub_name+'/run.sh', 'w') as f:
            f.write(self.slurm.__str__())

    def submit(self):
        return self.slurm.sbatch()
    
    def submit_dependence(self, dependence:int):
        #dependency="afterok:12345"
        self.add_option("dependency", "afterok:"+str(dependence))
        return self.slurm.sbatch()

    def read_config(self, config_path: str):
        if config_path.endswith(".json"):
            with open(config_path, "r") as f:
                config = json.load(f)
        elif config_path.endswith((".yaml", ".yml")):
            with open(config_path, "r") as f:
                config = yaml.safe_load(f)
        else:
            raise ValueError("Unsupported file format. Use .json, .yaml, or .yml")
        self.options.update(config['slurm'])
        try:
            self.module_name = config['script']['module']
        except:
            pass
        try:
            activation_cmd = ''
            env_type = config['env']['type']
            env_path = config['env']['path']
            if env_type == 'venv':
                activation_cmd = '. ' + env_path + r'/bin/activate'
            if env_type == 'conda':
                activation_cmd = '. activation ' + env_path
            self.env_cmd = activation_cmd
        except:
            pass
        self.slurm = Slurm(**self.options)
        self.slurm.set_shell("/bin/bash -l")
        self.add_cmds()

    def add_option(self, key, value):
        self.options[key] = value
        self.slurm = Slurm(**self.options)
        self.slurm.set_shell("/bin/bash -l")
        self.add_cmds()

    def set_path(self, path:str):
        self.path = path

    def add_array_from_path(self):
        incar_dirs = []
        for dirpath, dirnames, filenames in os.walk(self.path):
            if 'INCAR' in filenames:
                incar_dirs.append(dirpath)
        array_job_len = len(incar_dirs)-1
        results = 'jobs=('+' '.join(f'"{s}"' for s in incar_dirs)+')'
        self.options["array"] = "0-"+str(array_job_len)
        self.options["output"] = "output_%A_%a.out"
        self.options["error"] = "output_%A_%a.err"
        self.env_array = [results, "cd ${jobs[$SLURM_ARRAY_TASK_ID]}"]
        self.slurm = Slurm(**self.options)
        self.slurm.set_shell("/bin/bash -l")
        self.add_cmds()

    def submit_with_hold(self):
        options_org = self.options
        self.add_option("--hold"," ")
        job_id = self.slurm.sbatch(convert=False)
        self.options = options_org
        return job_id

    def get_start_time(aelf, jobid, timeout=20):
        for _ in range(timeout):
            result = subprocess.run(["squeue", "--start", "-j", str(jobid), "--format=%A %S"],
                                    capture_output=True, text=True)
            lines = result.stdout.strip().splitlines()
            if len(lines) >= 2:
                parts = lines[1].split()
                if len(parts) >= 2:
                    start_time = datetime.strptime(parts[1], "%Y-%m-%dT%H:%M:%S")
                    now = datetime.now()
                    return int((start_time - now).total_seconds())/3600
            else: return 0
            time.sleep(1)
        raise TimeoutError("Failed to read estimated launch time.")
