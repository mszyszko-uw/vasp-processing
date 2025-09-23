import subprocess
import re
from simple_slurm import Slurm
import time
from datetime import datetime

def get_available_resources(printing=False):
    result = subprocess.run(["scontrol", "show", "node"], capture_output=True, text=True)
    raw_nodes = result.stdout.strip().split("\n\n")

    resources = []

    for node_info in raw_nodes:
        node_name = re.search(r'NodeName=(\S+)', node_info)
        cpu_total = re.search(r'CPUTot=(\d+)', node_info)
        cpu_alloc = re.search(r'CPUAlloc=(\d+)', node_info)
        mem_total = re.search(r'RealMemory=(\d+)', node_info)
        mem_alloc = re.search(r'AllocMem=(\d+)', node_info)

        if node_name and cpu_total and cpu_alloc and mem_total and mem_alloc:
            name = node_name.group(1)
            free_cpu = int(cpu_total.group(1)) - int(cpu_alloc.group(1))
            free_mem = int(mem_total.group(1)) - int(mem_alloc.group(1))  # w MB

            resources.append({
                'node': name,
                'free_cpu': free_cpu,
                'free_mem_MB': free_mem
            })
    if printing:
        for res in resources:
            print(f"{res['node']}: {res['free_cpu']} CPU free, {res['free_mem_MB']} MB RAM free")

    return resources


def SLURM_get_start_time(jobid, timeout=20):
    for _ in range(timeout):
        result = subprocess.run(["squeue", "--start", "-j", str(jobid), "--format=%A %S"],
                                capture_output=True, text=True)
        lines = result.stdout.strip().splitlines()
        if len(lines) >= 2:
            parts = lines[1].split()
            if len(parts) >= 2:
                try:
                    start_time = datetime.strptime(parts[1], "%Y-%m-%dT%H:%M:%S")
                    now = datetime.now()
                    return int((start_time - now).total_seconds())/3600
                except:
                    time.sleep(1)
        else: return 0
        time.sleep(1)
    raise TimeoutError("Failed to read estimated launch time.")

def SLURM_make_times_report():
    nodes = [1, 2]
    cpus_per_node = [1, 2, 4, 8, 16]
    mems = ["1G", "2G", "4G", "8G"]
    times = ["01:00:00", "02:00:00", "04:00:00", "08:00:00"]
    for node in nodes:
        for cpus in cpus_per_node:
            for mem in mems:
                for time in times:
                    slurm = Slurm(
                            job_name=str(node)+"-"+str(cpus)+"-"+str(mem),
                            nodes=node,
                            ntasks=node*cpus,
                            mem=mem,
                            time=time,
                            partition="core"  # From file?
                            )
                    job_id = slurm.sbatch()
                    wating_time = SLURM_get_start_time(job_id)
                    subprocess.run(["scancel", str(job_id)], check=True)
                    print(node, cpus, mem, time, wating_time)

if __name__ == "__main__":
    res = get_available_resources(printing=True)
#    for r in res:
#        print(f"{r['node']}: {r['free_cpu']} CPU free, {r['free_mem_MB']} MB RAM free")

