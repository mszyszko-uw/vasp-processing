import subprocess
import re

def get_available_resources():
    # Pobierz dane z `scontrol show node`
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

    return resources

# Przykładowe użycie:
if __name__ == "__main__":
    res = get_available_resources()
    for r in res:
        print(f"{r['node']}: {r['free_cpu']} CPU free, {r['free_mem_MB']} MB RAM free")

