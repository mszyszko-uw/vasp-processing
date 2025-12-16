import os
from pathlib import Path
import math

#def run_dry():


# Exctracts NKPTS and NBANDS from OUTCAR - for use after dry run
#This is a copy from Preprocessing/toolkit.py
def get_NKPTS_NBANDS(file):
    nkpts, nbands = None, None
    with open(file, "r") as f:
        for line in f:
            if "NKPTS" in line and "NBANDS" in line: 
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "NKPTS":
                        nkpts = int(parts[i+2]) 
                    if "NBANDS" in p:
                        if p == 'NBANDS':
                            nbands = int(parts[i+2])
                        elif p =='NBANDS=':
                            nbands = int(parts[i+1])
                break
    return nkpts, nbands
def divisors(n):
    divs = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)
def get_parallelization(nkpts, nbands, MAX_CPUs):    
    # Possible NCORE: common divisors of nbands and cpu_per_socket
    possible_ncore = [d for d in divisors(nbands) if d in divisors(MAX_CPUs)]
    possible_kpar = divisors(nkpts)

    best = None
    for ncore in possible_ncore:
        for kpar in possible_kpar:
            if ncore * kpar <= MAX_CPUs:
                # heurystics: highest number of CPU
                score = ncore * kpar
                if best is None or score > best[2]:
                    best = (ncore, kpar, score)
    if best:
        return score
    else:
        return None

def find_last_out(root_folder):
    latest_file = None
    latest_time = 0
    for outcar_path in Path(root_folder).rglob("OUTCAR"):
        mod_time = outcar_path.stat().st_mtime
        if mod_time > latest_time:
            latest_time = mod_time
            latest_file = outcar_path
    if latest_file:
        return latest_file.parent
    else:
        return None

def set_parallelization(step_folder='./', MAX_CPUs=32):
		folder = find_last_out("./"+ step_folder)

#		os.path.join(folder, '00_t_dry')	
		NKPTS, NBANDS = get_NKPTS_NBANDS(os.path.join(folder, 'OUTCAR'))
		NCPUs = get_parallelization(NKPTS, NBANDS, MAX_CPUs)
		return NCPUs;


if __name__ =="__main__":
	if(set_parallelization()):
		print(set_parallelization())
