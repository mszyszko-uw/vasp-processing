import subprocess

def run_local(steps_def, step, vasp_command="mpiexec vasp_std", folder='./'):
	command_to_run = ""

	slurm_commands = steps_def[step]['cmd']
	for line in slurm_commands .splitlines():	
		if line.strip().startswith("module"):
			command_to_run += ";"+line

	command_to_run += '; cd '+ folder
	command_to_run += "; python3 Preprocessing/toolkit.py --step t_" + step + " --part dry"
	
	command_to_run += ";" + vasp_command 
	command_to_run += ";" + "cd .." 
	try:	
		subprocess.run(command_to_run)
	except:
		print("Cannot run:\n " + command_to_run+"\n")