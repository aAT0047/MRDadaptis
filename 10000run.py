import os
import subprocess
from tqdm import tqdm

work_dir = "/home/cloudam/GSDcreator/sh_files"
script_count = 20

for i in tqdm(range(1, script_count + 1), desc="Executing scripts", unit="script"):
    script_name = "A_stableCallerPaperSimFlowShell_{}.sh".format(i)
    script_path = os.path.join(work_dir, script_name)

    # Check if the script exists
    if os.path.isfile(script_path):
        try:
            # Execute the shell script
            subprocess.check_call(['bash', script_path])
            print("Successfully executed {}".format(script_name))
        except subprocess.CalledProcessError:
            print("Error in executing {}".format(script_name))
    else:
        print("{} does not exist!".format(script_name))

print("\nAll scripts executed!")
