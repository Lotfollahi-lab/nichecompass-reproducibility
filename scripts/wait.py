import os
import time

# Define the command to check whether the first script is still running
CHECK_COMMAND = "pgrep -f train_nichecompass_benchmarking_models.py"
CHECK_COMMAND_2 = "pgrep -f train_nichecompass_reference_model.py"

# Wait for the first script to finish
while True:
    if (os.system(CHECK_COMMAND) != 0) & (os.system(CHECK_COMMAND_2) != 0):
        print("First script finished. Launching follow up script.")
        break
    else:
        print("First script still running...")
    time.sleep(1800)