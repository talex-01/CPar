import re

# Open the input file and read the contents
with open('timesin', 'r') as file:
    lines = file.readlines()

# Prepare a list to store the formatted output
formatted_lines = []

# Loop through each line and extract the relevant information
num_threads = None
for line in lines:
    # Check for the line indicating the number of threads
    match_threads = re.match(r"^Running with OMP_NUM_THREADS=(\d+)", line)
    if match_threads:
        num_threads = match_threads.group(1)  # Capture the number of threads
    # Check for the line containing the time elapsed
    if 'seconds time elapsed' in line and num_threads:
        # Format and append the output line with all details from the time line
        formatted_lines.append(f"NUM_THREADS={num_threads} : {line.strip()}\n")
        num_threads = None  # Reset for the next block

# Write the formatted content to a new file
with open('formatted_output.txt', 'w') as output_file:
    output_file.writelines(formatted_lines)

print("Formatted output has been written to 'formatted_output.txt'")
