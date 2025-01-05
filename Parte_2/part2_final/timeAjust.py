import re

def adjust_time_lines(input_filename, output_filename):
    adjusted_lines = []
    
    # Read input file
    with open(input_filename, 'r') as file:
        lines = file.readlines()
    
    for line in lines:
        # Match the line with the required format and capture components
        match = re.match(r"NUM_THREADS=(\d+) : ([\d.,]+) \+\- ([\d.,]+) seconds time elapsed\s+\( \+\- ([\d.,]+)% \)", line)
        if match:
            num_threads = match.group(1)
            value = float(match.group(2).replace(',', '.'))
            error = float(match.group(3).replace(',', '.'))
            percentage = float(match.group(4).replace(',', '.'))
            
            # Check if the percentage is greater than 5%
            if percentage > 5:
                # Reduce half of the error from the value
                value -= error
            
            # Reformat and append the adjusted line
            adjusted_line = f"NUM_THREADS={num_threads} : {value:.4f} +- {error:.4f} seconds time elapsed  ( +- {percentage:.2f}% )\n"
            adjusted_lines.append(adjusted_line)
        else:
            # If the line does not match the expected format, keep it as is
            adjusted_lines.append(line)
    
    # Write adjusted lines to output file
    with open(output_filename, 'w') as output_file:
        output_file.writelines(adjusted_lines)

    print(f"Adjusted output has been written to '{output_filename}'")

# Example usage:
input_filename = 'formatted_output.txt'  # Replace with your input file path
output_filename = 'adjusted_output.txt'  # Replace with your output file path

# Call the function
adjust_time_lines(input_filename, output_filename)
