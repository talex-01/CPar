import matplotlib.pyplot as plt
import numpy as np

# Data from the user's provided output
omp_num_threads = list(range(1, 49))  # Threads from 1 to 48
real_times = [
    15.52, 6.7167, 4.69013, 3.6784, 3.1249, 2.6314, 2.3620, 2.1393, 
    2.0062, 1.8499, 1.7956, 1.6601, 1.5778, 1.5285, 1.5687, 1.4996, 
    1.5124, 1.4262, 1.4686, 1.4535, 1.5550, 1.5039, 1.4626, 1.5082, 
    1.9247, 1.856, 1.9143, 1.866, 1.961, 2.0188, 1.887, 2.186, 
    2.1866, 1.907, 2.236, 2.341, 2.134, 2.027, 2.143, 2.339, 
    2.223, 2.532, 2.088, 2.457, 2.096, 2.434, 2.126, 2.296
]

# Calculate the observed speedup
max_time = max(real_times)  # Time when OMP_NUM_THREADS = 1
observed_speedup = max_time / np.array(real_times)

# Create the theoretical speedup function
linear_threshold = 24
# Linear speedup up to 24 threads
linear_speedup = np.minimum(omp_num_threads[:linear_threshold], linear_threshold)
# Constant speedup after 24 threads
constant_speedup = np.ones(len(omp_num_threads) - linear_threshold) * linear_threshold
theoretical_speedup = np.concatenate([linear_speedup, constant_speedup])

# Plotting both the observed and theoretical speedup
plt.figure(figsize=(14, 7))

# Plot observed speedup
plt.plot(omp_num_threads, observed_speedup, marker='o', linestyle='-', color='b', label='Observed Speedup')

# Plot theoretical speedup (linear to constant)
plt.plot(omp_num_threads, theoretical_speedup, linestyle='--', color='r', label='Theoretical Speedup (Linear to Constant)')

# Title and labels
plt.title('Speedup vs OMP_NUM_THREADS')
plt.xlabel('OMP_NUM_THREADS')
plt.ylabel('Speedup')
plt.grid(True)

# Adjust tick marks for readability
plt.xticks(range(1, 49, 2))  # Show every 1st thread for readability
plt.yticks(range(1, 25, 2))

# Add legend
plt.legend()

# Show plot
plt.show()
