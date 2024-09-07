import numpy as np
import scipy.stats as stats

def calculate_mean_and_ci(data, confidence=0.95):
    mean = np.mean(data)
    std_dev = np.std(data, ddof=1)
    n = len(data)
    
    # Calculate the confidence interval
    z = stats.norm.ppf((1 + confidence) / 2)
    margin_of_error = z * (std_dev / np.sqrt(n))
    
    ci_lower = mean - margin_of_error
    ci_upper = mean + margin_of_error
    
    return mean, (ci_lower, ci_upper)

def relative_ci_change(ci, mean):
    return [(ci[0] - mean) / mean, (ci[1] - mean) / mean]

# Data
A = [1, 2, 4, 5, 6, 0]
B = [1, 3, 4, 5, 7, 0]

# Calculate mean and CI for A
mean_A, ci_A = calculate_mean_and_ci(A)
relative_ci_A = relative_ci_change(ci_A, mean_A)

# Calculate mean and CI for B
mean_B, ci_B = calculate_mean_and_ci(B)
relative_ci_B = relative_ci_change(ci_B, mean_B)

print(f"Mean and CI for A: {mean_A:.2f}, {ci_A}")
print(f"Relative CI for A: {relative_ci_A}")
print()
print(f"Mean and CI for B: {mean_B:.2f}, {ci_B}")
print(f"Relative CI for B: {relative_ci_B}")

# Check if B is within 5% of A's CI
def is_within_5_percent(relative_ci_A, relative_ci_B):
    return all(abs(b) <= 0.05 for b in relative_ci_B)

if is_within_5_percent(relative_ci_A, relative_ci_B):
    print("Dãy B gần giống với dãy A theo CI.")
else:
    print("Dãy B không gần giống với dãy A theo CI.")
