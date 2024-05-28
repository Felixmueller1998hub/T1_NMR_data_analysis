import numpy as np
import matplotlib.pyplot as plt

# Instrumental specifications
beta_inst = 0.05 * (np.pi / 180)  # Convert instrumental broadening to radians

# Observed FWHM values and their errors (convert to radians)
FWHMSPE25 = np.array([0.268]) * (np.pi / 180)
FWHMSPE10 = np.array([0.294]) * (np.pi / 180)
FWHMCPE25 = np.array([0.183]) * (np.pi / 180)
FWHMCPE10 = np.array([0.200]) * (np.pi / 180)

er_obs = (np.array([0.004, 0.005, 0.002, 0.002]) + 0.026) * (np.pi / 180)

# Combine FWHM values into an array
FWHM_obs = np.concatenate((FWHMSPE25, FWHMSPE10, FWHMCPE25, FWHMCPE10))

# Function to calculate the sample broadening and its error
def calculate_sample_broadening(FWHM_obs, FWHM_inst, er_obs):
    beta_sample = np.sqrt(FWHM_obs**2 - FWHM_inst**2)
    beta_sample_error = np.sqrt((FWHM_obs * er_obs / beta_sample)**2)
    return beta_sample, beta_sample_error

# Calculate sample broadening and errors
beta_sample, beta_sample_error = calculate_sample_broadening(FWHM_obs, beta_inst, er_obs)

# Displaying the results
sample_names = ["SPE 25°C", "SPE -10°C", "CPE 25°C", "CPE -10°C"]
for i, sample in enumerate(sample_names):
    print(f"Sample {sample}:")
    print(f"  Observed FWHM: {FWHM_obs[i] * (180 / np.pi):.6f}°")
    print(f"  Sample Broadening: {beta_sample[i] * (180 / np.pi):.6f}°")
    print(f"  Sample Broadening Error: {beta_sample_error[i] * (180 / np.pi):.6f}°")

# Calculate grain size and errors for each sample
lam = 0.15418e-9
K = 0.9

def calculate_grain_size_and_error(theta, beta_sample, beta_sample_error):
    grain_size = (K * lam) / (beta_sample * np.cos(theta))
    grain_size_error = grain_size * np.sqrt((beta_sample_error / beta_sample)**2 + (np.tan(theta) * 0.5 * (np.pi / 180) / np.cos(theta))**2)
    return grain_size, grain_size_error

thetaSPE25 = np.array([19.60]) / 2 * (np.pi / 180)
thetaSPE10 = np.array([19.64]) / 2 * (np.pi / 180)
thetaCPE25 = np.array([19.63]) / 2 * (np.pi / 180)
thetaCPE10 = np.array([19.67]) / 2 * (np.pi / 180)

theta_values = [thetaSPE25, thetaSPE10, thetaCPE25, thetaCPE10]

# Calculating grain size and errors
grain_sizes = []
grain_size_errors = []

for i, theta in enumerate(theta_values):
    grain_size, grain_size_error = calculate_grain_size_and_error(theta, beta_sample[i], beta_sample_error[i])
    grain_sizes.append(grain_size[0] * 1e9)  # Convert to nm
    grain_size_errors.append(grain_size_error[0] * 1e9)  # Convert to nm
    print(f"Grain Size for Sample {sample_names[i]}: {grain_size[0]*1e9:.2f} nm ± {grain_size_error[0]*1e9:.2f} nm")

# Plotting the results
plt.figure(figsize=(6, 4))
plt.bar(sample_names, grain_sizes, yerr=grain_size_errors, capsize=5, color='skyblue')
plt.xlabel('Sample')
plt.ylabel('Grain Size (nm)')
plt.title('Grain Size Using Scherrers equation')
plt.show()