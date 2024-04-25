import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

k = 8.6173e-5  # Boltzmann constant in eV/K

def Arrhenius(T, tau_0, E_a):
    return tau_0 * np.exp(E_a / (k * T))

def VTF(T, tau_0, E_a, T0):
    return tau_0 * np.exp(E_a / (k *(T-T0)))

def DSE(T, n, V):
    return (n*V)/(k*T)

def abragam(Delta_Hf, Delta_RL):
    return (0.1 / Delta_Hf) * np.tan(np.pi / 2 * (Delta_Hf / Delta_RL) ** 2)

def d_tau_c_d_Delta_Hf(Delta_Hf, Delta_RL):
    term1 = -1 / Delta_Hf**2 * np.tan(np.pi / 2 * (Delta_Hf / Delta_RL)**2)
    term2 = (np.pi / Delta_RL * 2 * Delta_Hf / Delta_RL * (1 / np.cos(np.pi / 2 * (Delta_Hf / Delta_RL)**2))**2) / Delta_Hf
    return term1 + term2

def d_tau_c_d_Delta_RL(Delta_Hf, Delta_RL):
    term = -np.pi * Delta_Hf / Delta_RL**3 * 2 * Delta_Hf / Delta_RL * (1 / np.cos(np.pi / 2 * (Delta_Hf / Delta_RL)**2))**2
    return term

def calculate_error_in_tau_c(Delta_Hf, Delta_RL, sigma_Delta_Hf, sigma_Delta_RL):
    sigma_tau_c = np.abs(np.sqrt((d_tau_c_d_Delta_Hf(Delta_Hf, Delta_RL) * sigma_Delta_Hf)**2 +
                          (d_tau_c_d_Delta_RL(Delta_Hf, Delta_RL) * sigma_Delta_RL)**2))
    return sigma_tau_c


#Take these values from HB_modelfit, where both are given as ouput
Delta_RL1 = 5387.79 # Single value for Delta_RL
sigma_RL1 = 80.7  # Uncertainty in Delta_RL
Delta_RL2 = 6262.12 # Single value for Delta_RL
sigma_RL2 = 119.1  # Uncertainty in Delta_RL

# Data import from an excel sheet
folder = r'C:\Users\felix\OneDrive\Desktop\Master Thesis'
filename1 = r'PEO_LiTSFI_18.1wArgyrodite.xlsx'
filename2 = r'PEO_LiTSFI_18.1wArgyrodite_dry.xlsx'
sheet_name = r'LithiumFWHM'

Title = 'PEO:LiTFSI 18:1 / 10wt.% $Li_6PS_5Cl$'
Dataset1 = 'solvent'
Dataset2 = 'Dry'

# Load data for both datasets
data1 = pd.read_excel(f'{folder}/{filename1}', sheet_name=sheet_name)
data2 = pd.read_excel(f'{folder}/{filename2}', sheet_name=sheet_name)

# Prepare data for first dataset
T1 = data1['Temperature (°C)'] + 273    #convert from celcius to kelvin
PW1 = 116.642 * data1['FWHM PEO (ppm)']     #convert from ppm to Hz
errordata1 = data1['Error (ppm)'] * 116.642*1.96    #ppm to Hz and 95%CI
p1 = 24
errordata1 = errordata1[:p1]
PW1 = PW1[:p1]
T1 = T1[:p1]
tau_array1 = abs(abragam(PW1, Delta_RL1))

# Prepare data for second dataset
T2 = data2['Temperature (°C)'] + 273
PW2 = 116.642 * data2['FWHM PEO (ppm)']
errordata2 = data2['Error (ppm)'] * 116.642*1.96
p2 = 24
errordata2 = errordata2[:p2]
PW2 = PW2[:p2]
T2 = T2[:p2]
tau_array2 = abs(abragam(PW2, Delta_RL2))

#ERRORS
# Calculate error in tau_c for each Delta_Hf point
errors_tau_c1 = calculate_error_in_tau_c(PW1, Delta_RL1, errordata1, sigma_RL1)
print("Errors in tau_c data1:", errors_tau_c1)
errors_tau_c2 = calculate_error_in_tau_c(PW2, Delta_RL2, errordata2, sigma_RL2)
print("Errors in tau_c data2:", errors_tau_c2)
errors_tau_c1 = abs(errors_tau_c1)
errors_tau_c2 = abs(errors_tau_c2)


#Fit Arrhenius region
n = 24
m = 14
r1 = T1[m:n]
r2 = T1[m+1:n]
tau_r1 = tau_array1[m:n]
tau_r2 = tau_array2[m+1:n]
error1 = errors_tau_c1[m:n]
error2 = errors_tau_c2[m+1:n]


#FIT VTF region
st = 5
r3 = T1[st:m+1]
r4 = T1[st:m+3]
error3 = errors_tau_c1[st:m+1]  # Assuming error1 has been calculated as shown previously
error4 = errors_tau_c2[st:m+3]  # Assuming error2
tau_r3 = tau_array1[st:m+1]
tau_r4 = tau_array2[st:m+3]

# Initial guess for parameters tau_0 and E_a
initial_guess = [1e-9, 0.5]  # Example values; adjust based on expected range


# Fit Arrhenius equation to the first dataset
popt1, pcov1 = curve_fit(Arrhenius, r1, tau_r1, p0=initial_guess, sigma=error1, absolute_sigma=True, maxfev = 8000)
tau_0_fitted1, E_a_fitted1 = popt1

# Fit Arrhenius equation to the second dataset
popt2, pcov2 = curve_fit(Arrhenius, r2, tau_r2, p0=initial_guess, sigma=error2, absolute_sigma=True, maxfev = 8000)
tau_0_fitted2, E_a_fitted2 = popt2

# Fit VTF equation to the first dataset
popt3, pcov3 = curve_fit(VTF, r3, tau_r3, p0=[1e-10, 0.1, 203], sigma=error3, absolute_sigma=True, maxfev = 8000)
tau_0_fitted3, E_a_fitted3, T03 = popt3

# Fit VTF equation to the second dataset
popt4, pcov4 = curve_fit(VTF, r4, tau_r4, p0=[1e-10, 0.1, 203], sigma=error4, absolute_sigma=True, maxfev = 8000)
tau_0_fitted4, E_a_fitted4, T04 = popt4


# Generate temperatures for plotting fitted curves
T_r1 = np.linspace(min(r1), max(r1), 100)
T_r2 = np.linspace(min(r2), max(r2), 100)
T_r3 = np.linspace(min(r3), max(r3), 100)
T_r4 = np.linspace(min(r4), max(r4), 100)


# Plotting both datasets and fits
fig, ax1 = plt.subplots(dpi=180)
#data
plt.errorbar(1000/(T1), (tau_array1), yerr=errors_tau_c1, markerfacecolor = 'none', markeredgecolor = 'red', fmt='o', label=f'{Dataset1}', ecolor='grey', capsize=3)
plt.errorbar(1000/(T2), (tau_array2), yerr=errors_tau_c2, markerfacecolor = 'none', markeredgecolor = 'maroon', fmt='o', label=f'{Dataset2}', ecolor='grey', capsize=3)
ax1.plot(1000/(T_r1), (Arrhenius(T_r1, *popt1)), 'y--', label='Arrhenius fit')
ax1.plot(1000/(T_r2), (Arrhenius(T_r2, *popt2)), 'r--',  label='Arrhenius fit')
#second region
ax1.plot(1000/(T_r3), (VTF(T_r3, *popt3)), 'y-', label='VTF fit')
ax1.plot(1000/(T_r4), (VTF(T_r4, *popt4)), 'r-',  label='VTF fit')

ax1.set_title(f'{Title}')
ax1.set_xlabel('1000/T ($K^{-1}$)')
ax1.set_ylabel(r'log($\tau_c$) (s)')
plt.yscale('log')
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

# Choose specific Celsius values from 80°C to -20°C, every 10 degrees
celsius_values = np.arange(80, -46, -20)
# Convert Celsius to the corresponding 1/kbT values
corresponding_1_over_kbT = 1000 / ((celsius_values + 273) )

ax2.set_xticks(corresponding_1_over_kbT)
ax2.set_xticklabels([f'{int(t)}°C' for t in celsius_values])
ax2.set_xlabel('Temperature (°C)')

ax1.legend()
plt.show()

perr1 = np.sqrt(np.diag(pcov1))
perr2 = np.sqrt(np.diag(pcov2))
perr3 = np.sqrt(np.diag(pcov3))
perr4 = np.sqrt(np.diag(pcov4))

# Output fitted parameters
print(f"1 Arrhenius region - Fitted tau_0: {tau_0_fitted1:.2e} s, Estimated Activation Energy: {E_a_fitted1:.4f}+- {perr1[1]:.4f} eV")
print(f"2 Arrhenius region- Fitted tau_0: {tau_0_fitted2:.2e} s, Estimated Activation Energy: {E_a_fitted2:.4f}+-{perr2[1]:.4f} eV")
print(f"1 VTF region - Fitted tau_0: {tau_0_fitted3:.2e} s, Estimated Activation Energy: {E_a_fitted3:.4f}+-{perr3[1]:.4f} eV, T_0: {T03:.4f} K")
print(f"2 VTF region- Fitted tau_0: {tau_0_fitted4:.2e} s, Estimated Activation Energy: {E_a_fitted4:.4f}+-{perr4[1]:.4f} eV,  T_0: {T04:.4f} K")
