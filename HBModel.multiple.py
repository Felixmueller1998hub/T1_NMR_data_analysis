# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:43:47 2023
FWHM fitting
@author: felix
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os as os
import pandas as pd

# Define the equation you want to fit
def HB_model(T, Delta_v_r, B, E_a, D):
    return Delta_v_r * (1 + ((Delta_v_r / B) - 1) * np.exp(-E_a / (8.6173*10**(-5) * T)))**(-1) + D

def HB_VTF(T, Delta_v_r, B, E_a, D, T_0):
    return Delta_v_r * (1 + ((Delta_v_r / B) - 1) * np.exp(-E_a / (8.6173*10**(-5) *(T-T_0))))**(-1) + D


# Data
# specify the folder
folder = r'C:\Users\felix\OneDrive\Desktop\Master Thesis'

# specify the filename of the datafile
filename = r'PEO_LiTSFI_13.1wArgyrodite.xlsx'
filename2 = r'PEO_LiTSFI_13.1wArgyrodite_dry.xlsx'
filename3 = r'PEO_LiTSFI_13.1nofiller.xlsx'
 
# specify the sheet name
sheet_name = r'LithiumFWHM'
n = 28
 
# open the folder
os.chdir(folder)
 
# import the data
data = pd.read_excel(filename,sheet_name = sheet_name)
data2 = pd.read_excel(filename2,sheet_name = sheet_name)
data3 = pd.read_excel(filename3,sheet_name = sheet_name)

T = data['Temperature (°C)'] + 273
PW = 116.642*data['FWHM PEO (ppm)']
PW_error = 2.326*116.642*data['Error (ppm)']
T2 = data2['Temperature (°C)'] + 273
PW2 = 116.642*data2['FWHM PEO (ppm)']
PW_error2 = 2.326*116.642*data2['Error (ppm)']
T3 = data3['Temperature (°C)'] + 273
PW3 = 116.642*data3['FWHM PEO (ppm)']
PW_error3 = 2.326*116.642*data3['Error (ppm)']

T = T[:n]
PW = PW[:n]
PW_error = PW_error[:n]
T2 = T2[:n]
PW2= PW2[:n]
PW_error2 = PW_error2[:n]

# Create a figure with two subplots, one on top of the other
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=180, figsize=(8, 8))

# For the first dataset
popt, pcov = curve_fit(HB_model, T, PW,  p0=[6000, 1, 0.4, 20], maxfev=10000, bounds=([4000, 10**-13, 0.3, 1], [8000, 100, 0.6, 200]))

# For the second dataset
popt2, pcov2 = curve_fit(HB_model, T2, PW2,  p0=[6000, 0.1, 0.001617*250, 20], maxfev=10000, bounds=([4000, 10**-13, 0.4, 1], [8000, 100, 1.5, 200]))

# For the third dataset
popt3, pcov3 = curve_fit(HB_VTF, T3, PW3,  p0=[6000, 0.1, 0.01, 200, 200], maxfev=10000, bounds=([5000, 10**-13, 0.01, 1, 0], [8000, 100, 2, 200, 250]))
popt4, pcov4 = curve_fit(HB_model, T3, PW3,  p0=[6000, 0.1, 0.001617*250, 20], maxfev=10000, bounds=([4000, 10**-13, 0.4, 1], [8000, 100, 1.5, 200]))


# Fit and plot the first dataset
T_fit = np.linspace(min(T), max(T), 100)
PW_fit = HB_model(T_fit, *popt)
ax1.errorbar(T, PW, yerr=PW_error, fmt='o', label='CPE Solvent 13:1', ecolor='grey', capsize=5, elinewidth=2, color='red')
ax1.plot(T_fit, PW_fit, 'g--', label='H-B. fit')  # Make the fit curve dashed
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Peak width (Hz)")
ax1.set_title('$PEO_{13}:LiTFSI$ PEO line width')
ax1.set_xlim(210,360)
ax1.text(350, 2900, '$^{7}Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax1.text(270, 2600, f'$E_a$ = {popt[2]:.2f}± 0.02 eV', fontsize=12, color='black', verticalalignment='top')

ax1.legend()

# Fit and plot the second dataset
T_fit2 = np.linspace(min(T2), max(T2), 100)
PW_fit2 = HB_model(T_fit2, *popt2)
ax2.errorbar(T2, PW2, yerr=PW_error2, fmt='o', label='CPE Dry 13:1', ecolor='gray', capsize=5, elinewidth=2, color='maroon', markerfacecolor = 'none')
ax2.plot(T_fit2, PW_fit2, 'g--', label='H-B. fit')  # Make the fit curve dashed
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Peak width (Hz)")
#ax2.set_title('$PEO_{13}:LiTFSI$ Dry Method')
ax2.set_xlim(210,360)
#ax2.text(350, 3200, '$7Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax2.text(275, 2600, f'$E_a$ = {popt2[2]:.2f}± 0.02 eV', fontsize=12, color='black', verticalalignment='top')

ax2.legend()

# Fit and plot the second dataset
T_fit3 = np.linspace(min(T3), max(T3), 100)
PW_fit3 = HB_VTF(T_fit3, *popt3)
PW_fit4 = HB_model(T_fit3, *popt4)
ax3.errorbar(T3, PW3, yerr=PW_error3, fmt='o', label='SPE 13:1', ecolor='gray', capsize=5, elinewidth=2, color='orange')
ax3.plot(T_fit3, PW_fit3, 'r--', label='VTF H-B. fit')  # Make the fit curve dashed
ax3.plot(T_fit3, PW_fit4, 'b--', label='Arrhenius H-B. fit')  # Make the fit curve dashed
ax3.set_xlabel("Temperature (K)")
ax3.set_ylabel("Peak width (Hz)")
#ax3.set_title('$PEO_{13}:LiTFSI$ no filler')
ax3.set_xlim(210,360)
#ax3.text(350, 2800, '$7Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax3.text(280, 2600, f'$B$ = {popt3[2]:.3f}± 0.005', fontsize=12, color='black', verticalalignment='top')

ax3.legend()


# Adjust layout
plt.tight_layout()
plt.show()

# Print the fit results for both datasets
print("Fitted parameters for the first dataset:")
print(f"Delta_v_r: {popt[0]:.4f} ± {np.sqrt(np.diag(pcov))[0]:.4f}")
print(f"B: {popt[1]:.10f} ± {np.sqrt(np.diag(pcov))[1]:.10f}")
print(f"E_a: {popt[2]:.4f} ± {np.sqrt(np.diag(pcov))[2]:.4f}")
print(f"D: {popt[3]:.4f} ± {np.sqrt(np.diag(pcov))[3]:.4f}")

print("\nFitted parameters for the second dataset:")
print(f"Delta_v_r: {popt2[0]:.4f} ± {np.sqrt(np.diag(pcov2))[0]:.4f}")
print(f"B: {popt2[1]:.10f} ± {np.sqrt(np.diag(pcov2))[1]:.10f}")
print(f"E_a: {popt2[2]:.4f} ± {np.sqrt(np.diag(pcov2))[2]:.4f}")
print(f"D: {popt2[3]:.4f} ± {np.sqrt(np.diag(pcov2))[3]:.4f}")

print("\nFitted parameters for the third dataset:")
print(f"Delta_v_r: {popt3[0]:.4f} ± {np.sqrt(np.diag(pcov3))[0]:.4f}")
print(f"B: {popt3[1]:.4f} ± {np.sqrt(np.diag(pcov3))[1]:.4f}")
print(f"E_a: {popt3[2]:.4f} ± {np.sqrt(np.diag(pcov3))[2]:.4f}")
print(f"D: {popt3[3]:.4f} ± {np.sqrt(np.diag(pcov3))[3]:.4f}")
print(f"t_0: {popt3[4]:.4f} ± {np.sqrt(np.diag(pcov3))[4]:.4f}")