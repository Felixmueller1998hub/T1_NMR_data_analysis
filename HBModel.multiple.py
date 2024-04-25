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


# Data
# specify the folder
folder = r'C:\Users\felix\OneDrive\Desktop\Master Thesis'

# specify the filename of the datafile
filename = r'PEO_LiTSFI_18.1wArgyrodite.xlsx'
filename2 = r'PEO_LiTSFI_18.1wArgyrodite_dry.xlsx'
filename3 = r'PEO_LiTSFI_18.1nofiller.xlsx'
 
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
popt3, pcov3 = curve_fit(HB_model, T3, PW3,  p0=[6000, 0.1, 1, 20], maxfev=10000, bounds=([100, 10**-13, 0.5, 1], [8000, 100, 2, 200]))


# Fit and plot the first dataset
T_fit = np.linspace(min(T), max(T), 100)
PW_fit = HB_model(T_fit, *popt)
ax1.errorbar(T, PW, yerr=PW_error, fmt='o', label='18.1 w/LPSC Solvent Method', ecolor='grey', capsize=5, elinewidth=2, color='olive')
ax1.plot(T_fit, PW_fit, 'g--', label='H-B. fit')  # Make the fit curve dashed
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Peak width (Hz)")
ax1.set_title('$PEO_{13}:LiTFSI$ Solvent Method')
ax1.set_xlim(210,360)
ax1.text(350, 2900, '$7Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax1.text(270, 2600, f'E_a = {popt[2]:.2f}± 0.02 eV', fontsize=12, color='black', verticalalignment='top')

ax1.legend()

# Fit and plot the second dataset
T_fit2 = np.linspace(min(T2), max(T2), 100)
PW_fit2 = HB_model(T_fit2, *popt2)
ax2.errorbar(T2, PW2, yerr=PW_error2, fmt='o', label='18.1 w/LPSC Dry Method', ecolor='gray', capsize=5, elinewidth=2, color='darkgreen', markerfacecolor = 'none')
ax2.plot(T_fit2, PW_fit2, 'g--', label='H-B. fit')  # Make the fit curve dashed
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Peak width (Hz)")
ax2.set_title('$PEO_{13}:LiTFSI$ Dry Method')
ax2.set_xlim(210,360)
ax2.text(350, 3200, '$7Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax2.text(275, 2600, f'E_a = {popt2[2]:.2f}± 0.03 eV', fontsize=12, color='black', verticalalignment='top')

ax2.legend()

# Fit and plot the second dataset
T_fit3 = np.linspace(min(T3), max(T3), 100)
PW_fit3 = HB_model(T_fit3, *popt3)
ax3.errorbar(T3, PW3, yerr=PW_error3, fmt='o', label='18.1 no filler', ecolor='gray', capsize=5, elinewidth=2, color='blue')
ax3.plot(T_fit3, PW_fit3, 'b--', label='H-B. fit')  # Make the fit curve dashed
ax3.set_xlabel("Temperature (K)")
ax3.set_ylabel("Peak width (Hz)")
ax3.set_title('$PEO_{13}:LiTFSI$ no filler')
ax3.set_xlim(210,360)
ax3.text(350, 2800, '$7Li$ \n 116.6 MHz', fontsize=12, color='black', ha='center')
# Add E_a value as annotation
ax3.text(280, 2600, f'E_a = {popt3[2]:.2f}± 0.03 eV', fontsize=12, color='black', verticalalignment='top')

ax3.legend()


# Adjust layout
plt.tight_layout()
plt.show()

# Print the fit results for both datasets
print("Fitted parameters for the first dataset:")
print(f"Delta_v_r: {popt[0]:.4f} ± {np.sqrt(np.diag(pcov))[0]:.4f}")
print(f"B: {popt[1]:.4f} ± {np.sqrt(np.diag(pcov))[1]:.4f}")
print(f"E_a: {popt[2]:.4f} ± {np.sqrt(np.diag(pcov))[2]:.4f}")
print(f"D: {popt[3]:.4f} ± {np.sqrt(np.diag(pcov))[3]:.4f}")

print("\nFitted parameters for the second dataset:")
print(f"Delta_v_r: {popt2[0]:.4f} ± {np.sqrt(np.diag(pcov2))[0]:.4f}")
print(f"B: {popt2[1]:.4f} ± {np.sqrt(np.diag(pcov2))[1]:.4f}")
print(f"E_a: {popt2[2]:.4f} ± {np.sqrt(np.diag(pcov2))[2]:.4f}")
print(f"D: {popt2[3]:.4f} ± {np.sqrt(np.diag(pcov2))[3]:.4f}")

print("\nFitted parameters for the third dataset:")
print(f"Delta_v_r: {popt3[0]:.4f} ± {np.sqrt(np.diag(pcov3))[0]:.4f}")
print(f"B: {popt3[1]:.4f} ± {np.sqrt(np.diag(pcov3))[1]:.4f}")
print(f"E_a: {popt3[2]:.4f} ± {np.sqrt(np.diag(pcov3))[2]:.4f}")
print(f"D: {popt3[3]:.4f} ± {np.sqrt(np.diag(pcov3))[3]:.4f}")
#%%

def calculate_r_squared(T, PW, popt, model):
    """
    Calculate the R-squared statistic for model evaluation.
    
    Parameters:
    - T (array): Independent variable data (e.g., Temperature).
    - PW (array): Observed dependent variable data (e.g., Peak Width).
    - popt (list): Optimal parameters from the curve fitting.
    - model (function): Model function used for curve fitting.
    
    Returns:
    - float: R-squared value.
    """
    # Predict values using the fitted model
    PW_pred = model(T, *popt)
    
    # Calculate the residual sum of squares (SS_res)
    SS_res = np.sum((PW - PW_pred)**2)
    
    # Calculate the total sum of squares (SS_tot)
    SS_tot = np.sum((PW - np.mean(PW))**2)
    
    # Calculate R-squared
    R_squared = 1 - (SS_res / SS_tot)
    
    return R_squared

# Calculate R-squared for the first dataset
R_squared1 = calculate_r_squared(T, PW, popt, HB_model)
print("R-squared for Dataset 1:", R_squared1)

# Calculate R-squared for the second dataset
R_squared2 = calculate_r_squared(T2, PW2, popt2, HB_model)
print("R-squared for Dataset 2:", R_squared2)

# Calculate R-squared for the third dataset
R_squared3 = calculate_r_squared(T3, PW3, popt3, HB_model)
print("R-squared for Dataset 3:", R_squared3)