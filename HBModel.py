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
filename = r'PEO_LiTSFI_13.1wArgyrodite.xlsx'
 
# specify the sheet name
sheet_name = r'LithiumFWHM'
n = 29
 
# open the folder
os.chdir(folder)
 
# import the data from excel sheet
data = pd.read_excel(filename,sheet_name = sheet_name)
 
#convert data units
T = data['Temperature (°C)'] + 273
PW = 116.642*data['FWHM PEO (ppm)']
PW_error = 2.326*116.642*data['Error (ppm)']

T = T[:n]
PW = PW[:n]
PW_error = PW_error[:n]


# Fit the data
popt, pcov = curve_fit(HB_model, T, PW, p0=(6000, 1, 0.001617*250, 20), maxfev = 10000, bounds=([4000,10**-9,0.1,1], [8000, 100, 1.5,200]))

# Extract the fitted parameters
Delta_v_r_fit, B_fit, E_a_fit, D_fit = popt

# Create a temperature range for the fitted curve
T_fit = np.linspace(min(T), max(T), 100)

# Calculate the fitted curve
PW_fit = HB_model(T_fit, *popt)

fig, ax1 = plt.subplots(dpi=180)

# Plot the data with error bars
ax1.errorbar(T, PW, yerr=PW_error, fmt='o', label='Data', ecolor='gray', capsize=5, elinewidth=2)

# Plot the fitted curve
ax1.plot(T_fit, PW_fit, label='H-B. fit', color='red')

# Add secondary x-axis for Celsius
ax2 = ax1.twiny()
ax2.set_xlim((min(T)-273.15, max(T)-273.15))
ax2.set_xlabel("Temperature (°C)")

# Set labels and title
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Peak width (Hz)")
ax1.set_title('$PEO_{13}:LiTFSI$ Peak Width with HB fit')

plt.show()

# Extract the fitted parameters and their standard deviations
Delta_v_r_fit, B_fit, E_a_fit, D_fit = popt
Delta_v_r_std, B_std, E_a_std, D_std = np.sqrt(np.diag(pcov))

# Print the values and their uncertainties
print(f"Delta_v_r: {Delta_v_r_fit:.4f} ± {Delta_v_r_std:.4f}")
print(f"B: {B_fit:.4f} ± {B_std:.4f}")
print(f"E_a: {E_a_fit:.4f} ± {E_a_std:.4f}")
print(f"D: {D_fit:.4f} ± {D_std:.4f}")