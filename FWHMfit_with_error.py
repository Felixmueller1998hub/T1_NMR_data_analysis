# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:35:03 2023

@author: felix
"""
import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from sklearn.linear_model import LinearRegression

def calculate_and_display_slope(T, Y, m, n, pos, neg):
    selected_T = T[m:n].values.reshape(-1, 1)
    selected_Y = Y[m:n].values

    regressor = LinearRegression()
    regressor.fit(selected_T, selected_Y)
    
    # Create an extended range for T for plotting
    T_min, T_max = selected_T.min(), selected_T.max()
    T_range = T_max - T_min
    extended_T_min = T_min - pos * T_range
    extended_T_max = T_max + neg * T_range
   
    extended_T = np.linspace(extended_T_min, extended_T_max, 500).reshape(-1, 1)
    fit_line = regressor.predict(extended_T)
    
    # Plotting the linear fit line extended beyond the selected points
    ax1.plot(extended_T, fit_line, linestyle=':', color='black')


# Data
# specify the folder
folder = r'C:\Users\felix\OneDrive\Desktop\Master Thesis'

# specify the filename of the datafile
filename1 = r'PEO_LiTSFI_18.1nofiller.xlsx'
filename2 = r'PEO_LiTSFI_18.1wArgyrodite.xlsx'
filename3 = r'PEO_LiTSFI_18.1wArgyrodite_dry.xlsx'
 
# specify the sheet name
sheet_name = r'LithiumFWHM'
 
# open the folder
os.chdir(folder)
 
# import the data
data1 = pd.read_excel(filename1,sheet_name = sheet_name)
data2 = pd.read_excel(filename2,sheet_name = sheet_name)
data3 = pd.read_excel(filename3,sheet_name = sheet_name)

#X and Y
T1 = data1['Temperature (°C)'] + 273
PW1 = 116.642*data1['FWHM PEO (ppm)']
T2 = data2['Temperature (°C)'] + 273
PW2 = 116.642*data2['FWHM PEO (ppm)']
T3 = data3['Temperature (°C)'] + 273
PW3 = 116.642*data3['FWHM PEO (ppm)']
p = 14

# Assuming constant error for demonstration; replace with your actual error data
PW_error1 = 2.326*116.642*data1['Error (ppm)']  
PW_error2 = 2.326*116.642*data2['Error (ppm)']  
PW_error3 = 2.326*116.642*data3['Error (ppm)']  

fig, ax1 = plt.subplots(dpi=180)

# Plot the data with specified markers and styles
# Triangle markers for all, with the third dataset having open markers
ax1.errorbar(T1[:p], PW1[:p], yerr = PW_error1[:p], fmt='o', ecolor='gray', capsize=5, elinewidth=2, color='orange', label='no filler')
ax1.errorbar(T2[:p], PW2[:p], yerr = PW_error2[:p],fmt='o', ecolor='grey', label='solvent method', color = 'red', capsize=5, elinewidth=2,)
ax1.errorbar(T1[:p], PW3[:p], yerr = PW_error3[:p],fmt='o', ecolor='grey', label='dry method', markerfacecolor='none', color= 'maroon', capsize=5, elinewidth=2,)

ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

# Choose specific Celsius values from 80°C to -20°C, every 10 degrees
celsius_values = np.arange(80, 16, -10)
# Convert Celsius to the corresponding 1/kbT values
corresponding_1_over_kbT = (celsius_values + 273) 

ax2.set_xticks(corresponding_1_over_kbT)
ax2.set_xticklabels([f'{int(t)}°C' for t in celsius_values])
ax2.set_xlabel('Temperature (°C)')


# Set labels and title
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Peak width (Hz)")
ax1.set_title('$PEO_{18}:LiTFSI$  Peak width at narrow regime')
ax1.set_ylim(49,350)
ax1.legend(loc='upper right')  # Renaming legend and changing location
# Legend and function text adjustments
ax1.text(345, 240, '$Li^7$ \n 116.6 MHz', fontsize=12, color='black', ha='center')


#calculate_and_display_slope(T1, PW1, 15, 20,0.4,0)
#calculate_and_display_slope(T1, PW1, 21, 27,0,0.6)
#calculate_and_display_slope(T2, PW2, 14, 25,0.3,0)
#calculate_and_display_slope(T3, PW3, 14, 25,0.3,0)


plt.show()




