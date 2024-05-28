# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:26:32 2023
@author: felix
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from lmfit import Model
from sklearn.linear_model import LinearRegression


def bpp_300(T, K_D, Ea, tau_0, beta):
    omega_L= larmor_frequency_hz # You can adjust this value based on your specific case
    tau_c = tau_0 * np.exp(Ea / (8.617e-5 * T))
    numerator_1 = 2 * tau_c
    denominator_1 = 1 + (omega_L * tau_c) ** beta

    numerator_2 = 4 * tau_c
    denominator_2 = 1 + (2 * omega_L * tau_c) ** beta

    result = K_D * (numerator_1 / denominator_1 + numerator_2 / denominator_2)
    return result

def richards(T, K_D, tau_0, Ea):
    omega= larmor_frequency_hz
    tau = tau_0 * np.exp(Ea / (8.617e-5 * T))
    term1 = tau * np.log(1 + (1 / (omega * tau)**2))
    term2 = 4 * tau * np.log(1 + (1 / (2 * omega * tau)**2))
    R1 = K_D * (term1 + term2)
    return R1


def fit_richards(Y_data, T, error, plot=False):
    mod_2d = Model(richards)
    
    # make model parameters
    params = mod_2d.make_params()
    params['K_D'].set(value=10**7)
    params['Ea'].set(value=0.2, min=0.1, max=0.8)
    params['tau_0'].set(value= 1e-13, min=1e-20, max=1e-9)
    
    result = mod_2d.fit(Y_data, params, T=T, weight = 1/error)
    
    return result.best_fit, result.best_values, result.fit_report(show_correl=False)

def fit_bpp(Y_data, T, error, plot=False):
    mod_bpp = Model(bpp_300)
    
    # make model parameters
    params = mod_bpp.make_params()
    params['K_D'].set(value=10**7)
    params['Ea'].set(value=0.1, min=0.1, max=0.7)
    params['tau_0'].set(value= 5e-14, min=1e-19, max=1e-12)
    params['beta'].set(value = 0, min = 1, max = 2)
    
    result = mod_bpp.fit(Y_data, params, T=T, weight = 1/error)
    
    return result.best_fit, result.best_values, result.fit_report(show_correl=False)

#%%
# Data
# specify the folder
folder = r'C:\Users\felix\OneDrive\Desktop\Master Thesis'
 
# specify the filename of the datafile
filename = r'PEO_LiTSFI_13.1wArgyrodite.xlsx'
color = 'maroon'
 
# specify the sheet name
sheet_name = r'LithiumT1-P'
 
# open the folder
os.chdir(folder)
 
# import the data
data = pd.read_excel(filename,sheet_name = sheet_name)
 

T = data['Temperature (K)']
T1_lithium = data['T1 (s)']
T1_comp2 = data ['T1 comp 2 (s)']
comp1_error = data['T1 stderror']
comp2_error = data['comp2 stderror']

#%%
# Calculate 1/(T1)
Y_lithium = 1/T1_lithium
Y_comp2 = 1/T1_comp2

#transform the error to (1/y)
trans_error1 = comp1_error/T1_lithium**2
trans_error2 = comp2_error/T1_comp2**2


# Larmor frequency for 7Li at 300MHz and conversion to rad/s
omega_larmor = 116.642
larmor_frequency_hz = omega_larmor * 1e6
larmor_frequency_rad_per_s = larmor_frequency_hz * 2 * np.pi

# Fit the model with the filtered data
n=5
fit, fitted_params, report = fit_richards(Y_comp2[:n], T[:n], comp2_error[:n])
m=6
fit2, fitted_params2, report2 = fit_richards(Y_lithium[m:], T[m:], comp1_error[m:])

# Extend the temperature range for plotting the fit
T_extended = np.linspace(min(T), max(T), 500)  # Generates 500 points evenly spaced across the range

# Calculate the extended fits
extended_fit = richards(T_extended, fitted_params['K_D'], fitted_params['tau_0'], fitted_params['Ea'])
extended_fit2 = richards(T_extended, fitted_params2['K_D'], fitted_params2['tau_0'], fitted_params2['Ea'])

# Extract fitted parameters
K_D_fit = fitted_params['K_D']
#beta_fit = fitted_params['beta']
Ea_fit = fitted_params['Ea']
tau_0_fit = fitted_params['tau_0']

# Plotting code...
fig, ax1 = plt.subplots(figsize=(8, 6), dpi=180)


# Data points and fitted curve for adjusted BPP equation
ax1.errorbar(1000/(T), (Y_lithium), yerr=trans_error1, markerfacecolor='none',markeredgecolor=color, ecolor = color, fmt='o', capsize=3, label='1. component')
ax1.errorbar(1000/(T), (Y_comp2), yerr=trans_error2, markerfacecolor='black',markeredgecolor=color, ecolor = color, fmt='o', capsize=3, label='2. component')

ax1.plot(1000/T_extended, extended_fit, '--', color='black', label='2D Richards fit')
ax1.plot(1000/T_extended, extended_fit2, '--', color=color, label = '2D Richards fit')

ax1.set_yscale('log')
ax1.set_xlabel('1000/T [1/K]')
ax1.set_ylabel('log(1/T1) [$s^{-1}$]')
ax1.legend(loc='lower left')  # Renaming legend and changing location
ax1.minorticks_on()

# Adding a secondary x-axis at the top with temperature in Celsius
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

# Choose specific Celsius values from 80°C to -20°C, every 10 degrees
celsius_values = np.arange(80, -21, -10)
# Convert Celsius to the corresponding 1/kbT values
corresponding_1_over_kbT = 1000 / ((celsius_values + 273.15))

ax2.set_xticks(corresponding_1_over_kbT)
ax2.set_xticklabels([f'{int(t)}°C' for t in celsius_values])
ax2.set_xlabel('Temperature (°C)')

plt.show()
print(report, report2)
