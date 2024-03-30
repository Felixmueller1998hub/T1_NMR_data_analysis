# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:39:29 2024

@author: felix
#%%
"""
#%%
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import pandas as pd

# Define the first custom model function
def custom_model(vdlist, Amplitude, Constant, Coefficient, T1):
    return Amplitude * (Constant + Coefficient * np.exp(-vdlist / abs(T1)))

# Define the sum of the two custom models
def custom_model_sum(vdlist, Amplitude1, Constant1, Coefficient1, T1_1, Amplitude2, Constant2, Coefficient2, T1_2):
    return (
        custom_model(vdlist, Amplitude1, Constant1, Coefficient1, T1_1) +
        custom_model(vdlist, Amplitude2, Constant2, Coefficient2, T1_2)
    )

# Create a fitting function for the sum of the two models
def fit_amplitude_data_sum(vdlist, amplitude_data, plot=False):
    mod_custom_sum = Model(custom_model_sum)

    params = mod_custom_sum.make_params()

    # Parameters for custom_model
    params.add('Amplitude1', value=2*10**5, min = 0)
    params.add('Constant1', value=0.01, min = 0)
    params.add('Coefficient1', value=-1, max = 0)
    params.add('T1_1', value=20, min=1, max=50)

    # Parameters for custom_model
    params.add('Amplitude2', value=2*10**5, min = 0)
    params.add('Constant2', value=2, min = 0)
    params.add('Coefficient2', value=-2.3, max = 0)
    params.add('T1_2', value=3, min=0.1, max = 6)


    result = mod_custom_sum.fit(amplitude_data, params, vdlist=vdlist)
    
    # Print the fit report
    print(result.fit_report())

    return result.best_fit, result.best_values

# Replace {foldertitle} with the actual folder title
foldertitle = "13.1wArgyro80C"

data_set = '13.1wArgyro'

#%% import data into T1data frame and vdlist array

# File paths
t1_data_path = rf'C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\T1data.txt'
vdlist_path = rf'C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\vdlist.txt'

# Read data into DataFrame
T1data = pd.read_csv(t1_data_path, delimiter=' ', header=None, names=['FWHM_PEO', 'FWHM_Argyrodite', 'Height_PEO', 'Height_Argyrodite', 'error_heightPEO', 'error_heightArg', 'error_fwhmPEO', 'error_fwhmArg', 'center_PEO', 'center_Argyrodite', 'amplitude_PEO', 'amplitude_Argyrodite'])

# Convert the 'Amplitude_PEO' column from T1data DataFrame to a NumPy array
height_data_PEO = T1data['Height_PEO'].values
height_data_Argyrodite = T1data['Height_Argyrodite'].values

#give errors
PEO_error = T1data['error_heightPEO']
Amp_error = T1data['error_heightArg']

#what to plot
height_data = height_data_Argyrodite

# Read data into array
vd_list = np.loadtxt(vdlist_path)

# Fit the sum of the two models
fit_result_sum, best_fit_params_sum = fit_amplitude_data_sum(vd_list, height_data, PEO_error)


#%%
    # Access the best-fit parameters for each component
Amplitude_best_fit_sum1 = best_fit_params_sum['Amplitude1']
Constant_best_fit_sum1 = best_fit_params_sum['Constant1']
Coefficient_best_fit_sum1 = best_fit_params_sum['Coefficient1']
T1_best_fit_sum1 = best_fit_params_sum['T1_1']

Amplitude_best_fit_sum2 = best_fit_params_sum['Amplitude2']
Constant_best_fit_sum2 = best_fit_params_sum['Constant2']
Coefficient_best_fit_sum2 = best_fit_params_sum['Coefficient2']
T1_best_fit_sum2 = best_fit_params_sum['T1_2']

# Plotting the data and the fitted curve for the sum of the two models
plt.figure(figsize=(10, 6), dpi=130)  # Adjust figure size as needed
plt.plot(vd_list, height_data, 'bo', label='PEO Data')
plt.errorbar(vd_list, height_data_PEO, yerr=PEO_error, fmt='none', capsize=5)
plt.plot(vd_list, fit_result_sum , 'r-', label='Best Fit Curve (Sum Model)')

# Plot the individual components with T1 values
plt.plot(vd_list, custom_model(vd_list, best_fit_params_sum['Amplitude1'], best_fit_params_sum['Constant1'], best_fit_params_sum['Coefficient1'], best_fit_params_sum['T1_1']),'--', label=f'Component 1 (T1={T1_best_fit_sum1:.3f})')
plt.plot(vd_list, custom_model(vd_list, best_fit_params_sum['Amplitude2'], best_fit_params_sum['Constant2'], best_fit_params_sum['Coefficient2'], best_fit_params_sum['T1_2']),'--', label=f'Component 2 (T1={T1_best_fit_sum2:.3f})')

plt.xscale('log')

plt.xlabel('time (log scale)')
plt.ylabel('Amplitude')
plt.title(f'Fit of Parameters - 2 Component Model\nFolder Title: {foldertitle}')
plt.legend()
plt.show()