#%%
#import
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model

#%%
# defining fitting functions

def custom_model(vdlist, Amplitude, Constant, Coefficient, T1):
    return Amplitude * (Constant + Coefficient * np.exp(-vdlist / abs(T1)))

def fit_amplitude_data(vdlist, amplitude_data, error, plot=False):
    # Create a custom model
    mod_custom = Model(custom_model)

    # Set initial parameter values and constraints
    params = mod_custom.make_params()
    params['Amplitude'].set(value=2*10**5)                             #close to the max value
    params['Constant'].set(value=0.01)                                 #sets
    params['Coefficient'].set(value= -1.5)                             #sets the initial plateau length 
    params['T1'].set(value=2, min = 0.1, max = 8)                                     #sets the angle and curve between two plateaus

    # Perform the fit
    result = mod_custom.fit(amplitude_data, params, vdlist=vdlist, weight = 1/error)

     # Print the fit report
    print(result.fit_report())

    return result.best_fit, result.best_values, result.fit_report

#%% import data into T1data frame and vdlist array
# Replace {foldertitle} with the actual folder title
data_set = '18.1wArgyro'

foldertitle = "18.1wArgyro15C(2)"

# File paths
t1_data_path = rf'C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\T1data.txt'
vdlist_path = rf'C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\vdlist.txt'

# Read data into DataFrame
T1data = pd.read_csv(t1_data_path, delimiter=' ', header=None, names=['FWHM_PEO', 'FWHM_Argyrodite', 'Height_PEO', 'Height_Argyrodite', 'error_heightPEO', 'error_heightArg', 'error_fwhmPEO', 'error_fwhmArg', 'center_PEO', 'center_Argyrodite', 'amplitude_PEO', 'amplitude_Argyrodite'])

#give errors
PEO_error = T1data['error_heightPEO']
Amp_error = T1data['error_heightArg']

# Read data into array
vd_list = np.loadtxt(vdlist_path)

#%%

# Convert the 'Amplitude_PEO' column from T1data DataFrame to a NumPy array
height_data_PEO = T1data['Height_PEO'].values
height_data_Argyrodite = T1data['Height_Argyrodite'].values

# Fit the PEO amplitude data using the custom model
fit_result, best_fit_params, report = fit_amplitude_data(vd_list, height_data_PEO, PEO_error, plot=False)

# Access the best-fit parameters
Amplitude_best_fit = best_fit_params['Amplitude']
Constant_best_fit = best_fit_params['Constant']
Coefficient_best_fit = best_fit_params['Coefficient']
T1_best_fit = best_fit_params['T1']

# Plotting the data with logarithmic x-axis
plt.figure(figsize=(10, 6), dpi=130)  # Adjust figure size as needed
plt.plot(vd_list, height_data_PEO, 'bo', label='PEO Data')
plt.plot(vd_list, fit_result , 'r-', label=f'Best Fit Curve with T1 = {T1_best_fit:.2f}')
plt.errorbar(vd_list, height_data_PEO, yerr=PEO_error, fmt='none', capsize=3)
plt.xscale('log')  # Set the x-axis to be logarithmic

plt.xlabel(r'$log_{10}(\tau)$ [s]')
plt.ylabel('Amplitude')
plt.title(f'Fit of Parameters - \nFolder Title: {foldertitle}')
plt.legend()
plt.show()