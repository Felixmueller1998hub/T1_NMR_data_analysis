# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 13:26:40 2024
This code opens up multiple text files convertet from Bruker to txt file 
and deconvolutes the data into two peaks more peaks may be added
@author: felix
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import PseudoVoigtModel
from lmfit.models import ConstantModel
import sys
 
#%%

# ALL THE FUNCTIONS
# Function to plot 3D mesh or surf plot
def plot_3d(data,chemical_shift, vd_list, title):
    fig = plt.figure(figsize=(12, 6), dpi=130)
    ax = fig.add_subplot(111, projection='3d')
    
    #cut the data into the right shape
    low= int(np.round(32.68875*(400-upper_limit)))
    high= int(np.round(32.68875*(400-lower_limit)))
    cut_data = data.iloc[low:high]
    
    x = chemical_shift[low:high]
    y = vd_list
    X, Y = np.meshgrid(x, y)
    Z = np.transpose(cut_data)
    
    surf = ax.plot_surface(X, np.log(Y), Z, cmap='viridis', edgecolor='k')

    ax.set_xlabel('Chemical shift [ppm]')
    ax.set_ylabel(r'$log_{10}(\tau)$ [s]')
    ax.set_zlabel('Intensity')
    ax.set_title(title)
    #ax.set_zlim(bottom=0, top=Zlim)
    
    #adjust layout
    ax.invert_xaxis()
    ax.auto_scale_xyz(X, np.log(Y), Z)   
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)  # Add colorbar for reference
   
    plt.tight_layout()
    plt.show()
    
    #ADD MORE PEAKS IN THIS MODEL if needed
def fit_Li6_LITFSI_with_Argyrodite(measurements, plot = False):
    # Model for PEO Peak
    mod_PEO = PseudoVoigtModel(prefix='PEO_')
    pars = mod_PEO.guess(measurements['intensity'], x=measurements['chemical shift'])
    pars['amplitude'].set(value=10**8, min = 10**4, max = 10**10)
    pars['center'].set(value= PEO_center, min = PEO_center - PEO_plmn, max = PEO_center + PEO_plmn)
    pars['sigma'].set(value = (PEO_sigmamax+PEO_sigmamin)/2 , min = PEO_sigmamin, max = PEO_sigmamax) 
    pars['alpha'].set(value=1, min = 0, max = 1)
    
    # Model for Argyrodite
    mod_Argyrodite = PseudoVoigtModel(prefix = 'Argyrodite_')
    pars.update(mod_Argyrodite.make_params())
    pars['Argyrodite_amplitude'].set(value=10**6, min = 10**4, max = 10**13)
    pars['Argyrodite_center'].set(value= Arg_center, min = Arg_center - Arg_plmn, max = Arg_center + Arg_plmn)
    pars['Argyrodite_sigma'].set(value= (Arg_sigmamax+Arg_sigmamin)/2 , min = Arg_sigmamin, max = Arg_sigmamax)
    pars['Argyrodite_alpha'].set(value=1, min = 0, max = 1)
    
   # Calculate background using the custom function
    background_values = background(measurements['chemical shift'], lower_limit, upper_limit) 
   
    # Model for constant background
    mod_background = ConstantModel(prefix='background_')
    pars.update(mod_background.make_params())
    pars['background_c'].set(value=background_values.mean())
    
    # add all the diferent models of the different contributions to one fit function
    mod = mod_PEO + mod_Argyrodite + mod_background
   
    # perform the fit
    out = mod.fit(measurements['intensity'], pars, x=measurements['chemical shift'])
    
    # Accessing standard errors of parameters
    std_errors = np.array([out.params[param].stderr for param in out.params])
    
    if plot == True:
        print(out.fit_report(min_correl=0.25))
        plt.figure()
        out.plot_fit()
    return out.best_fit, out.best_values, out.fit_report(show_correl=False), std_errors
    
def mixed(x, center, A, sigma, alpha):
    lorentzian_part = (A / np.pi) * (sigma / ((x - center) ** 2 + sigma ** 2))
    gaussian_part = (A / (sigma*np.sqrt(2*np.pi)))* np.exp(-0.5 * ((x - center) / sigma) ** 2)
    return (1 - alpha)*lorentzian_part + alpha * gaussian_part

def background(x, lower_limit, upper_limit):
    # Select data points outside the specified limits
    background_data = measurements.loc[(measurements['chemical shift'] < lower_limit) | (measurements['chemical shift'] > upper_limit)]
    # Calculate the average intensity of the selected data
    avg_intensity = background_data['intensity'].mean()
    # Create an array with the same length as x, filled with the average intensity
    background_values = np.full_like(x, avg_intensity)
    return background_values
 
def import_NMR(filename_measurement, start_point, end_point, lower_limit = -500, upper_limit = 500):
   
    # import the data
    measurements = pd.read_csv(filename_measurement,sep=" ", skiprows=3, header=None)
    
    # rename the column
    measurements.columns = ['intensity']
   
    # generate the chemical shift column form the start point and number of datapoints.
    measurements['chemical shift'] = np.linspace(start_point,end_point,len(measurements))
 
    # optional: only consider part of the data between the lower and upper limit
    ppm_low = measurements['chemical shift'] >= lower_limit-0.001
    ppm_high = measurements['chemical shift'] <= upper_limit+0.001
    measurements = measurements.loc[ppm_low & ppm_high]
    measurements = measurements.reset_index(drop=True)
    return measurements

def make_dataframe(T1_dataframe, chemical_shift):
    PEO_fit = pd.DataFrame()
    Arg_fit = pd.DataFrame()
    
    for i in range (19, -1, -1):
        PEO_center = T1_dataframe.at[i,'center_PEO']
        PEO_amplitude = T1_dataframe.at[i,'amplitude_PEO']
        PEO_sigma = T1_dataframe.at[i,'FWHM_PEO']/2
        PEO_alpha = T1_dataframe.at[i,'alpha_PEO']
        
        Arg_center = T1_dataframe.at[i,'center_Argyrodite']
        Arg_amplitude = T1_dataframe.at[i,'amplitude_Argyrodite']
        Arg_sigma = T1_dataframe.at[i,'FWHM_Argyrodite']/2
        Arg_alpha = T1_dataframe.at[i,'alpha_Argyrodite']
        
        Arg_fit[i] = mixed(chemical_shift, Arg_center, Arg_amplitude, Arg_sigma, Arg_alpha)
        PEO_fit[i] = mixed(chemical_shift, PEO_center, PEO_amplitude, PEO_sigma, PEO_alpha)
    
    return PEO_fit, Arg_fit


#%%
#data set
data_set = '13.1wArgyro'

#Folder title
foldertitle= '13.1wArgyro-10C(2)'

# Load vdlist from the specified path "D:\DATA\NMR\300MHz\240222-PEO-18-1_wArgyrodite\1"
vdlist_path = r"C:\Users\felix\NMRdeconvoluted_txt\18.1wArgyro\18.1wArgyro15C(2)\vdlist.txt"
vd_list = np.flip(np.loadtxt(vdlist_path))


################################### Adjust parameters as you see fit ########################################################### 

#Peak centers, start wide and then go narrow
PEO_center = -37.6
PEO_plmn = 2                #0.1 final value at lower temperatures more leeway is given 
PEO_sigmamax = 13
PEO_sigmamin = 0.3

Arg_center = -35
Arg_plmn = 2                #0.1 final value
Arg_sigmamax = 2
Arg_sigmamin = 0.2

##############################################################################################################################

#%%
# Create empty DataFrame to save rawdata for later 3D plotting
intensityData = pd.DataFrame()

# Create an empty DataFrame to save fit parameters
T1data = pd.DataFrame(columns=['FWHM_PEO', 'FWHM_Argyrodite', 'Height_PEO', 'Height_Argyrodite', 'error_heightPEO', 'error_heightArg', 'error_fwhmPEO', 'error_fwhmArg','alpha_PEO', 'alpha_Arg',  'amplitude_PEO', 'amplitude_Argyrodite'])

#FOR LOOP AND FLAG
# Flag to indicate whether to continue asking for user input
ask_user_input = True

for i in range(len(vd_list), 0, -1):
    ################################## Import the data ###############################################
    filename_measurement = rf"C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\{i}.txt"
    print(filename_measurement)
    # read start_point and end_point from header
    with open(filename_measurement, 'r') as file:
        start_point = float(file.readline().strip().split(": ")[1])
        end_point = float(file.readline().strip().split(": ")[1])

    measurements = import_NMR(filename_measurement, start_point, end_point, lower_limit=-400, upper_limit=400)
    intensityData[i-1] = measurements['intensity']
    # Limits
    lower_limit = -60
    upper_limit = -20
    
    # Fit the Data
    fit, fitted_params, report, std_errors = fit_Li6_LITFSI_with_Argyrodite(measurements, plot=False)
    
    # Calculate FWHM and Maximum Amplitude
    fwhm_PEO = 2*fitted_params['sigma']
    fwhm_Argyrodite = 2*fitted_params['Argyrodite_sigma']
    
    amplitude_PEO = fitted_params['amplitude']
    amplitude_Argyrodite = fitted_params['Argyrodite_amplitude']
    
    height_PEO = fitted_params['alpha']*(amplitude_PEO/(fitted_params['sigma']*np.pi))+ (1-fitted_params['alpha'])*(amplitude_PEO/(fitted_params['sigma']*np.sqrt(2*np.pi)))
    height_Argyrodite = amplitude_Argyrodite/(fitted_params['Argyrodite_sigma']*np.pi)
    
    alpha_PEO = fitted_params['alpha']
    alpha_Argyrodite = fitted_params['Argyrodite_alpha']
    
    error_heightPEO = height_PEO * np.sqrt((std_errors[0]/amplitude_PEO)**2 + (std_errors[2]/ fitted_params['sigma'])**2)
    error_heightArg = height_Argyrodite * np.sqrt((std_errors[6]/amplitude_Argyrodite)**2 + (std_errors[8]/fitted_params['Argyrodite_sigma'])**2 )
    
    error_fwhmPEO = 2*std_errors[2]
    error_fwhmArg = 2*std_errors[8]
    
    # Save parameters to DataFrame using loc
    T1data.loc[i - 1] = [fwhm_PEO, fwhm_Argyrodite, height_PEO, height_Argyrodite, error_heightPEO, error_heightArg, error_fwhmPEO, error_fwhmArg, alpha_PEO, alpha_Argyrodite, amplitude_PEO, amplitude_Argyrodite]

    print(report)

     # Plot the data and fit in one figure
    plt.figure(figsize = (10, 10/1.62), dpi= 130)
    plt.plot(measurements['chemical shift'], measurements['intensity'],label = 'data', linewidth = 0.6)
    plt.plot(measurements['chemical shift'], fit,'-',label = 'fit')
    plt.plot(measurements['chemical shift'], mixed(measurements['chemical shift'], fitted_params['center'], fitted_params['amplitude'], fitted_params['sigma'], fitted_params['alpha']),'--',label = 'PEO')
    plt.plot(measurements['chemical shift'], mixed(measurements['chemical shift'], fitted_params['Argyrodite_center'], fitted_params['Argyrodite_amplitude'], fitted_params['Argyrodite_sigma'], fitted_params['Argyrodite_alpha']),'--',label = 'Argyrodite')
    # Flip the x-axis and put limits
    plt.xlim(lower_limit, upper_limit)
    plt.title(rf'NMR spectra #{i}')
    plt.gca().invert_xaxis()
    plt.legend(frameon=False,fontsize = 12)     
    plt.xlabel('Chemical shift [ppm]',fontsize = 16)
    plt.ylabel('Intensity',fontsize = 16)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.tight_layout()
    plt.show()
    
    # Prompt the user for input only if ask_user_input flag is True
    if ask_user_input:
        user_input = input("Please enter 'y' if deconvolution is right, 'n' to stop, or 'p' to pass all (y/n/p): ")
        if user_input.lower() == 'n':
            print("Stopping the run.")
            sys.exit()
        elif user_input.lower() == 'y':
            print("Continuing with the next iteration.")
            continue
        elif user_input.lower() == 'p':
            print("Passing all remaining breaks.")
            ask_user_input = False  # Set the flag to False to skip further user input
            continue
    
#%%
# Convert the 'Amplitude_PEO' column from T1data DataFrame to a NumPy array
height_data_PEO = T1data['Height_PEO'].values
PEO_error = T1data['error_heightPEO'].values
height_data_Argyrodite = T1data['Height_Argyrodite'].values
Argyrodite = T1data['error_heightArg'].values

chemical_shift = measurements['chemical shift']

#plot data
plt.figure(figsize=(10, 6), dpi=130)  # Adjust figure size as needed
plt.plot(vd_list, height_data_PEO, 'bo', label='PEO Data')
plt.errorbar(vd_list, height_data_PEO, yerr=PEO_error)
plt.xscale('log')  # Set the x-axis to be logarithmic
plt.xlabel('time (log scale)')
plt.ylabel('value at maximum')
plt.title('Fit of Parameters - PEO Amplitude Data')
plt.legend()
plt.show()

#save T1data into a txt-file to be processed at a different time
output_file_path = rf"C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\T1data.txt"
np.savetxt(output_file_path, T1data)

#save vdlist into same folder
output_file_path = rf'C:\Users\felix\NMRdeconvoluted_txt\{data_set}\{foldertitle}\vdlist.txt'
#np.savetxt(output_file_path, vd_list)

#%%
######################### make 3D plot of raw data and make 3D plot of the deconvoluted fits #############################
Zlim = 2.1*10**7

#plotting raw data
chemical_shift = measurements['chemical shift']
plot_3d(intensityData, chemical_shift, vd_list,r'LPSC (10 wt.%)$-PEO_{18}LiTFSI$ at 15°C 3D-plot')

#make lorenzian PEO_dataframe Argyrodite_dataframe
dfPEO, dfArgyrodite = make_dataframe(T1data,chemical_shift)
plot_3d(dfPEO,chemical_shift, vd_list, 'Deconvoluted PEO peak')
plot_3d(dfArgyrodite,chemical_shift, vd_list, 'Deconvoluted Argyrodite peak')