'''
This script will use the g flux and g flux error provided by Gaia to determine the amplitude dilution of each system for which flux is provided.

AMPLITUDE DILUTION = FLUX OF VISUAL PAIRS / TOTAL FLUX OF THE SYSTEM

AMPLITUDE DILUTION ERROR = AMPLITUDE DILUTION * SQRT((VISUAL PAIR FLUX ERROR/VISUAL PAIR FLUX)**2 + (TOTAL SYSTEM FLUX ERROR/TOTAL SYSTEM FLUX)**2)

Either of these numbers can be multiplied by 100 to give the percentage value. 
 
The amplitude dilution is a measure of the amount of system flux contributed by the visual pairs in the system, and is provided in units of percentage. 

INPUTS: 
- .csv file containing the g flux and g flux error of the systems. 
- Each of the systems should have a list of g flux values and g flux error values, depending on how many stars are found in the system by gaia. 

OUTPUT: 
- Two columns added to the .csv file, one with the amplitude dilution in percentage and one with the amplitude dilution error in percentage. 
- Histogram of the amplitude dilution values for the systems.

NB: THIS CODE NEEDS TO BE RUN WHENEVER A NEW GAIA SEARCH IS PERFORMED, AS THE OUTPUT .CSV FILE FROM GAIA WILL OVERWRITE THE ONE WITH THE AMPLITUDE DILUTION AND AMPLITUDE DILUTION ERROR VALUES. 
'''

# Import modules
import pandas as pd 
import numpy as np 
import ast 
import matplotlib.pyplot as plt

# Paths and files 
path = '/users/jess/ss21b_asteroseismology/'
gaia = pd.read_csv(path + 'gaia_search_results_oscillating_stars.csv', delimiter='|') 

# Data 
g = gaia['g_flux'].map(ast.literal_eval) 
g_err = gaia['g_flux_err'].map(ast.literal_eval) 
kic = gaia['kic'] 

# Initializing lists to save output
amp_dil = []
amp_dil_err = []

# Function for amplitude dilution 
def amplitude_dilution(fluxes):
    
    total_flux = np.sum(fluxes) # Calculate the system total flux
    companions_flux = np.sum(fluxes[1:]) # Calculate the flux of the visual pairs in the system
    ad = (companions_flux / total_flux) # Calculate the amplitude dilution 
    ad_percentage = ad * 100 # Convert to a percentage 
    amp_dil.append(ad_percentage) # Save the resulting amplitude dilution

    return total_flux, companions_flux, ad

def amplitude_dilution_error(fluxes, total_flux, companion_flux, ad): 
    
    total_flux_error = []
    # Add all the flux errors in quadrature 
    for j in fluxes:
        total_flux_error.append(j**2) # Square each of the flux error values
    total_flux_error = np.sum(total_flux_error) # Sum over the flux errors
    total_flux_error = np.sqrt(total_flux_error) # Square root the result

    comp_fluxes = fluxes[1:] # Create a list of stars WITHOUT the primart source target (i.e. only visual pairs) 
    companion_flux_error = [] 
    for j in comp_fluxes: # Same process as above
        companion_flux_error.append(j**2)
    companion_flux_error = np.sum(companion_flux_error) 
    companion_flux_error = np.sqrt(companion_flux_error)  

    # Compute the amplitude dilution error 
    ad_error = ad * np.sqrt((companion_flux_error / companion_flux)**2 + (total_flux_error / total_flux)**2) * 100
    amp_dil_err.append(ad_error) # Save the resulting amplitude dilution error

# Initializing list to catch stars with no flux
no_flux = []
no_flux_kic = []

# Looping through each of the systems
for i in range(len(g)):
    if len(g[i]) == 1: # If there is only one star in the system (amplitude dilution will be 0%) 
        amp_dil.append(0) 
        amp_dil_err.append(0)
        no_flux.append(0) 
    else: # If at least one of the stars identified does not have a flux value
        if -999 in g[i]:
            no_flux.append(1)
            np_flux_kic.append(kic[i]) 
            
            g[i].remove(-999) # Remove any stars with no flux value
            # Calculate the amplitude dilution 
            total_flux, companion_flux, ad = amplitude_dilution(g[i])             
            
            g_err[i].remove(-999) # Remove any stars with no flux error value 
            # Calculate the amplitude dilution error
            amplitude_dilution_error(g_err[i], total_flux, companion_flux, ad) 
          
        else: # If all the stars have available flux and flux error values
            no_flux.append(0)
           
            # Calculate the amplitude dilution 
            total_flux, companion_flux, ad = amplitude_dilution(g[i]) 

            # Calculate the amplitude dilution error 
            amplitude_dilution_error(g_err[i], total_flux, companion_flux, ad) 

med_amp_dil = np.median(amp_dil) # Calculate the median value of amplitude dilution 

# Print information for stars with no fluxes or flux errors
print('Of the ' + str(len(no_flux)) + ' systems, ' + str(sum(no_flux)) + ' systems have at least one star with no flux value.') 
print('These systems have the following KIC IDs: ' + str(no_flux_kic))

# Determine which systems (and their KIC IDs) have an amplitude dilution >= 10% 
big_kic = []
for i in range(len(amp_dil)):
    if amp_dil[i] >= 10:
        ind = amp_dil.index(amp_dil[i]) 
        big_kic.append(kic[ind]) 

# Save the resulting amplitude dilution and amplitude dilution error into the .csv file 
gaia['amp_dil_per'] = amp_dil 
gaia['amp_dil_err_per'] = amp_dil_err
gaia.to_csv(path + 'gaia_search_results_oscillating_stars.csv', sep='|', index=False) 

# Plot a histogram of the amplitude dilutions for this sample and save it
plt.hist(amp_dil, bins=30, color='#7bccc4', edgecolor='#43a2ca') 
plt.title('Median Amplitude Dilution = ' + str(round(med_amp_dil, 2)) + '%.\nThe Following KIC IDs have Amplitude Dilution >= 10%:\n' + str(big_kic))
plt.yscale('log') 
plt.ylabel('N') 
plt.xlabel('Amplitude Dilution [%]') 
plt.savefig(path + 'amp_dil_osc_stars.png') 
