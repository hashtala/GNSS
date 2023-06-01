
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 23:51:06 2023

@author: gt759
"""


import GNSS
import ca_generator
import numpy as np
import matplotlib.pyplot as plt



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

plt.rc('font', **font)

CA1 = ca_generator.PRN_NRZ(5)
#CA2 = ca_generator.PRN_NRZ(6)
#CA3 = ca_generator.PRN_NRZ(17)



bitspeed =  1.023e6 #actual speed 
N = 20
carrier_freq = N*bitspeed
bits_per_sinusoid = 15
A = np.linspace(1, 20, 5) #amplification of noise 
gps_object = GNSS.GNSS_GPS(carrier_freq , bitspeed, bits_per_sinusoid)

GPS_signal_modulated = gps_object.modulate_BPSK(CA1, return_time = 0); #generate signal 




    
number_repeat = 25 #we will repeat same noise 5 times to evaluate C/N
#GPS_signal_modulated += CW_noise
jammer_to_noise = []
C_N0_array = []


for i in range(0, A.size):
    print("Step " + str(i) + " out of " + str(A.size))

    for j in range(number_repeat):
                          
        
        temp_gps = np.copy(GPS_signal_modulated)
        CW_noise = A[i]*gps_object.get_CW(carrier_freq +2e6, 5*bitspeed) #slowly increase noise, notice that by calculating noise again, we get different sequence (random noise)
         
        Signal_power = gps_object.get_power_discrete_sequence(temp_gps) #calcuate signal power
        Jammer_power = gps_object.get_power_discrete_sequence(CW_noise) #calculate noise power
        J_S = 20*np.log10(Jammer_power/Signal_power) #calculate J/S 
        jammer_to_noise.append(J_S) #store calcualte J_S value
        
        temp_gps += CW_noise #add noise
    
        results = gps_object.retreive_all_satellites(temp_gps) #demodulate signal and retreive CA / correlation function
        try:
                
            retreived  = results[0] #we expect to only have single CA code in signal so directly import it form results
            xcorr = retreived[2]  #do really need to look at correlation function but need to calcualte C_N0 
            N0 = (np.sum(xcorr) - np.max(xcorr) ) / xcorr.size #compute average value of noise
            C_N0 = 20*np.log10(xcorr.max() / N0) #compute in dB
            C_N0_array.append(C_N0) #append value in dB

        except:
            
            C_N0_array.append(-np.inf) #append value in dB
    
  
    
plt.plot(C_N0_array, linewidth = 2.5)
plt.xlabel("Time (ms)") #time is ms since 1 frame takes 1 ms... can be calcualted using linear space...s
plt.ylabel(" C/N0 (dB)")
plt.grid()
plt.ylim(20, 35)
#plt.legend()


