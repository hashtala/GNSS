# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 15:06:19 2023

@author: gt759
"""

import GNSS
import numpy as np
import matplotlib.pyplot as plt
import ca_generator


CA1 = ca_generator.PRN_NRZ(5) 


bitspeed =  1.023e6 #actual speed 
N = 20 #carrier speed w.r.t bitspeed



carrier_freq = N*bitspeed
f_interference = 0.98*carrier_freq
bits_per_sinusoid = 15

Ratio_Coh_time_period = 5
Ns = Ratio_Coh_time_period*(1/bitspeed)*len(CA1) #coherent integration time
CA_circular = np.zeros(0)



for i in range(int(Ns / ((1/bitspeed)*len(CA1))) ):
    CA_circular = np.concatenate((CA_circular, CA1))  #this would be an actual demodulated signalos hermanios
 
    
gps_object = GNSS.GNSS_GPS(carrier_freq,bitspeed, bits_per_sinusoid) #initialize gps object


time_axis, modulated_signal = gps_object.modulate_BPSK(CA_circular) #get modulated C/A sequence 

dt = gps_object.dt
t = 1/gps_object.bitspeed
tim_ax = np.linspace(0, (1/bitspeed)*len(CA1)*Ratio_Coh_time_period, N*bits_per_sinusoid*len(CA1)*Ratio_Coh_time_period ) 

    

CWI_narrowband = 1*np.cos(2*np.pi*f_interference*tim_ax) #Time domain representatino of Cosine noise
modulated_signal += CWI_narrowband;
#RBW_actual = 1/(gps_object.dt * modulated_signal.size) #actual bandwidth of a signal
#noise_BW =  gps_object.get_gaussian_time_domain(carrier, f_3db)
#demodulated_signal  = gps_object.demodulate_BPSK_I_Q(modulated_signal, 0)


CWI_narrowband_demodulated = gps_object.demodulate_BPSK_I_Q(modulated_signal, 0)

print("Performing Autocorrelation")

xcorr_signal = np.correlate(CA1, CA_circular, "same") / len(CA1) #normalized autocorrelation output
#xcorr_CWI = np.correlate(CA1, CWI_narrowband_demodulated, "same") / len(CA1) #normalized autocorrelation output
plt.figure()
plt.plot(xcorr_signal)
#plt.plot(xcorr_CWI)


#N_power = abs(np.sum(xcorr_CWI**2))

#print( 20*np.log10(C_power / N_power) )
'''
correlator_signal_lines = abs(np.fft.fft(xcorr_signal) / len(xcorr_signal))
freq_axis = np.linspace(0, RBW_actual*len(xcorr_signal), len(xcorr_signal))
plt.figure()
plt.plot(freq_axis/1e6, correlator_signal_lines)
plt.xlabel("Freq MHz")
'''






