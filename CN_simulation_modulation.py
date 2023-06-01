# -*- coding: utf-8 -*-
"""
Created on Mon May 15 21:29:53 2023

@author: gt759
"""


import numpy as np
import matplotlib.pyplot as plt
import GNSS
import ca_generator

bitspeed =  1.023e6 #actual speed 
carrier_freq = 10*bitspeed
bits_per_sinusoid = 10
Tc = 1e-3
Ratio_Coh_time_period = 8

CA_seq_single = ca_generator.PRN_NRZ(5) 
CA_unused = ca_generator.PRN_NRZ(6) 

T_dur = 1e-3 / len(CA_seq_single) #duration of chip
samples_square = 20; #sample per chip 
T_d = 8*1e-3 #coherent integration time

dt = T_dur / samples_square




Ps = -129 #dbm
Pj = -105 #dbm
N0 = -150 #dbm

Ps_linear = np.sqrt(1e-3*np.power(10, Ps/10))
Pj_linear = np.sqrt(1e-3*np.power(10, Pj/10))
N0_linear = np.sqrt(1e-3*np.power(10, N0/10))

GPS_Object = GNSS.GNSS_GPS(carrier_freq, bitspeed )
CA_seq = ca_generator.PRN_NRZ(5) 
CA_modulated_single = GPS_Object.modulate_BPSK(CA_seq, return_time = 0)


CA2 = ca_generator.PRN_NRZ(12)  
CA2_modulated_single = GPS_Object.modulate_BPSK(CA2, return_time = 0)

CA_modulated = np.zeros(0)
CA2_modulated = np.zeros(0)

for i in range(int(T_d/Tc) - 1):
    CA_modulated = np.concatenate((CA_modulated, CA_modulated_single)) #to get Td / Tc times correlation peaks ! in Time domain
    CA2_modulated = np.concatenate((CA2_modulated, CA2_modulated_single))



timespace = np.linspace(0, GPS_Object.dt*(len(CA_modulated) - 1), len(CA_modulated)) #generate time space for CWI
N0 = N0_linear*1.78*np.random.uniform(0, 1, len(CA_modulated))


CN_s = []
freq_interference = np.linspace(0, 5*1e3, 2000)
for i in range(len(freq_interference)):
    
    f_j = carrier_freq + freq_interference[i] 

    CWI = np.sin(2*np.pi*f_j*timespace + np.pi/3.3)

    #thermal_N0 = N0_linear*np.random.uniform(0, 1)
    total_sig = Pj_linear*CWI + Ps_linear*CA_modulated + N0



    auto_peak = np.abs(np.sum(total_sig*CA_modulated)) 
    cross_peak = np.abs(np.sum(total_sig*CA2_modulated))

    
    #print("performing auto correlation")
    #autocorr = np.correlate(total_sig, CA_modulated, "same")
    #print("performing cross correlation")
    #crosscorr = np.correlate(total_sig, CA2_modulated, "same")
    
    #plt.plot(autocorr)
    #plt.plot(crosscorr)


    a = auto_peak**2/cross_peak**2 
    a = 10*np.log10(a)
    #print(a)
    CN_s.append(a)
    
    if(i % 100 == 0):
        print(i)




plt.plot(freq_interference / 1e3, CN_s, label =" Td = 3ms")










