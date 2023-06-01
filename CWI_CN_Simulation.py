# -*- coding: utf-8 -*-
"""
Created on Tue May  9 21:57:53 2023

@author: gt759
"""


import GNSS
import numpy as np
import matplotlib.pyplot as plt
import ca_generator
import matplotlib

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

def log(x):
    return np.log10(x)

bitspeed =  1.023e6 #actual speed 
Tc = 1e-3
CA_seq_single = ca_generator.PRN_NRZ(5) 
CA_unused = ca_generator.PRN_NRZ(20) 

Ps = -129 #dbm
Pj = -105 #dbm
N0 = -175 #dbm


T_dur = 1e-3 / len(CA_seq_single) #duration of chip
samples_square = 25; #sample per chip 
T_d = 7*1e-3 #coherent integration time

dt = T_dur / samples_square




dt = T_dur / samples_square


temp = np.ones(samples_square)
CA_sampled = np.zeros(0)
CA_unused_sampled = np.zeros(0)
for i in range(len(CA_seq_single)):
    temp = temp*CA_seq_single[i]
    CA_sampled = np.concatenate((CA_sampled, temp))
    
    
    temp2 = temp*CA_unused[i]
    CA_unused_sampled = np.concatenate((CA_unused_sampled, temp2))
    
    
        
CA_sampled_long = np.zeros(0) #this will store sampled CA code 8 times or 10times depending on coherent integration time
#This variable will also be used as a replica

CA_unused_sampled_long = np.zeros(0)

for i in range(int(T_d / Tc)):
    CA_sampled_long = np.concatenate((CA_sampled_long, CA_sampled))
    CA_unused_sampled_long = np.concatenate((CA_unused_sampled_long, CA_unused_sampled))
    
    
time_space = np.linspace(0, dt*(len(CA_sampled_long) - 1), len(CA_sampled_long))

 


Ps_linear = 1e-3*np.power(10, Ps/10)
Pj_linear = 1e-3*np.power(10, Pj/10)
N0_linear = 1e-3*np.power(10, N0/10)
AGC = 55
AGC_linear = 1e-3*np.power(10, AGC/10)


A_sig = np.sqrt(Ps_linear) #amplitude should be power square root
A_jammer = np.sqrt(Pj_linear)
A_sig = np.sqrt(Ps_linear) #amplitude should be power square root
A_noise = np.sqrt(N0_linear)

sig = A_sig*CA_sampled_long 

N0 = A_noise*AGC_linear*np.random.uniform(0, 1, len(CA_sampled_long))

auto_cor = (T_d / Tc)*np.abs(np.sum((sig)*CA_sampled_long)) / len(CA_sampled_long)  #this is correlator output power for given integration time
noise_corr = np.abs(np.sum((N0)*CA_sampled_long)) / len(CA_sampled_long)
freq_interference = np.linspace(0, 5*1e3, 1500)
C_N = []
for i in range(len(freq_interference)):
    
    
    f_j = freq_interference[i]
    #NB_noise = np.cos(2*np.pi*f_j*time_space)

    NB_noise = np.exp(2j*np.pi*f_j*time_space)  #plus nosie de la thermalos

    
    
    
    noise = A_jammer*NB_noise # +
    
    
    
    #print("performing auto-correlation")
    #auto_correlation = abs(np.sum(total_sig*CA_sampled_long))
    #auto_correlation = np.correlate(total_sig, CA_sampled_long, "same" )
    #print("performing cross-correlation")
    #cross_correlation = np.correlate(total_sig, CA_unused_sampled_long, "same" )
    #plt.plot(auto_correlation)
    #plt.plot(cross_correlation)
    #argmax = np.argmax(auto_correlation)
    #a = 10*np.log10( auto_correlation[argmax] / abs(cross_correlation[argmax]))
     
    CWI_power_output = np.abs(np.sum(noise*CA_sampled_long)) / len(CA_sampled_long) #abs is I and Q output bascisally...
    #cross_peak = np.abs(np.sum(total_sig*CA_unused_sampled_long)) 
    
    #b =  10*np.log10(auto_peak / cross_peak)
    C_N.append(auto_cor**2 /(noise_corr**2 + (CWI_power_output**2) ))
    
    #C_N.append(CWI_power_output**2)
    #C_N.append(auto_cor**2)
   
    if(i % 100 == 0): 
        
        print(i)

plt.plot(freq_interference / 1e3, 10*log(C_N), label ="Simulation")
plt.legend()
















