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
samples_per_chip =  25
f_CA = 1.023*1e6 #CA code 

f_sampling = samples_per_chip*f_CA #setting sampling rate
dt = 1/f_sampling   #time difference between samples 
COH_INT = 12; #20 times more than C/A sequence length
T_d = COH_INT*((1/f_CA)*len(CA1))  #coherent integration time of 20ms
T_a = 1 / f_CA #chip duration

ACF_nonzero = np.zeros(len(CA1)) #this is array of non-zero auto correlation values -65/1023 at 12.5% of probability adn etc...



index1 = int(len(CA1)*(12.5/100))
index2 = int(len(CA1)*(12.5/100)) + index1 

ACF_nonzero[: index1] = -65/1023 
ACF_nonzero[index1: index2] = 63/1023 
ACF_nonzero[index2: ] = 1/1023 


N_c = int(np.round(T_d / dt)) #total sample points 
N_chip= int(np.round(T_a / dt)) #total sample points 

AFC_TD = np.zeros(2*N_c) #this will be output from correlator. # we wil not directly use it to calcuate powert ....

j = 0; #couter
k = 1; #this watches chip duration
d = 1 #this watches C/A sequence time 


    
R_tau = np.zeros(0)
AFC_temp = np.zeros(samples_per_chip)
n = 0
for i in range(0, COH_INT):   
    
    for j in range(1, 2*len(CA1) +1):  #since correlation does make it twice as long 
        rand = int(np.random.uniform(0, len(CA1) - 1))
        A = ACF_nonzero[rand]; #so since CA is finite in length, non-alighment ACF is not zero, and we randomly set amplitude from either 3 of possible valiues 
        #note that AFC_nonzero already has terms set accorind to probability of what non-alighment autocorrelation evaluates to
              
                             
        for n in range(0,samples_per_chip):
                
            if j == 1:
                AFC_temp[n] = 1 - (abs(n*dt) / T_a) #right half of triangle since it is first chip half
                
            
       
            elif j == 2*len(CA1): 
                
                AFC_temp[n] =  (abs(n*dt) / T_a) #left hand of the triangle since it is last chip half
                
            else:
                
                 #AFC_temp = np.zeros(samples_per_chip)
                             
                if(j % 2 == 1):
                    
                    AFC_temp[n] =   A*((abs(n*dt) / T_a)) #left hand of the triangle

                else:
                    AFC_temp[n] = A*(1 - (abs(n*dt) / T_a))  #left hand of the triangle
                    
        R_tau= np.concatenate((R_tau, AFC_temp))
                     
        #break
        
'''     
plt.figure()
timespace = np.linspace(0, len(R_tau)*dt/2, len(R_tau))
plt.plot(timespace*1e3, R_tau, linewidth = 1.5)
plt.xlabel("Time (ms)")
plt.ylabel("Normalized Correlator Output (Analytical)")
plt.grid()
'''



plt.figure()
R_tau_fft= abs(np.fft.hfft(R_tau))/ len(R_tau)

#need to normalize tau

R_tau_fft = R_tau_fft / np.sum(R_tau_fft)


RBW = 1 / (dt * len(R_tau_fft))
freq_axis = np.linspace(0, RBW*len(R_tau_fft), len(R_tau_fft))
T_a = (1/f_CA)

sinc_2 = T_a*(np.sinc(2*np.pi*freq_axis*T_a)**2)
plt.plot(freq_axis/1e3, R_tau_fft)
plt.plot(freq_axis/1e3, sinc_2)
plt.xlabel(" Freq KHz")


'''

freq_res =  1 / (len(R_tau)*dt/2)
freq_axis = np.linspace(0, freq_res*len(R_tau), len(R_tau))

correlator_SDF = (np.fft.fft(R_tau) / len(R_tau))
#correlator_SDF[ int(len(R_tau) / 2) : ] = 0
correlator_SDF = correlator_SDF /( np.sum(correlator_SDF)**2)
plt.figure()

analytical_sinc = T_a*(np.sinc(np.pi*freq_axis*T_a))

plt.plot(freq_axis/1e3, 10*np.log10(np.abs(correlator_SDF)))
plt.plot(freq_axis/1e3, 10*np.log10(np.abs(analytical_sinc)))

plt.xlabel("Freq (KHz)")
plt.ylabel("Normalized Correlator Output (Analytical)")
plt.grid()
'''       

