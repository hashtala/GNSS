# -*- coding: utf-8 -*-
"""
Created on Tue May  9 19:48:43 2023

@author: gt759
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 15:06:19 2023

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


CA2 = ca_generator.PRN_NRZ(12) 

bitspeed =  1.023e6 #actual speed 
Tc = 1e-3
Ratio_Coh_time_period = 8
CA_seq_single = ca_generator.PRN_NRZ(5) 
T_dur = 1e-3 / len(CA_seq_single) #duration of chip
samples_square = 20; #sample per chip 
T_d = 8*1e-3 #coherent integration time

dt = T_dur / samples_square


Ps = -129 #dbm
Pj = -105 #dbm
N0 = -180 #dbm

Ps_linear = 1e-3*np.power(10, Ps/10)
Pj_linear = 1e-3*np.power(10, Pj/10)
N0_linear = 1e-3*np.power(10, N0/10)



temp = np.ones(samples_square)
SEQ_CA = np.zeros(0)

for i in range(len(CA_seq_single)): #what this does is takes bit and puts 20 samples for each bit (20 for example, in reality samples_square)
    temp = temp*CA_seq_single[i]
    SEQ_CA = np.concatenate((SEQ_CA, temp)) 

    
SEQ_CA2_single = np.zeros(0)

for i in range(len(CA2)): #what this does is takes bit and puts 20 samples for each bit (20 for example, in reality samples_square)
    temp = temp*CA2[i]
    SEQ_CA2_single = np.concatenate((SEQ_CA2_single, temp)) 
 

           
SEQ_CA_Single = SEQ_CA #lets save single sequence

for i in range(int(T_d/Tc) - 1):
    SEQ_CA = np.concatenate((SEQ_CA, SEQ_CA_Single)) #to get Td / Tc times correlation peaks ! in Time domain
    
f_j = 2000


timespace = np.linspace(0, dt*(len(SEQ_CA) - 1), len(SEQ_CA)) #generate time space for CWI


SEQ_CA = Ps_linear*SEQ_CA
CWI = Pj_linear*np.exp(1j*2*np.pi*f_j*timespace)
thermal_N0 = N0_linear*np.random.uniform(0, 1)



#SEQ_CA = Ps_linear*SEQ_CA
#CWI = np.cos(2*np.pi*f_j*timespace)
#thermal_N0 = np.random.uniform(0, 1)


R_s =  CWI +  SEQ_CA # + thermal_N0 



autocorrelation = np.correlate(R_s, SEQ_CA_Single) / len(R_s)
croscorrelation = np.correlate(R_s, SEQ_CA2_single) /  len(R_s)

Ps = np.sum(autocorrelation**2)
Pn = np.sum(croscorrelation**2)


print(10*np.log10(Ps/Pn))



plt.plot(autocorrelation, label = "auto-correlation")
plt.plot(croscorrelation, label = "cross-correlation")
