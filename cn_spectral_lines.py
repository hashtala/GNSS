import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import ca_generator

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)



def sinc(x):

    x += 0.0000001 #this is to avoid division by 0 and instead we just return value close to 1, 
    #not the best solution but fastest to implement ...
    return np.sin(x) / (x)
    

def log(x):
    return np.log10(x)

def CN_analytical(Ps, Pj, N0, Td, df, c_w):
    Ps_linear = 1e-3*np.power(10, Ps/10)
    Pj_linear = 1e-3*np.power(10, Pj/10)
    N0_linear = 1e-3*np.power(10, N0/10)



    return (Ps_linear) / ( N0_linear +  Pj_linear*Td*(c_w**2)*(sinc(np.pi*df*Td))**2)
    

def NP_analytical(Ps, Pj, N0, Td, df, c_w):
    Ps_linear = 1e-3*np.power(10, Ps/10)
    Pj_linear = 1e-3*np.power(10, Pj/10)
    N0_linear = 1e-3*np.power(10, N0/10)


    return Pj_linear*(c_w**2)*(sinc(np.pi*df*Td))**2
  
    return Pj_linear*Td*(c_w**2)*(sinc(np.pi*df*Td))**2



def SP_analytical(Ps, Pj, N0, Td, df, c_w):
    Ps_linear = 1e-3*np.power(10, Ps/10)
    Pj_linear = 1e-3*np.power(10, Pj/10)
    N0_linear = 1e-3*np.power(10, N0/10)



    
    return Ps_linear


def get_Cn_Lines(satellite_ID, return_freq = 1, singlesided = 0):



    CA_seq = ca_generator.PRN_NRZ(satellite_ID) 
    T_dur = 1e-3 / len(CA_seq)
    samples_square = 20;
    
    
    dt = T_dur / samples_square
    
    
    temp = np.ones(samples_square)
    SEQ_CA = np.zeros(0)
    
    for i in range(len(CA_seq)):
        temp = temp*CA_seq[i]
        SEQ_CA = np.concatenate((SEQ_CA, temp))
            
        
    
    
    for i in range(4):
        SEQ_CA = np.concatenate((SEQ_CA, SEQ_CA)) #to make resolution bandwith go down and see peaks 
    
    #SEQ_CA = SEQ_CA**2    
    N_samples = len(SEQ_CA)
    res_bw = 1  / ( dt*N_samples)
    
    #time_space = np.linspace(0, dt*N_samples, N_samples)
    #plt.figure()
    #plt.plot(time_space*1e3, SEQ_CA)
    #plt.xlabel("Time ms")
    freq_axis = np.linspace(0, ( N_samples - 1)*res_bw, N_samples )
    FFT_Seq = np.fft.fft(SEQ_CA) /( N_samples)

    if(singlesided == 1):
        FFT_Seq[ int(len(FFT_Seq) / 2 ):] = 0
        FFT_Seq = 2*FFT_Seq
        
    '''
    plt.figure()
    plt.plot(freq_axis/1e3, 20*np.log10(abs(FFT_Seq)  + 0.000001 )) #so that log(0) is not evaluated
    plt.xlabel(" Freq KHz")
    
    
    sinc_envelope = np.sqrt(T_dur*(sinc(np.pi*freq_axis*T_dur)**2))
    plt.plot(freq_axis/1e3,  20*np.log10(sinc_envelope))
    
    energy_sinc = np.sum(sinc_envelope**2)
    energy_spectral_lines = np.sum(abs(FFT_Seq)**2)
    
    print(energy_sinc)
    print(energy_spectral_lines)
    '''
    if return_freq == 0:
        return abs(FFT_Seq)
    return freq_axis,  abs(FFT_Seq)



'''
freq, FFT_PRN1 = get_Cn_Lines(1, singlesided= 0)
All_sattelite_FFT = np.zeros((32, len(FFT_PRN1)))
All_sattelite_FFT[0, :] = FFT_PRN1

for i in range(1, 32):
    All_sattelite_FFT[i, :] = get_Cn_Lines(i + 1, return_freq = 0, singlesided= 0)
    continue


worst_lines = np.zeros(32)
worst_line_freqs = np.arange(1, 33, 1)

for i in range(32):
    argmax =  np.argmax(All_sattelite_FFT[i,:])
    worst_line_freqs[i] = freq[argmax] / 1e3  #find what frequency it corresponds and divide by 1000 to store it as KHz
    worst_lines[i] = 20*np.log10((All_sattelite_FFT[i, argmax]))
    
    
    

plt.figure()
linspace = np.arange(1, 33, 1)
plt.bar(linspace, worst_lines)
plt.xlabel(" PRN ID")
plt.ylabel(" Worst Case Amplitude")
plt.ylim(-25, -18)


plt.figure()
linspace = np.arange(1, 33, 1)
plt.bar(linspace, worst_line_freqs)
plt.xlabel(" PRN ID")
plt.ylabel(" Worst Case Frequency KHz ")
'''



# now we would like to evaluate C/N ratio given the NBI interference

#simulation settings
Pg = -129 #dbm
#Pj = [-105] #dbm
Pj = np.linspace(-130, -85, 200)
N0 = -179 #dbm
Td = 7*1e-3 #8 ms

Threshold_CN = 30 #30 dBHz

freq_interference = np.linspace(0, 25*1e3, 5000)
 
PJ_thresholds = (Pj[-1])*np.ones(len(freq_interference)) #ok this variable will store what ...
#Pj value results in lower C/N0 than threshold, length is same as frequency 






freq_axis, FFT_PRN1 = get_Cn_Lines(9, singlesided= 1) #get frequnecy lines
RB  = freq_axis[1] - freq_axis[0]
RB_CWI = freq_interference[1] - freq_interference[0]
cw_index = int(1e3 / RB) #what this really means is htat if I want to get cw at second spectral line it will be 2*cw_index = term of array
CN = []



for i in range(len(freq_interference)):
    
    # 1. find closes frequency line (KHz multiple)
    # 2. find corresponding cw
    # 3. find corresponding df 
    # 4. get C/N...
    
    
    freq = freq_interference[i] / 1e3 #1 KHz between lines
    freq_rounded = np.round(freq) #temp now is an index with closest values
    df = abs(freq_rounded*1e3 - freq_interference[i]) 
    cw = FFT_PRN1[int(cw_index*freq_rounded)]

    for P_j in Pj: #not the best way to name variables but ok
        
        

        #cn_temp = SP_analytical(Pg, Pj, N0, Td, df, cw) 
        #cn_temp = NP_analytical(Pg, Pj, N0, Td, df, cw) 
        cn_temp = CN_analytical(Pg, P_j, N0, Td, df, cw) 
    
    
    #print(cn_temp)
    
        cn_temp = 10*log(cn_temp)
        if(cn_temp < Threshold_CN):
            PJ_thresholds[i] = P_j #- Pg #Pj /Ps in dB basically...
            break
        else:
            #print("Still higher than " + str(Threshold_CN))
            continue
        #CN.append(cn_temp) #we can still append C/Ns
    
    
    
    
CN = np.array(CN) 
#plt.figure()    
plt.plot(freq_interference / 1e3, PJ_thresholds, label = "C/N0 Threshold " + str(Threshold_CN) +"dB")
#plt.plot(freq_interference / 1e3, PJ_thresholds, label = "PRN9")
plt.ylim(Pj[0], Pj[-1] + 5)
plt.plot()
plt.grid()
plt.xlabel("CWI ΔFreq (KHz)")
plt.ylabel(" Pj (dBm)")
plt.legend()
#plt.ylim(int(Pj[0] + 15), Pj[-1] - 10 )

'''
plt.plot(freq_interference / 1e3, CN, label =" Analytical")

plt.legend()

plt.plot()
plt.grid()
plt.xlabel("CWI ΔFreq (KHz)")
plt.ylabel(" C/N0 (dB) ")
'''



