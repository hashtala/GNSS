 # -*- codin_g: utf-8 -*-
"""
Created on Tue Feb 21 18:56:21 2023

@author: gt759
"""

import ca_generator 
import numpy as np
import matplotlib.pyplot as plt
import sys

class GNSS_GPS(object):
    
    
    def __init__(self, carrier_freq = 1575.42e6, bitrate = 1.023e6, bits_per_sinusoid = 20):
        self.ca_size = 1023 #for GPS this will not change
        self.carrier_frequency = carrier_freq #1.5GHz 
        self.bitspeed = bitrate #1 megaherci  
        self.bits_per_sinusoid = bits_per_sinusoid; #this variable is used in many funtions and better be declared as global
        self.dt = (1/self.carrier_frequency)/self.bits_per_sinusoid #sets time step dt, one of the most important variables
        self.Q = 0
        self.I = 1
        
        
        '''
        When we set carrier freq, bitspeed, and bits per sinosoid, we have well defined time step (dt) and sequence size  (N)
        Notice that f_req = 1 / ( dt * N) where f_res is freqnecy bin spacing when FFT is performed on modulated signal
        Therefore we also need to have "built-int" frequency resolution that will be used globally in this object (just like dt)
        '''
        carrier_bitrate_ratio = int(np.round( self.carrier_frequency / self.bitspeed ))
        self.freq_res = 1 / (self.dt * (carrier_bitrate_ratio*self.bits_per_sinusoid*self.ca_size)) 
        
    def get_time_step(self):
        return self.dt
    
    def get_freq_res(self):
        return self.freq_res
    
    def sinusiodal_gaussian(self, t, alpha, t0, f):  
        return np.exp(-(alpha**2)*((t - t0)**2))*np.cos(2*np.pi*f*(t - t0)) 
    
    def get_gaussian_time_domain(self, f_center, f_3db, plot = 0, return_non_zero = 0):
        
        '''
        return non-zero flag means following:
            for adequate FFT freq reslotlion, we need many sample points, if we generate gaussian pulse, many of poitns will be 0s 
            but having signal with small part of actual signal and alot of zeros is incovnient to have as a noise. 
            so if this flag is set, only useful part of signal will be returned.
            NOTICE: dt (time interval between samples) will be still same as it is for signal so no worries...
        '''
        
        
        
        f_3db = f_3db / 2 #as argument, it is frequency from -3dB to -3dB
        alpha =  (2*np.pi*f_3db)**2 / (4*np.log(1/np.sqrt(2)));
        alpha = np.sqrt(-1*alpha)  #cleaner, easier to debug
                    
        t0 = (1/alpha)*3 #to shift signal nicely
            
        N_points = int(1/(self.dt * self.freq_res))  #this thing determines frequnecy resolution when FFT will be performed  
        
        t_linspace = np.linspace(0, N_points*self.dt - self.dt, N_points)
        sig = self.sinusiodal_gaussian(t_linspace, alpha, t0,  f_center)
        
        if(plot == 1):
            stopIndex = int( (2*t0) / self.dt )  
            plt.figure()
            font = {'family' : 'normal',
                    'weight' : 'normal',
                    'size'   : 20}

            plt.rc('font', **font)
            plt.plot(t_linspace[:stopIndex]*1e6, sig[:stopIndex], linewidth = 3)
            plt.xlabel("Time (us)")
            plt.ylabel( "Volts")
            plt.legend()
           
        if(return_non_zero):
            stopIndex = int( (2*t0) / self.dt ) 
            return sig[:stopIndex]
        else:
            
            return sig
       
    def get_power_discrete_sequence(self, discrete_seq_td):
        return np.sum(discrete_seq_td**2) / discrete_seq_td.size
        
       
    def get_gaussian_frequency_domain_normalized(self, f_center, f_3db,  plot = 0):
        
        sig = self.get_gaussian_time_domain(f_center, f_3db,  plot = 0)
        
        FFT_sig = np.fft.fft(sig)
        
        #FFT_sig = 2*np.abs(FFT_sig) / FFT_sig.size
        #FFT_sig[int(sig.size / 2):] = 0 #single sided 
        
        if(plot == 1):
            #plot normalized
            freq_linspace = np.linspace(0, FFT_sig.size*self.freq_res - self.freq_res, FFT_sig.size)
            FFT_sig_normalized = FFT_sig / FFT_sig.max()
            plt.figure()
            font = {'family' : 'normal',
                    'weight' : 'normal',
                    'size'   : 20}

            plt.rc('font', **font)
            plt.plot(freq_linspace/1e6, (FFT_sig_normalized), linewidth = 3)
            plt.xlabel("Frequency MHz")
            plt.ylabel( "Nomalized")
            plt.legend()
            
            
        return FFT_sig / FFT_sig.max()
        
        
    def get_CW(self, carrier_freq, bitspeed):
        '''
        return continuous wave of same size as modulated C/A code, wave has well defined spectrum density
        '''


        gaussian_pulse_FR = self.get_gaussian_frequency_domain_normalized(carrier_freq, bitspeed, plot = 0) #10 KHz frequency resolution

        random_phase = np.random.uniform(0, 2*np.pi, gaussian_pulse_FR.size) #generate random uniform number for phase randomization
        
        gaussian_pulse_FR = gaussian_pulse_FR*np.exp(1j*random_phase) #make random phased so that it does not evaluate to dirac delta :) 
        
        
        
        CW_time_domain = np.fft.ifft(gaussian_pulse_FR) #perform ifft after gaussian filter multiplication 
        freq_res = 1 / (self.dt * gaussian_pulse_FR.size)
        CW_dt = 1/(freq_res * CW_time_domain.size ) #check that dt is still same as object global dt (self.dt)
        
        if( abs(CW_dt - self.dt) > 0.01*self.dt): #1 percent accuracy 
           raise Exception("Wrong time steps, CCW time step is different than global time step")
         
        return CW_time_domain.real / np.max(CW_time_domain.real)
                
        
        
    def modulate_BPSK(self, ca_code, return_time = 0):
        
        '''    
        Parameters
        ----------
        ca_code : array
            C/A code sequence.
        carrier_freq : float
            carrier frequency.
        bitspeed : bloat
            chiprate for C/A core.
    
        Returns
        -------
        array
            Linear space (time)
        array
            BPSK modulated sequence.
        '''
        
        #this takes sequence of 0s and 1s and turns into -1s and 1s correspondingly
        dt = self.dt  #this just says if we have t0 period for sinusoid, dt shuold be t0/20
        t = 1/self.bitspeed #duration of each bit it seconds
    
        total_signal = np.array([])
        total_linspace = np.array([])
    
        for i in range(0, len(ca_code)):
            
            linspace = np.linspace(i*t, (i + 1)*t - dt, int(np.round(t/dt)))
            sin_temp = ca_code[i]*np.sin(2*np.pi*self.carrier_frequency*linspace)
            total_signal = np.concatenate((total_signal,sin_temp))
            total_linspace = np.concatenate((total_linspace,linspace)) 
            
        if(return_time == 1):
            return np.array(total_linspace), np.array(total_signal)
        else:
            return  np.array(total_signal)
    
    

    


    def demodulate_BPSK_I_Q(self, signal, I_Q, doppler = 0):
        '''
        Parameters
        ----------
        signal : array
            BPSK modulated CA code.
        carrier_freq : float
            carrier frequency.
        bitspeed : float
            bits per second for C/A code
        I_Q : int 
            demodulation I or Q - 0 for I and 1 for Q
        Returns : array
        -------
        Demodulated BPSK
        
        I demodulated (sine) - extremely simple mixer 
        '''
        
        if(I_Q != 1 and I_Q != 0):
            raise Exception(" I or Q is not chosen properl, use 0 for I and 1 for Q ")
        
        length = len(signal)
        
        linspace = np.linspace(0, length*self.dt, length) #generates linear space based on length of singal and dt
        #notice that dt will be same as for signal since it is controlled using carrier frequency and bits_per_sinusoid
        if(I_Q == 0):
            carrier = np.sin(2*np.pi*(self.carrier_frequency + doppler)*linspace)
        else:
            carrier = np.cos(2*np.pi*(self.carrier_frequency + doppler)*linspace)
        
        step_down_signal = carrier * signal
        del carrier
        
        bit_sequence = [];
        integration_num_bits = int((len(signal)/self.bits_per_sinusoid)*(self.bitspeed /self.carrier_frequency))
        
        integration_index = int( (1/self.bitspeed) / self.dt)
        for i in range(integration_num_bits):
            bit_signal = np.sum(step_down_signal[(i)*integration_index : (i + 1)*integration_index]) / integration_index
            if(bit_signal > 0):
                bit_sequence.append(1) 
            else:
                bit_sequence.append(-1)     
        
        
        
        return np.array(bit_sequence)
    

    
    def retrieve_CA_signal_max_correlation(self, ca_code):
        '''
        Parameters
        ----------
        ca_code_sequence : array
            CA_code_Sequence, most likely demodulated .
    
        Returns
        -------
        max_lookup_index : INT
            returns index of CA code that it matched.
    
        function takes demodulated signal (CA code) and tries to look up matching sequence. returns index of sequence taps
    
        '''
        
        CA_codes_lookup = []
        for i in range(1, 33):
            CA_codes_lookup.append(ca_generator.PRN_NRZ(i))
            
        #correlated_signals = []
        length = len(ca_code)
        max_lookup_index = 0;
        max_correlation = 0;
        max_correlation_array = np.array([])
        for i in range(len(CA_codes_lookup)):
            lookup = CA_codes_lookup[i]
            xcorr = abs(np.correlate(lookup, ca_code, "same" ) / length) #normalize and abs because I and Q, it might be flipped in peak
            #correlated_signals.append(xcorr)
            maxVal = np.max(xcorr)
            if (maxVal > max_correlation):
                max_correlation = maxVal
                max_lookup_index = i + 1 #imito ro iq waweulia 1-it dictionaryshi :)) 
                max_correlation_array = xcorr
                
        return max_lookup_index, max_correlation_array
    
    def get_array_size_MB(self, modualted_CA_array):   
        return len(modualted_CA_array) * (modualted_CA_array.itemsize) / 1e6

    def retreive_all_satellites(self, signal, threshold = 0.20):
        
        results = []
        
        
        I = self.demodulate_BPSK_I_Q(signal, 0)
        Q = self.demodulate_BPSK_I_Q(signal, 1)
        


        #the above if-else statement just checks where we have highest correlation I or Q channel and appends corresponding information if theshold peak ahs been detected
        
        CA_codes_lookup = []
        for i in range(1, 33):
            CA_codes_lookup.append(ca_generator.PRN_NRZ(i))
            
        length = len(I) #could be Q... does not matter
        for i in range(0, len(CA_codes_lookup)):
            #print("Lookup ID - " + str(i + 1))
            lookup = CA_codes_lookup[i]
            xcorr_I = abs(np.correlate(lookup, I, "same" ) / length) #normalize and abs because I and Q, it might be flipped in peak
            xcorr_Q = abs(np.correlate(lookup, Q, "same" ) / length) #normalize and abs because I and Q, it might be flipped in peak
            
            #correlated_signals.append(xcorr)           
            if(xcorr_I.max() > xcorr_Q.max()):
              if(xcorr_I.max() > threshold):
                  results.append([self.I, i + 1, xcorr_I])
            else:
                if(xcorr_I.max() > threshold):
                    results.append([self.Q, i + 1, xcorr_Q])       
        
        
        return results





    