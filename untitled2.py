import numpy as np



def sqrt(x): #define square root so that you do not 
    return np.sqrt(x)

def get_impedance(width):
    
    h = 0.1
    er = 3
    pi = np.pi
    
    A = 60/sqrt(er)
    
    B = np.log(1 + (8*h/(pi*width)) * ( (16*h / (pi*width))  +      sqrt( (16*h / (pi*width))**2 + 6.27  ))) #impelemting analytical formula
    
    
    return A*B


#function that computes impedance difference
def imepdace_difference(Z, Z_target): 
    return Z_target - Z


w_current = 0.5 #initialize width to random value
target_impedance = 50; #target impedance

lr = 100e-5 #set learning rate (this is hyper parameter, play around with this)
iteration = 2000 #iteration parameter





for i in range(iteration): #enter loop 
    
    Z_current = get_impedance(w_current) #get current value of impedance
    
    diff = imepdace_difference(Z_current, target_impedance)  #calculate how far we are from target (smaller the better)
            
    
    w_current = w_current - lr*(diff) #move in the opposite direction of derivative (this is actually gradient descent)
    
    if( i % 200 == 0):   # we check information every 200 iterations  
        print("difference between target and current impedance " + str(abs(diff)) ) #print information
    
    

    
print("Obtained Z value " + str(get_impedance(w_current)) + " with width of " + str(w_current))  

    
    





