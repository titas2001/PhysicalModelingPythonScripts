# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from IPython import get_ipython  # same as clear in matlab
get_ipython().magic('reset -sf') # same as clear in matlab

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('svg')


import numpy as np
import matplotlib.pyplot as plt
# %matplotlib qt  # makes plots in a separate window
import math
from pylab import *
from drawnow import figure, drawnow
plt.close('all')
import scipy
from scipy import interpolate, signal


EXCITER = 0 # set 0 for flow, 1 for pressure


#%%  Definitions

fs = 44100              # sampling rate
k = 1./fs
B = 142000.             # bulk modulus of the air
rho = 1.225             # air density
c = math.sqrt(B/rho) 
h = c*k                 # calculate grid spacing
L = 1.                  # scaling lenght
gamma = c/L             # scaling
durration = 1.          # synthesised sound lenght in s
dur = int(math.floor(fs*durration))

N = math.floor(1/h) 
h = 1/N 
lambdaSq = c**2 * k**2 / h**2  #sanity check
# lambda = gamma * k /h 

# intialise states of the system

PsiNext = np.zeros(N) 
Psi = np.zeros(N) 
PsiPrev = np.zeros(N) 


# initialise output
out = np.zeros(dur) 

# Create mock shape function
# random shape of a vocal tract, just to check shape function, later will be replaced with a shape of an instrument
shapeSampled = 0.001*np.array([34, 20, 12, 14, 16, 20, 26, 30, 34, 38, 34, 30, 26, 32]) 
# shapeSampled = np.linspace(1,1.1,20) 
# shapeSampled = np.ones(20)
 
# Interpolate shape function to match the tube in array length

M = len(shapeSampled)
x = np.linspace(0.01,0.1,M) 
# f = interpolate.interp1d(x, shapeSampled)

xnew = np.linspace(0.01,0.1,N)
# S = f(xnew) # shape function is from now just "S"
#%%

S = np.interp(xnew,x,shapeSampled)


#%%

avgS = np.zeros(N)

i=0
for i in range(N-1):
    avgS[i] = (S[i+1]+2*S[i]+S[i-1])/4
avgS[0] = S[0].copy()
avgS[N-1] = S[N-1].copy() 

a1 = 1/(2*(0.8216**2) * gamma)  # page 253 bilbao
a2 = L/(0.8216 * math.sqrt(avgS[0]*S[0]/math.pi))

exciterP = np.zeros(dur)
exciterV = np.zeros(dur)
# exciter[0] = 1 
f0 = 200 
# t = np.linspace(0.,dur/(f0/200.))





maximus = 0
#%% Flow input  

    
if EXCITER==0:
    flowVec = PsiNext.copy()
    def draw_fig():
        plt.ylim(-15,15)
        plt.plot(flowVec)
        
    print("flow")
    t = 0
    for t in range(int(dur/(f0/f0))):
        exciterV[t] =  math.cos(2*math.pi*f0*(t/fs))# input is a cos wave
        
    n = 0
    # Loop for flow input
    for n in range(dur):
        l=1
        for l in range(N-1):    
            PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
                lambdaSq * 0.5*((S[l+1]+S[l])/avgS[l]) * Psi[l+1] + \
                lambdaSq * 0.5*((S[l-1]+S[l])/avgS[l]) * Psi[l-1]
                
              
        # flow as an input
        PsiNext[0] = 2*lambdaSq*Psi[1] + 2*(1-lambdaSq)*Psi[0] - PsiPrev[0] + (2*h*lambdaSq*S[0]/avgS[0])*exciterV[n] 
        
        num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (h*((a1/k) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]  
        den = 1 + h*((a1/k) + a2)*(lambdaSq*S[N-1]/avgS[N-1]) 
        PsiNext[N-1] = num/den 
        
        # ##Plot the flow in the tube
        # i = 0
        # for i in range(N-1):
        #     flowVec[i] = (Psi[i+1].copy() - Psi[i].copy())/(h)
            
        # flowVec[N-1] = (Psi[N-1].copy() - Psi[N-2].copy())/(h)
            
        # drawnow(draw_fig)
        
        
        # out[n] = PsiNext[N-1]
        
        # output as flow
        out[n] = (Psi[N-3] + Psi[N-1])/(2*h)
        
        PsiPrev  = Psi.copy() # .copy() prevents from setting Psi = PsiPrev
        Psi = PsiNext.copy()  # .copy() prevents from setting PsiNext = Psi

 
#%% Pressure input
if EXCITER==1:
    pressureVec = PsiNext.copy()
    def draw_fig():
        plt.ylim(-3200,3200)
        plt.plot(pressureVec)
    
    print("pressure")
    t = 0
    for t in range(int(dur/(f0/f0))):
        exciterP[t] = 2000 * math.cos(2*math.pi*f0*(t/fs)) # input is a cos wave
        
    n = 0    
    # Loop for pressure input
    
    for n in range(dur):
        l=1
        for l in range(N-1):    
            PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
                lambdaSq * 0.5*((S[l+1]+S[l])/avgS[l]) * Psi[l+1] + \
                lambdaSq * 0.5*((S[l-1]+S[l])/avgS[l]) * Psi[l-1]
                
              
            
        # pressure as an input
        PsiNext[0] = PsiPrev[0].copy() - ((2.*k)/rho)*exciterP[n]
        
        num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (h*((a1/k) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]  
        den = 1 + h*((a1/k) + a2)*(lambdaSq*S[N-1]/avgS[N-1]) 
        PsiNext[N-1] = num/den 
            
        # ##Plot the pressure in the tube
        i = 0
        for i in range(N):
            pressureVec[i] = (PsiNext[i] - PsiPrev[i])/(2*k)
             
        # drawnow(draw_fig)
        
        
        # out[n] = PsiNext[N-1]
        
        # output as pressure
        out[n] = rho*(PsiNext[N-1]-PsiPrev[N-1])/(2*k)
        
        
        PsiPrev  = Psi.copy() # .copy() prevents from setting Psi = PsiPrev
        Psi = PsiNext.copy()  # .copy() prevents from setting PsiNext = Psi
    

#%% Plot

# plt.plot(exciter)
    
plt.plot(out)

#%% Save array as an audio file
from scipy.io.wavfile import write

scaled = np.int16(out/np.max(np.abs(out)) * 32767)
write('test.wav', int(fs), scaled)
