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
from scipy import interpolate


EXCITER = 0 # set 0 for flow, 1 for pressure

#%%  Definitions

fs = 44100.
k = 1./fs
B = 142000.              # bulk modulus of the air
rho = 1.225             # air density
c = math.sqrt(B/rho) 
h = c*k                 # calculate grid spacing
L = 1                   # scaling lenght
gamma = c/L             # scaling
durration = 1.           # synthesised sound lenght in s
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
shapeSampled = np.ones(10)
 
# Interpolate shape function to match the tube in array length

M = len(shapeSampled)
x = np.arange(0,M,1) 
f = interpolate.interp1d(x, shapeSampled)

xnew = np.arange(0, 1, h)
S = f(xnew) # shape function is from now just "S"
#%%

avgS = np.zeros(N)

i=0
for i in range(N-1):
    avgS[i] = (S[i+1]+S[i-1])/2
avgS[0] = S[0] # just made them to work, its not correct though
avgS[N-1] = S[N-1] # just made them to work, its not correct though

a1 = 1/(2*(0.8216**2) * gamma)  # page 253 bilbao
a2 = L/(0.8216 * math.sqrt(avgS[0]*S[0]/math.pi))

exciterP = np.zeros(dur)
exciterV = np.zeros(dur)
# exciter[0] = 1 
f0 = 200 
# t = np.linspace(0.,dur/(f0/200.))
t = 0

    
for t in range(int(dur/(f0/10))):
    exciterV[t] =  math.cos(2*math.pi*f0*(t/fs))


def draw_fig():
    plt.ylim(-0.5,0.5)
    plt.plot(PsiNext)
    
#%% Flow input  
if EXCITER==0:
    print("flow")
    for t in range(int(dur/(f0/10))):
        exciterV[t] =  math.cos(2*math.pi*f0*(t/fs))# input is a cos wave
        
    n = 0
    # Loop for flow input
    for n in range(dur):
        l=1
        for l in range(N-1):    
            PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
                lambdaSq * (S[l+1]/avgS[l]) * Psi[l+1] + \
                lambdaSq * (S[l]/avgS[l]) * Psi[l-1]
                
              
        # flow as an input
        PsiNext[0] = 2*lambdaSq*Psi[1] + 2*(1-lambdaSq)*Psi[0] - PsiPrev[0] + (2*h*lambdaSq*S[0]/avgS[0])*exciterV[n] 
        
        num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (h*((a1/k) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]  
        den = 1 + h*((a1/k) + a2)*(lambdaSq*S[N-1]/avgS[N-1]) 
        PsiNext[N-1] = num/den 
            
    
        # drawnow(draw_fig)
        
        # out[n] = PsiNext[N-1]
        
        # output as flow
        out[n] = (Psi[N-3] + Psi[N-1])/(2*h)
        
        PsiPrev  = Psi.copy() 
        Psi = PsiNext.copy()   

 
#%% Pressure input
if EXCITER==1:
    print("pressure")
    for t in range(int(dur/(f0/10))):
        exciterP[t] = 2000 * math.cos(2*math.pi*f0*(t/fs)) # input is a cos wave
        
    n = 0    
    # Loop for pressure input
    
    for n in range(dur):
        l=1
        for l in range(N-1):    
            PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
                lambdaSq * (S[l+1]/avgS[l]) * Psi[l+1] + \
                lambdaSq * (S[l]/avgS[l]) * Psi[l-1]
                
              
            
        # pressure as an input
        PsiNext[0] = PsiPrev[0].copy() - ((2.*k)/rho)*exciterP[n]
        
        num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (h*((a1/k) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]  
        den = 1 + h*((a1/k) + a2)*(lambdaSq*S[N-1]/avgS[N-1]) 
        PsiNext[N-1] = num/den 
            
    
        # drawnow(draw_fig)
        
        # out[n] = PsiNext[N-1]
        
        # output as pressure
        out[n] = rho*(PsiNext[N-1]-PsiPrev[N-1])/(2*k)
        
        
        PsiPrev  = Psi.copy() 
        Psi = PsiNext.copy()   
    

#%% Plot

# plt.plot(exciter)

plt.plot(out)

#%% Save array as an audio file
from scipy.io.wavfile import write

scaled = np.int16(out/np.max(np.abs(out)) * 32767)
write('test.wav', int(fs), scaled)
