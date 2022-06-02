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
import math
from pylab import *
from drawnow import figure, drawnow
plt.close('all')



#%%  Definitions

fs = 44100
k = 1/fs
B = 142000;             # bulk modulus of the air
rho = 1.225;            # air density
c = math.sqrt(B/rho);
h = c*k;                # calculate grid spacing
L = 1;                  # scaling lenght
gamma = c/L;            # scaling
durration = 1;          # synthesised sound lenght in s
dur = fs*durration;

N = math.floor(1/h);
h = 1/N;
lambdaSq = c**2 * k**2 / h**2; #sanity check
# lambda = gamma * k /h;


# intialise states of the system

PsiNext = np.zeros(N);
Psi = np.zeros(N);
PsiPrev = np.zeros(N);

# intialise states of the system

PsiNext = np.zeros(N);
Psi = np.zeros(N);
PsiPrev = np.zeros(N);


# initialise output
out = np.zeros(dur);

rangeI = np.arange(1,N-1) #range for inner update equation "l"
S = np.ones(N);
M = len(S);
rangeS = np.arange(1,N-1); #range of Shape "S" for inner update equation "S_l"
avgS = np.zeros(N)

i=0
for i in range(N-1):
    avgS[i] = (S[i+1]+S[i-1])/2
avgS[0] = S[0] # just made them to work, its not correct though
avgS[N-1] = S[N-1] # just made them to work, its not correct though

a1 = 1/(2*(0.8216**2) * gamma); # page 253 bilbao
a2 = L/(0.8216 * math.sqrt(avgS[0]*S[0]/math.pi));

exciter = np.zeros(dur);
# exciter[0] = 1;
f0 = 200;
# t = np.linspace(0.,dur/(f0/200.))
t = 0
for t in range(int(dur/(f0/10))):
    exciter[t] = math.cos(2*math.pi*f0*(t/fs));


#%% Loop
def draw_fig():
    plt.ylim(-0.5,0.5)
    plt.plot(PsiNext)
    

n = 0
for n in range(dur):
    l=1
    for l in range(N-1):    
        PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
            lambdaSq * (S[l+1]/avgS[l]) * Psi[l+1] + \
            lambdaSq * (S[l]/avgS[l]) * Psi[l-1]
            
          
    # flow as an input
    PsiNext[0] = 2*lambdaSq*Psi[1] + 2*(1-lambdaSq)*Psi[0] - PsiPrev[0] + (2*h*lambdaSq*S[0]/avgS[0])*exciter[n];
    
    # pressure as an input
    # PsiNext[0] = (1-2*lambdaSq)*Psi[0] - (k/rho)*exciter[n] + \
    #     lambdaSq * (S[1]/avgS[0]) * Psi[1] + \
    #     lambdaSq * (S[0]/avgS[0] * Psi[0])
    
    num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (h*((a1/k) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]; 
    den = 1 + h*((a1/k) + a2)*(lambdaSq*S[N-1]/avgS[N-1]);
    PsiNext[N-1] = num/den;
        

    # drawnow(draw_fig)
    
    out[n] = PsiNext[N-1];
    PsiPrev  = Psi.copy();
    Psi = PsiNext.copy();  
    

#%% Plot

plt.plot(exciter)
plt.plot(out)

#%% Audio
from scipy.io.wavfile import write

scaled = np.int16(out/np.max(np.abs(out)) * 32767)
write('test.wav', 44100, scaled)
