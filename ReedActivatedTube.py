#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Thu Jun 30 14:04:05 2022

# @author: titas
# """

from IPython import get_ipython  # same as clear in matlab
get_ipython().magic('reset -sf') # same as clear in matlab
get_ipython().run_line_magic('matplotlib', 'qt')

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

#%%  Definitions

fs = 44100.
deltaT = 1./fs
B = 142000.              # bulk modulus of the air
rho = 1.225             # air density
c = math.sqrt(B/rho) 
deltaX = c*deltaT                 # calculate grid spacing
L = 0.4126              # scaling lenght
gamma = c/L             # scaling
durration = 1.           # synthesised sound lenght in s
dur = int(math.floor(fs*durration))

N = int(math.floor(1/(deltaX)))
deltaX = 1.0/N 
lambdaSq = c**2 * deltaT**2 / deltaX**2  #sanity check
# lambda = gamma * deltaT /deltaX 

# intialise states of the system

PsiNext = np.zeros(N) 
Psi = np.zeros(N) 
PsiPrev = np.zeros(N) 

yNext = 0.0
y = 0.0
yPrev = 0.0

Pb=2000.0
Pm = Pb
ym = 4e-04
yc = 0.6*ym
sigma = 3000  # Damping
sigmaF = 0*1.
A = 9.856e-05 
b = 0.013 
k = 8.66e06
kc = 1e10
m = 0.05
alpha = 3
# initialise output
out = np.zeros(dur) 

# Create mock shape function
#shapeSampled = 0.0001*np.ones(10)
rleft = 0.0001
S0 = pi*rleft**2;       # surface area (left end)

shapeSampled = np.array([[0, pi*(0.011/2)**2], [0.03155, pi*(0.011/2)**2], [0.03165, pi*(0.0127/2)**2], [0.0826, pi*(0.015/2)**2], [L, pi*(0.015/2)**2]])
xS = shapeSampled[:,0]
yS = shapeSampled[:,1]
minDiff = 1
for i in range(len(xS)-1):
    if minDiff > xS[i+1] - xS[i]:
        minDiff = xS[i+1] - xS[i]
        
newX = np.arange(xS[0],xS[-1],minDiff)

# Interpolate shape function to match the tube in array length

M = len(shapeSampled)

f = interpolate.interp1d(xS, yS)

xnew = np.arange(0, 1, deltaX)
shapeOf = f(newX) # shape function is from now just "S"
interpolation = interpolate.interp1d(np.arange(len(shapeOf)), shapeOf, axis=0, fill_value='extrapolate')
S = interpolation(np.linspace(0, len(shapeOf), N))
# S = 0.0002*ones(N)
plt.figure(1)
plot(S)
#%%

avgS = np.zeros(N)
i = 1
print(i)
for i in range(N-1):
    avgS[i] = (S[i+1]+S[i])/2
avgS[0] = S[0] # just made them to work, its not correct though
avgS[N-1] = S[N-1] # just made them to work, its not correct though

a1 = 1/(2*(0.8216**2) * c)  # page 253 bilbao
a2 = L/(0.8216 * math.sqrt(avgS[0]*S[0]/math.pi))

Pin = 0
# exciter[0] = 1 
f0 = 200 
# t = np.linspace(0.,dur/(f0/200.))
t = 0

def draw_fig():
    plt.ylim(-0.5,0.5)
    plt.plot(PsiNext)
    
#%% Reed model 2012 
print("reed model 2012")

gammaArray = np.zeros(dur)
n = 0
# Loop for flow input
for n in range(dur):
    l=1
    for l in range(N-1):    
        PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + \
            lambdaSq * (S[l+1]/avgS[l]) * Psi[l+1] + \
            lambdaSq * (S[l]/avgS[l]) * Psi[l-1]


    # calculate update for reed mass spring (yNext)
    if y>yc:
        ydiff = (y - yc)**alpha
    else:
        ydiff = 0
        
    yNext = ((4 - (2*k*deltaT**2)/m)/(2+deltaT*sigma))*y + ((deltaT*sigma - 2)/(2 + deltaT*sigma))*yPrev - ((2*kc*deltaT**2)/(m*(2+deltaT*sigma)))*ydiff + ((2*deltaT**2)/(m*(2+deltaT*sigma)))*(Pb-Pin)
    # numerator = (4*m - 2*k*(deltaT**2))*y + (deltaT*m*sigma - 2*m + deltaT*kc*sigmaF*ydiff)*yPrev + 2*(Pm-Pin)*deltaT**2 - 2*(deltaT**2)*kc*ydiff
    # yNext = numerator/(m * (2+deltaT*sigma)+ kc*sigmaF*deltaT*ydiff)
    
    
    
    # calculate update for PsiNext[0] 
    h = ym - yNext
    Ur = A * ((y - yPrev)/deltaT)
    Gamma = (((b*h)**2)/deltaT) * (((2*deltaX)/S[0])*Ur + 4*PsiNext[1] - PsiNext[2] - 4*Psi[0] + PsiPrev[0]) - ((2*(b*h)**2)/rho) * Pb
    Lambda = ((b*h)**2 * deltaX)/(deltaT * S[0])
    gammaArray[n] = Gamma
    if Gamma>0:
        Uf = Lambda - math.sqrt(Lambda**2 + Gamma)
    else :
        Uf = (-1) *  Lambda + math.sqrt(Lambda**2 - Gamma)

    PsiNext[0] = ((2*deltaX)/(3*S[0])) * (Ur+Uf) + (4/3)*PsiNext[1] - (1/3)*PsiNext[1]
    
    
    # Calculate the radiating boundary
    num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (deltaX*((a1/deltaT) - a2)*(lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]  
    den = 1 + deltaX*((a1/deltaT) + a2)*(lambdaSq*S[N-1]/avgS[N-1]) 
    PsiNext[N-1] = num/den 
    # PsiNext[N-1] = Psi[N-1]

    # drawnow(draw_fig)
    
    # out[n] = PsiNext[N-1]
    
    # output as pressure
    
    Pin = rho*((3*PsiNext[0] - 4*Psi[0] + PsiPrev[0])/(2*deltaT))
    
    out[n] = Pin
    

    
    yPrev = y
    y = yNext
    
    PsiPrev  = Psi.copy() 
    Psi = PsiNext.copy()   

#%% Plot

# plt.plot(exciter)
plt.figure(2)
plt.plot(out)

#%% Save array as an audio file
from scipy.io.wavfile import write

scaled = np.int16(out/np.max(np.abs(out)) * 32767)
write('test.wav', int(fs), scaled)
