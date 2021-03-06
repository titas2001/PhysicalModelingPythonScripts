#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:28:28 2022

@author: titas
"""

from scipy.io.wavfile import write
from scipy import interpolate
import scipy
from drawnow import figure, drawnow
from pylab import *
import math
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import set_matplotlib_formats
from IPython import get_ipython  # same as clear in matlab
get_ipython().magic('reset -sf')  # same as clear in matlab
get_ipython().run_line_magic('matplotlib', 'qt')
import scipy.fftpack
# import mplcursors as mpc
# set_matplotlib_formats('svg')
# import time

# %matplotlib qt  # makes plots in a separate window
plt.close('all')

# %%  Definitions
duration = 1.         # synthesised sound lenght in s
fs = 44100.             # sample rate
deltaT = 1./fs          # calculate time step size
B = 142000.             # bulk modulus of the air
rho = 1.225             # air density
c = math.sqrt(B/rho)    # calculate speed of sound
deltaX = c*deltaT       # calculate grid spacing
L = 0.4126              # lenght of a resonating tube in meters
gamma = c/L             # scaling

dur = int(math.floor(fs*duration)) # simulation duration in samples

N = int(math.floor(1/deltaX))       # lenght of a resonating tube in samples
deltaX = 1.0/N                      # calculate grid spacing
lambdaSq = c**2 * deltaT**2 / deltaX**2  # courant number squared
# lambda = gamma * deltaT /deltaX
Pm = 3637.0   # pressure at the mouthpiece
yeq = 4.09e-04  # reed equilibrium position
yc = 0.6*yeq

A = 9.856e-05   # Effective reed surface
sigma = 3000        # Damping per unit area
sigmaF = 1.          # Impact Damping

b = 0.012           # Effective reed width
k = 1766.1952       # Reed stiffness
kc = 1971200.0      # Impact stiffness constant
m = 0.03e-03        # Reed mass
alpha = 2       # Impact exponent
outPos = 0 
#%% Intialise states of the system

PsiNext = np.zeros(N)
Psi = np.zeros(N)
PsiPrev = np.zeros(N)

yNext = 0.0
y = 0.0
yPrev = 0.0

# initialise output
out = np.zeros(dur)


#%% Tube shape function 

rleft = 0.0001
S0 = pi*rleft**2       # surface area (left end)

shapeSampled = np.array([[0, pi*(0.011/2)**2], [0.03155, pi*(0.011/2)**2],
                        [0.03165, pi*(0.0127/2)**2], [0.0826, pi*(0.015/2)**2], [L, pi*(0.015/2)**2]])
xS = shapeSampled[:, 0]
yS = shapeSampled[:, 1]
minDiff = 1
for i in range(len(xS)-1):
    if minDiff > xS[i+1] - xS[i]:
        minDiff = xS[i+1] - xS[i]

newX = np.arange(xS[0], xS[-1], minDiff)

# Interpolate shape function to match the tube in array length

M = len(shapeSampled)

f = interpolate.interp1d(xS, yS)

shapeOf = f(newX)  # shape function is from now just "S"
interpolation = interpolate.interp1d(np.arange(len(shapeOf)), shapeOf, axis=0, fill_value='extrapolate')
S = interpolation(np.linspace(0, len(shapeOf), N))

plt.figure(1)
plot(S)

avgS = np.zeros(N)
i = 1
print(i)
for i in range(N-1):
    avgS[i] = (S[i+1]+S[i])/2
avgS[0] = S[0]  # just made them to work, its not correct though
avgS[N-1] = S[N-1]  # just made them to work, its not correct though

#%%
# radiation variables a1 and a2
a1 = 1/(2*(0.8216**2) * c)  # page 253 bilbao
a2 = L/(0.8216 * math.sqrt((sum(S)/len(S))*S[1]/math.pi))

Pin = 0


def draw_fig():
    plt.ylim(-2, 2)
    plt.plot(PsiNext)


# %% Reed model 2017
print("reed model 2017")

gammaArray = np.zeros(dur)
hIsZero = np.zeros(dur)
n = 0

# Loop for flow input
for n in range(dur):
    l = 1
    for l in range(N-1):
        PsiNext[l] = 2*(1-lambdaSq)*Psi[l] - PsiPrev[l] + lambdaSq * (S[l+1]/avgS[l]) * Psi[l+1] + lambdaSq * (S[l]/avgS[l]) * Psi[l-1]

    # calculate update for reed mass spring (yNext)
    if y > yc:
        ydiff = (y - yc)**alpha
    else:
        ydiff = 0
    # ydiff = (0.5+0.5*math.tanh(10e10*(y - yc)))*((y - yc)**alpha)
    
    numerator = (4*m - 2*k*(deltaT**2))*y + (deltaT*m*sigma - 2*m + deltaT*kc*sigmaF*ydiff)*yPrev + 2*A*(Pm-Pin)*deltaT**2 - 2*(deltaT**2)*kc*ydiff
    yNext = numerator/(m * (2+deltaT*sigma)+ kc*sigmaF*deltaT*ydiff)
    # yNext = ((4 - (2*k*deltaT**2)/m)/(2+deltaT*sigma))*y + ((deltaT*sigma - 2)/(2 + deltaT*sigma))*yPrev - ((2*kc*deltaT**2)/(m*(2+deltaT*sigma)))*ydiff + ((2*deltaT**2)/(m*(2+deltaT*sigma)))*(Pm-Pin)

    # if yeq > y:
    #     h = yeq - y
    # else:
    #     h = 0
    #     hIsZero[n] = 1
        
    h = yeq - y
    Ur = A * ((yNext - yPrev)/(2*deltaT))
    Gamma = (((b*h)**2)/deltaT) * (((2*deltaX) / S[0])*Ur + 4*PsiNext[1] - PsiNext[2] - 4*Psi[0] + PsiPrev[0]) - ((2*(b*h)**2)/rho) * Pm
    Lambda = ((b*h)**2 * deltaX)/(deltaT * S[0])
    gammaArray[n] = Gamma
    if Gamma > 0:
        Uf = Lambda - math.sqrt(Lambda**2 + Gamma)
    else:
        Uf = (-1) * Lambda + math.sqrt(Lambda**2 - Gamma)
    # calculate update for PsiNext[0]
    PsiNext[0] = ((2*deltaX)/(3*S[0])) * (Ur+Uf) + (4/3)*PsiNext[1] - (1/3)*PsiNext[1]

    # Calculate the radiating boundary
    num = 2*(1-lambdaSq)*Psi[N-1] - PsiPrev[N-1] + (deltaX*((a1/deltaT) - a2) * (lambdaSq*S[N-1])/avgS[N-1])*PsiPrev[N-1] + 2*lambdaSq*Psi[N-2]
    den = 1. + deltaX*((a1/deltaT) + a2)*(lambdaSq*S[N-1]/avgS[N-1])
    PsiNext[N-1] = num/den
    # PsiNext[N-1] = Psi[N-1]

    # drawnow(draw_fig)  # draws PsiNext

    Pin = rho*((3*PsiNext[0] - 4*Psi[0] + PsiPrev[0])/(2*deltaT))

    out[n] = rho*((3*PsiNext[outPos] - 4*Psi[outPos] + PsiPrev[outPos])/(2*deltaT))

    

    yPrev = y
    y = yNext

    PsiPrev = Psi.copy()
    Psi = PsiNext.copy()

# %% Plot

plt.figure(2)
plt.plot(out)


yf = scipy.fftpack.fft(out)
xf = np.linspace(0.0, 1.0/(2.0*deltaT), dur//2)
fig, ax = plt.subplots()
ax.set_yscale("log", nonpositive='clip')
ax.plot(xf, 2.0/dur * np.abs(yf[:dur//2]))
# mpc.cursor()

# %% Save array as an audio file

scaled = np.int16(out/np.max(np.abs(out)) * 32767)
write('test.wav', int(fs), scaled)
