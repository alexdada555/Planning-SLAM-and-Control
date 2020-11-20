import time
import control.matlab as cont
import matplotlib as mpl
import numpy as np
import scipy as sp

M = 0.799406
g = 9.81

# Dimensions of Multirotor

L = 0.269/2 # along X-axis and Y-axis Distance from left and right motor pair to center of mass

# Mass Moment of Inertia as Taken from the CAD

Ixx = 6.609E+04/(1000*10000)
Iyy = 6.610E+04/(1000*10000)
Izz = 1.159E+05/(1000*10000)

# Motor Thrust and Torque Constants (To be determined experimentally)

Ktau =  7.708e-10
Kthrust =  1.812e-07
Kthrust2 = 0.0007326
Mtau = (1/44.22)
Ku = 515.5*Mtau

# Air resistance damping coeeficients

Dxx = 0.01212
Dyy = 0.01212
Dzz = 0.0648                          

# Equilibrium Input 

#W_e = sqrt(((M*g)/(4*(Kthrust))))*ones(4,1)
W_e = ((-4*Kthrust2) + sqrt((4*Kthrust2)^2 - (4*(-M*g)*(4*Kthrust))))/(2*(4*Kthrust))*ones(4,1)
U_e = (W_e/Ku)
#U_e = [176.1,178.5,177.2,177.6]'
#W_e = U_e*Ku
W_eV = [0;0;0;0;0;0;0;0;W_e] 

# Define Discrete-Time BeagleBone Dynamics

T = 0.01; # Sample period (s)- 100Hz
ADC = 3.3/((2^12)-1) # 12-bit ADC Quantization
DAC =  3.3/((2^12)-1) # 12-bit DAC Quantization
