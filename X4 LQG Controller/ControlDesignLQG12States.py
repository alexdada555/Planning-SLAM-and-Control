import sys
import time
import control.matlab as cont
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math

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

# Air resistance damping coeficients

Dxx = 0.01212
Dyy = 0.01212
Dzz = 0.0648                          

# Equilibrium Input 

W_e = ((-4*Kthrust2) + math.sqrt(math.pow((4*Kthrust2),2) - (4*(-M*g)*(4*Kthrust))))/(2*(4*Kthrust))*np.array([1,1,1,1])
U_e = (W_e/Ku)

# Define Discrete-Time BeagleBone Dynamics

T = 0.01; # Sample period (s)- 100Hz
ADC = 3.3/(math.pow((2),12)-1) # 12-bit ADC Quantization
DAC =  3.3/(math.pow((2),12)-1) # 12-bit DAC Quantization

# Define Linear Continuous-Time Multirotor Dynamics: x_dot = Ax + Bu, y = Cx + Du         

# A =  12x12 matrix
A = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 2*Kthrust*W_e[0]/M, 2*Kthrust*W_e[1]/M, 2*Kthrust*W_e[2]/M, 2*Kthrust*W_e[3]/M],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, L*2*Kthrust*W_e[0]/Ixx, L*2*Kthrust*W_e[1]/Ixx, -L*2*Kthrust*W_e[2]/Ixx, -L*2*Kthrust*W_e[3]/Ixx],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, L*2*Kthrust*W_e[0]/Iyy, -L*2*Kthrust*W_e[1]/Iyy, L*2*Kthrust*W_e[2]/Iyy, -L*2*Kthrust*W_e[3]/Iyy],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -2*Ktau*W_e[0]/Izz, 2*Ktau*W_e[1]/Izz, 2*Ktau*W_e[2]/Izz, -2*Ktau*W_e[3]/Izz],
            [0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau]])
 
# B = 12x4 matrix
B = np.array([[0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [Ku/Mtau, 0, 0, 0],
                [0, Ku/Mtau, 0, 0],
                [0, 0, Ku/Mtau, 0],
                [0, 0, 0, Ku/Mtau]])

# C = 4x12 matrix
C = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]])
     
# D = 4x4 matrix
D = np.array([[0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0]])

# Discrete-Time System

sysdt = cont.c2d(cont.ss(A,B,C,D),T,'zoh')  # Generate Discrete-Time System

Adt = sysdt.A
Bdt = sysdt.B
Cdt = sysdt.C
Ddt = sysdt.D

# System Characteristics

poles = np.linalg.eig(Adt)[0]
# System Unstable

fig, ax = plt.subplots()
xax = ([1,2,3,4,5,6,7,8,9,10,11,12])
ax.scatter(xax,poles,marker="x")
ax.grid(True)
ax.set_title('Discrete System Eigenvalues')

cntr = np.linalg.matrix_rank(cont.ctrb(Adt,Bdt))
# Fully Reachable

obs = np.linalg.matrix_rank(cont.obsv(Adt,Cdt))
# Fully Observable

# Discrete-Time Full Integral Augmaneted System 

Cr  = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]])    

r = 4;                              # number of reference inputs
n = A[0].size;                      # number of states
q = 4;                              # number of controlled outputs

Dr = np.zeros((q,4))

num = np.concatenate((Adt, np.zeros((n,r))),axis=1)
dnum = np.concatenate(((-Cr*Adt), np.eye(q)),axis=1)

Adtaug = np.concatenate((num,dnum))

Bdtaug = np.concatenate((Bdt,-Cr*Bdt))

Cdtaug = np.concatenate((C, np.zeros((r,r))),axis=1)

# Discrete-Time Full State-Feedback Control
# State feedback control design with integral control via state augmentation
# Z Phi Theta Psi are controlled outputs

Q = np.diag([1000,0,500,0,500,0,5000,0,0,0,0,0,5,20,20,0.4]) # State penalty
R = np.eye(4)*(math.pow(10,-3))  # Control penalty

RiccSol,polesCL,Kdtaug = cont.dare(Adtaug,Bdtaug,Q,R) # DT State-Feedback Controller Gains
Kdt = Kdtaug[:,0:n]    # LQR Gains
Kidt = Kdtaug[:,n:]  # Integral Gains

fig2, ax2 = plt.subplots()
ax2.scatter(polesCL.real,polesCL.imag,marker="x")
ax2.grid(True)
ax2.set_title('Discrete System Closed Loop Eigenvalues')

# Discrete-Time Kalman Filter Design x_dot = A*x + B*u + G*w, y = C*x + D*u + H*w + v

Gdt = 1e-1*np.eye(n)

Rw = np.diag([0.5,0.5,0.01,0.1,0.01,0.01,0.01,0.01,math.pow(10,-10),math.pow(10,-10),math.pow(10,-10),math.pow(10,-10)])   # Process noise covariance matrix
Rv = np.diag([500,math.pow(10,-5),math.pow(10,-5),math.pow(10,-5)])     # Measurement noise covariance matrix Note: use low gausian noice for Rv

Ldt,RiccKal,polesKal = cont.lqe(Adt,Gdt,Cdt,Rw,Rv)
 
plt.show()
