import numpy as np
import math

# Initialise Outputs

dX = np.zeros((16,1)).reshape(16)

# Motor Forces and Torques

F = np.zeros((4,1)).reshape(4)
T = np.zeros((4,1)).reshape(4)

# Mass of the Multirotor in Kilograms as taken from the CAD

M = 0.799406
g = 9.81

# Dimensions of Multirotor

L = 0.269/2 # along X-axis and Y-axis Distance from left and right motor pair to center of mass

#  Mass Moment of Inertia as Taken from the CAD
# Inertia Matrix and Diagolalisation CAD model coordinate system rotated 90 degrees about X

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

def Quad_Dynamics(t,X,U):                    
    
    [x,xdot,y,ydot,z,zdot,phi,p,theta,q,psi,r,w1,w2,w3,w4] = X 
        
    # Motor Dynamics: dX = [w1dot w2dot w3dot w4dot], U = Pulse Width of the pwm signal 0-1000%

    #print(U.shape)
    
    dX[12] = -(1/Mtau)*w1 + (Ku/Mtau)*U[0,0]
    dX[13] = -(1/Mtau)*w2 + (Ku/Mtau)*U[1,0]
    dX[14] = -(1/Mtau)*w3 + (Ku/Mtau)*U[2,0]
    dX[15] = -(1/Mtau)*w4 + (Ku/Mtau)*U[3,0]
    
    # Motor Forces and Torques
    
    F[0]= Kthrust*(w1*w1)# + Kthrust2*w1
    F[1]= Kthrust*(w2*w2)# + Kthrust2*w2
    F[2]= Kthrust*(w3*w3)# + Kthrust2*w3
    F[3]= Kthrust*(w4*w4)# + Kthrust2*w4
    
    T[0]= -Ktau*(w1*w1)
    T[1]=  Ktau*(w2*w2)
    T[2]=  Ktau*(w3*w3)
    T[3]= -Ktau*(w4*w4)
    
    Fn = sum(F)
    Tn = sum(T)
    
    # First Order Direvatives dX = [xdot ydot zdot phidot thetadot psidot]
    
    dX[0] = X[1]
    dX[2] = X[3]
    dX[4] = X[5]
    dX[6] = p 
    dX[8] = q
    dX[10] = r
    
    # Second Order Direvatives: dX = [xddot yddot zddot pdot qdot rdot]
    
    dX[1] = Fn/M*(math.cos(phi)*math.sin(theta)*math.cos(psi)) + Fn/M*(math.sin(phi)*math.sin(psi)) - (Dxx/M)*xdot
    dX[3] = Fn/M*(math.cos(phi)*math.sin(theta)*math.sin(psi)) - Fn/M*(math.sin(phi)*math.cos(psi)) - (Dyy/M)*ydot
    dX[5] = Fn/M*(math.cos(phi)*math.cos(theta)) - g - (Dzz/M)*zdot
    
    dX[7] = (L/Ixx)*(F[0]+F[1] - F[2]+F[3]) - (((Izz-Iyy)/Ixx)*(r*q))
    dX[9] = (L/Iyy)*(F[0]+F[2] - F[1]+F[3]) - (((Izz-Ixx)/Iyy)*(p*r))
    dX[11] = Tn/Izz - (((Iyy-Ixx)/Izz)*(p*q))

    return dX
