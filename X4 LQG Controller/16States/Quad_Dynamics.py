import numpy as np
import math
import cmath

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

# Mass Moment of Inertia as Taken from the CAD
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
    
    # Motor Forces and Torques
    F[0] = Kthrust*math.pow(w1,2) + Kthrust2*w1
    F[1] = Kthrust*math.pow(w2,2) + Kthrust2*w2
    F[2] = Kthrust*math.pow(w3,2) + Kthrust2*w3
    F[3] = Kthrust*math.pow(w4,2) + Kthrust2*w4
    
    T[0] = -Ktau*math.pow(w1,2)
    T[1] =  Ktau*math.pow(w2,2)
    T[2] =  Ktau*math.pow(w3,2)
    T[3] = -Ktau*math.pow(w4,2)
    
    Fn = sum(F)
    Tn = sum(T)
    
    # First Order Direvatives dX = [xdot ydot zdot phidot thetadot psidot]
    dX[0] = xdot
    dX[2] = ydot
    dX[4] = zdot
    dX[6] = p + q*(cmath.sin(phi).real*cmath.tan(theta).real) + r*(cmath.cos(phi).real*cmath.tan(theta).real)
    dX[8] = q*cmath.cos(phi).real - r*cmath.sin(phi).real
    dX[10] = q*(cmath.sin(phi).real/cmath.cos(theta).real) + r*(cmath.cos(phi).real/cmath.cos(theta).real)
    
    # Second Order Direvatives: dX = [xddot yddot zddot pdot qdot rdot]
    dX[1] = Fn/M*(cmath.cos(phi).real*cmath.sin(theta).real*cmath.cos(psi).real) + Fn/M*(cmath.sin(phi).real*cmath.sin(psi).real) - (Dxx/M)*xdot
    dX[3] = Fn/M*(cmath.cos(phi).real*cmath.sin(theta).real*cmath.sin(psi).real) - Fn/M*(cmath.sin(phi).real*cmath.cos(psi).real) - (Dyy/M)*ydot
    dX[5] = Fn/M*(cmath.cos(phi).real*cmath.cos(theta).real) - g - (Dzz/M)*zdot
    
    dX[7] = (L/Ixx)*(F[0]+F[1] - F[2]+F[3]) - (((Izz-Iyy)/Ixx)*(r*q))
    dX[9] = (L/Iyy)*(F[0]+F[2] - F[1]+F[3]) - (((Izz-Ixx)/Iyy)*(p*r))
    dX[11] = Tn/Izz - (((Iyy-Ixx)/Izz)*(p*q))

    # Motor Dynamics: dX = [w1dot w2dot w3dot w4dot], U = Pulse Width of the pwm signal 0-1000%
    dX[12] = -(1/Mtau)*w1 + (Ku/Mtau)*U[0]
    dX[13] = -(1/Mtau)*w2 + (Ku/Mtau)*U[1]
    dX[14] = -(1/Mtau)*w3 + (Ku/Mtau)*U[2]
    dX[15] = -(1/Mtau)*w4 + (Ku/Mtau)*U[3]
    
    return dX
