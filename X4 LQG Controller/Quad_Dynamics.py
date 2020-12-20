import numpy as np
import math

# Initialise Outputs

dX = np.zeros((16,1))
Y = np.zeros((4,1))

# Motor Forces and Torques

F = np.zeros((4,1))
T = np.zeros((4,1))

def Quad_Dynamics(X,U):

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

    # X = [x xdot y ydot z zdot phi p theta q psi r w1 w2 w3 w4 w5 w6]

    x = X[0]
    xdot = X[1]
    y = X[2]
    ydot = X[3]
    z = X[4]
    zdot = X[5]

    phi = X[6]
    p = X[7]
    theta = X[8]
    q = X[9]
    psi = X[10]
    r = X[11]

    w1 = X[12]
    w2 = X[13]
    w3 = X[14]
    w4 = X[15]

    # Motor Dynamics: dX = [w1dot w2dot w3dot w4dot], U = Pulse Width of the pwm signal 0-1000%

    dX[12] = -(1/Mtau)*w1 + (Ku/Mtau)*U[0]
    dX[13] = -(1/Mtau)*w2 + (Ku/Mtau)*U[1]
    dX[14] = -(1/Mtau)*w3 + (Ku/Mtau)*U[2]
    dX[15] = -(1/Mtau)*w4 + (Ku/Mtau)*U[3]

    # Motor Forces and Torques

    F[0]= Kthrust*(math.pow(w1,2)) + Kthrust2*w1
    F[1]= Kthrust*(math.pow(w2,2)) + Kthrust2*w2
    F[2]= Kthrust*(math.pow(w3,2)) + Kthrust2*w3
    F[3]= Kthrust*(math.pow(w4,2)) + Kthrust2*w4

    T[0]= -Ktau*(math.pow(w1,2))
    T[1]=  Ktau*(math.pow(w2,2))
    T[2]=  Ktau*(math.pow(w3,2))
    T[3]= -Ktau*(math.pow(w4,2))

    Fn = sum(F)
    Tn = sum(T)

    # First Order Direvatives dX = [xdot ydot zdot phidot thetadot psidot]

    dX[0] = X[1]
    dX[2] = X[3]
    dX[4] = X[5]
    dX[6] = p + q*(math.sin(phi)*math.tan(theta)) + r*(math.cos(phi)*math.tan(theta))
    dX[8] = q*math.cos(phi) - r*math.sin(phi)
    dX[10] = q*(math.sin(phi)/math.cos(theta)) + r*(math.cos(phi)/math.cos(theta))

    # Second Order Direvatives: dX = [xddot yddot zddot pdot qdot rdot]

    dX[1] = Fn/M*(math.cos(phi)*math.sin(theta)*math.cos(psi)) + Fn/M*(math.sin(phi)*math.sin(psi)) - (Dxx/M)*xdot
    dX[3] = Fn/M*(math.cos(phi)*math.sin(theta)*math.sin(psi)) - Fn/M*(math.sin(phi)*math.cos(psi)) - (Dyy/M)*ydot
    dX[5] = Fn/M*(math.cos(phi)*math.cos(theta)) - g - (Dzz/M)*zdot

    dX[7] = (L/Ixx)*((F[0]+F[1]) - (F[2]+F[3])) - (((Izz-Iyy)/Ixx)*(r*q))
    dX[9] = (L/Iyy)*((F[0]+F[2]) - (F[1]+F[3])) - (((Izz-Ixx)/Iyy)*(p*r))
    dX[11] = Tn/Izz - (((Iyy-Ixx)/Izz)*(p*q))

    # Measured States

    Y[0] = z
    Y[1] = phi
    Y[2] = theta
    Y[3] = psi

    return dX