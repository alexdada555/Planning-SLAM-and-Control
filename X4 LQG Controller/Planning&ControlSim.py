import sys
import os
import time
import math

import numpy as np
import control.matlab as cont
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

import Linear_Quadratic_Gaussian as lqg
import AStar as AS
import Quad_Dynamics as QD

def main():
    
    #Path Planning
    show_animation = True
    
    print(__file__ + " start!!")
    
    # start and goal position
    sx = 0.0  # [m]
    sy = 0.0  # [m]
    gx = 55.0  # [m]
    gy = 0.0  # [m]
    grid_size = 1.0  # [m]
    robot_radius = 4.0  # [m]
    
    # set obstacle positions
    ox, oy = [], []
    for i in range(-10, 60): #Bottom
        ox.append(i)
        oy.append(-10.0)
    for i in range(-10, 60): #Right
        ox.append(60.0)
        oy.append(i)
    for i in range(-10, 61): #Top
        ox.append(i)
        oy.append(60.0)
    for i in range(-10, 61): #Left
        ox.append(-10.0)
        oy.append(i)
    for i in range(-10, 40): #Wall 1
        ox.append(10.0)
        oy.append(i)
    for i in range(0, 50): #Wall 2
        ox.append(25.0)
        oy.append(60.0 - i)
    for i in range(-10, 50): #Wall 3
        ox.append(40.0)
        oy.append(i)
            
    if show_animation:  # pragma: no cover
        plt.plot(ox, oy, "ok")
        ax0 = plt.gca()
        start = plt.Circle((sx,sy),radius=1,color="green")
        goal = plt.Circle((gx,gy),radius=1,color="blue")
        ax0.add_patch(start)
        ax0.add_patch(goal)
        plt.grid(True)
        plt.axis("equal")
        ax0.set_xlabel('X-pos Meters(m)')
        ax0.set_ylabel('Y-pos Meters(m)')
        ax0.set_title('Floor Map')
        
    # Initialise A* Path Planning Class
    a_star = AS.AStarPlanner(ox, oy, grid_size, robot_radius)
    #  A* Plan Path
    rx, ry = a_star.planning(sx, sy, gx, gy)
    
    if show_animation:  # pragma: no cover
        plt.plot(rx, ry, ".y")
        plt.pause(0.001)
        #plt.show()
    
    # Import Control Config Files
    Adt = np.loadtxt("Adt.txt", delimiter=",")
    Bdt = np.loadtxt("Bdt.txt", delimiter=",")
    Cdt = np.loadtxt("Cdt.txt", delimiter=",")
    Ddt = np.loadtxt("Ddt.txt", delimiter=",")
    Kdt = np.loadtxt("Kdt.txt", delimiter=",")
    Kidt = np.loadtxt("Kidt.txt", delimiter=",")
    Ldt = np.loadtxt("Ldt.txt", delimiter=",")
    U_e = np.loadtxt("U_e.txt", delimiter=",")
    
    # Simulation Variables
    T = 0.01
    sFactor = 1
    Time = (len(rx) - 1)/sFactor
    kT = round(Time/T)
    
    X = np.zeros((Adt.shape[0],1))
    prevX = np.zeros((Adt.shape[0],1))
    Xreal = np.zeros((Adt.shape[0]+Bdt.shape[1],1))
    dX = np.zeros((Adt.shape[0]+Bdt.shape[1],1))
    
    Ref = np.zeros((4,1))
    Y = np.zeros((4,1))
    U = np.zeros((Bdt.shape[1],1))
    prevU = np.zeros((Bdt.shape[1],1))
    
    Xplot = np.zeros((Adt.shape[0],kT))
    Xrealplot = np.zeros((Adt.shape[0]+Bdt.shape[1],kT))
    Refplot = np.zeros((4,kT))
    Yplot = np.zeros((4,kT))
    Uplot = np.zeros((Bdt.shape[1],kT))
    
    t = np.arange(0.0,kT)*T
    
    # Initialise LQG Class
    LQG = lqg.LQG(Adt,Bdt,Cdt,Ddt,Kdt,Kidt,Ldt,U_e)
    
    #print(rx[100])
    # Simulation
    for k in range(0,kT):
        
        Ref[0,0] = rx[len(rx)-int(k*sFactor/(T*10000))-1]
        Ref[1,0] = ry[len(ry)-int(k*sFactor/(T*10000))-1]
        Ref[2,0] = 1
        
        E_x = Ref[0,0] - Xplot[0,k]
        E_y = Ref[1,0] - Xplot[2,k]
        targPsi = (math.pi/2) - math.atan(E_y/E_x)
        Ref[3,0] = 0  
        
        #Y = Xreal[[4,6,8,10]].reshape((4,1))
        Y = X[[0,2,4,10]]
        
        U = LQG.calculate(prevU,Y,Ref,True,False)
        
        X = Adt @ X + Bdt @ U  # Linear Dynamics
        
        #dX = QD.Quad_Dynamics(k,Xreal[:].reshape(16),U[:].reshape(4)) # Forward Euler Integration Nonlinear Dynamics
        #print(Xreal[:].shape)
        #Xreal[:] = Xreal[:] + T*dX[:].reshape((16,1))
        
        #Xy = odeint(QD.Quad_Dynamics,k,Xreal[:],args=(U[:],))
        
        #Xy = solve_ivp(QD.Quad_Dynamics,(k-1,k),Xreal[:].reshape(16),method='RK45',vectorized=True,args=(U[:],))
        #Xreal = Xy.y[:,1]
        
        prevU = U
        Refplot[:,k] = Ref[:,0]
        Yplot[:,k] = Y[:,0]
        #Yplot[:,k] = Y[:].reshape(4)
        Uplot[:,k] = U[:,0]
        Xplot[:,k] = X[:,0]
        #Xrealplot[:,k] = Xreal[:].reshape(16)
        
    # Plots
    fig, ax = plt.subplots()
    ax.plot(t,Xplot[0,:])
    ax.plot(t,Refplot[0,:],"-g")
    ax.plot(t,Xplot[2,:])
    ax.plot(t,Refplot[1,:],"-b")
    ax.plot(t,Xplot[4,:])
    ax.plot(t,Refplot[2,:],"-r")
    ax.grid(True)
    ax.set_xlabel('Time(s)')
    ax.set_ylabel('Meters(m)')
    ax.set_title('Position')
        
    fig2, ax2 = plt.subplots()
    ax2.plot(t,Xplot[10,:]*180/math.pi)
    ax2.plot(t,Refplot[3,:]*180/math.pi)
    ax2.grid(True)
    ax2.set_xlabel('Time(s)')
    ax2.set_ylabel('Degrees(d)')
    ax2.set_title('Heading')
    
    fig3, ax3 = plt.subplots()
    ax3.plot(t,Uplot[0,:])
    ax3.plot(t,Uplot[1,:])
    ax3.plot(t,Uplot[2,:])
    ax3.plot(t,Uplot[3,:])
    ax3.grid(True)
    ax3.set_xlabel('Time(s)')
    ax3.set_ylabel('Micro Seconds(ms)')
    ax3.set_title('Inputs PWM  Signal')
    
    for k in range(0,kT):
        UAV = plt.Circle((Xplot[0,k],Xplot[2,k]),radius=1,color="red")
        ax0.add_patch(UAV)
        
    '''
    rectangle = plt.Rectangle((10,10),width=2,height=4,color="red")
    '''
    plt.show()
    
if __name__ == '__main__':
    main()
    
