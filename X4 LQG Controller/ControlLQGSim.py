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
import Quad_Dynamics2 as QD2

def main():
    '''
    #Path Planning
    
    show_animation = True
    
    print(__file__ + " start!!")
    
    # start and goal position
    sx = -5.0  # [m]
    sy = -5.0  # [m]
    gx = 55.0  # [m]
    gy = 0.0  # [m]
    grid_size = 2.0  # [m]
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
    #for i in range(-10, 40):
    #    ox.append(20.0)
    #    oy.append(i)
    for i in range(-10, 50): #Wall 3
        ox.append(40.0)
        oy.append(i)
    
    if show_animation:  # pragma: no cover
        plt.plot(ox, oy, "ok")
        #plt.plot(sx, sy, "og")
        #plt.plot(gx, gy, "xb")
        ax0 = plt.gca()
        circle = plt.Circle((sx,sy),radius=2,color="blue")
        circle2 = plt.Circle((gx,gy),radius=2,color="green")
        ax0.add_patch(circle)
        ax0.add_patch(circle2)
        plt.grid(True)
        plt.axis("equal")
    
    a_star = AS.AStarPlanner(ox, oy, grid_size, robot_radius)
    rx, ry = a_star.planning(sx, sy, gx, gy)
    
    if show_animation:  # pragma: no cover
        plt.plot(rx, ry, "-r")
        plt.pause(0.001)
        #plt.show()
    '''
    # Import Control Config Files
    Adt = np.loadtxt("Adt.txt", delimiter=",")
    Bdt = np.loadtxt("Bdt.txt", delimiter=",")
    Cdt = np.loadtxt("Cdt.txt", delimiter=",")
    Ddt = np.loadtxt("Ddt.txt", delimiter=",")
    Kdt = np.loadtxt("Kdt.txt", delimiter=",")
    Kidt = np.loadtxt("Kidt.txt", delimiter=",")
    Ldt = np.loadtxt("Ldt.txt", delimiter=",")
    U_e = np.loadtxt("U_e.txt", delimiter=",")
    
    # Control Variables
    T = 0.01
    Time = 100
    kT = round(Time/T)
    
    X = np.zeros((Adt.shape[0],1))
    prevX = np.zeros((Adt.shape[0],1))
    Xreal = np.zeros((Adt.shape[0]+Bdt.shape[1],1))
    Xrealprev = np.zeros((Adt.shape[0]+Bdt.shape[1],1))
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
    
    # Simulation
    for k in range(0,kT):

        if (k == int(10/T)):
            Ref[0,0] = 1
        if (k == int(15/T)):
            Ref[1,0] = 30*math.pi/180
        if (k == int(20/T)):
            Ref[1,0] = 0
        if (k == int(25/T)):
            Ref[2,0] = 30*math.pi/180
        if (k == int(30/T)):
            Ref[2,0] = 0
        if (k == int(40/T)):
            Ref[3,0] = 90*math.pi/180
        
        #Y = Xreal[[4,6,8,10]].reshape((4,1))
        Y = X[[0,2,4,6]]
        
        U = LQG.calculate(prevU,Y,Ref,True)
        
        X = Adt @ X + Bdt @ U  # Fully Linear Dynamics
        
        #dX = QD2.Quad_Dynamics(k,Xreal[:].reshape(16),U[:].reshape(4)) # Forward Euler Integration Nonlinear Dynamics
        #print(Xreal[:].shape)
        #Xreal[:] = Xreal[:] + T*dX[:].reshape((16,1))
        
        #Xy = odeint(QD.Quad_Dynamics,k,Xreal[:],args=(U[:],))
        
        #Xy = solve_ivp(QD2.Quad_Dynamics,(0,1),Xreal[:].reshape(16),method='RK45',vectorized=True,args=(U[:],))
        #Xreal = Xy.y[:,1]
        #Xrealprev = Xreal

        prevU = U
        
        Refplot[:,k] = Ref[:,0]
        Yplot[:,k] = Y[:,0]
        #Yplot[:,k] = Y[:].reshape(4)
        Uplot[:,k] = U[:,0]
        Xplot[:,k] = X[:,0]
        #Xrealplot[:,k] = Xreal[:]
        
    # Plots
    t = np.arange(0.0,kT)*T
    
    fig, ax = plt.subplots()
    ax.plot(t,Xplot[0,:])
    ax.plot(t,Refplot[0,:])
    ax.grid(True)
    ax.set_title('z')
    
    fig2, ax2 = plt.subplots()
    ax2.plot(t,Xplot[2,:]*180/math.pi)
    ax2.plot(t,Refplot[1,:]*180/math.pi)
    ax2.plot(t,Xplot[4,:]*180/math.pi)
    ax2.plot(t,Refplot[2,:]*180/math.pi)
    ax2.plot(t,Xplot[6,:]*180/math.pi)
    ax2.plot(t,Refplot[3,:]*180/math.pi)
    ax2.grid(True)
    ax2.set_title('phi')

    fig3, ax3 = plt.subplots()
    ax3.plot(t,Uplot[0,:])
    ax3.plot(t,Uplot[1,:])
    ax3.plot(t,Uplot[2,:])
    ax3.plot(t,Uplot[3,:])
    ax3.grid(True)
    ax3.set_title('U')
    
    '''
    fig5, ax5 = plt.subplots()
    circle = plt.Circle((5,5),radius=2,color="red")
    rectangle = plt.Rectangle((10,10),width=2,height=4,color="red")
    ax5.add_patch(circle)
    ax5.add_patch(rectangle)
    plt.axis("scaled")
    ax5.grid(True)
    ax5.set_title('Shapes')
    
    
    fig6, ax6 = plt.subplots()
    ax6.plot(t,Xreal[0,:])
    ax6.grid(True)
    ax6.set_title('x')

    fig7, ax7 = plt.subplots()
    ax7.plot(t,Xreal[2,:])
    ax7.grid(True)
    ax7.set_title('y')

    fig8, ax8 = plt.subplots()
    ax8.plot(t,Xreal[4,:])
    ax8.grid(True)
    ax8.set_title('z')

    fig9, ax9 = plt.subplots()
    ax9.plot(t,Xreal[6,:]*180/math.pi)
    ax9.grid(True)
    ax9.set_title('phi')

    fig10, ax10 = plt.subplots()
    ax10.plot(t,Xreal[8,:]*180/math.pi)
    ax10.grid(True)
    ax10.set_title('theta')

    fig11, ax11 = plt.subplots()
    ax11.plot(t,Xreal[10,:]*180/math.pi)
    ax11.grid(True)
    ax11.set_title('psi')
    '''
    plt.show()
    
if __name__ == '__main__':
    main()
    
