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
    '''
    #Path Planning
    
    show_animation = True
    
    print(__file__ + " start!!")
    
    # start and goal position
    sx = 0.0  # [m]
    sy = 0.0  # [m]
    gx = 55.0  # [m]
    gy = 55.0  # [m]
    grid_size = 2.0  # [m]
    robot_radius = 1.0  # [m]
    
    # set obstacle positions
    ox, oy = [], []
    for i in range(-10, 60):
        ox.append(i)
        oy.append(-10.0)
    for i in range(-10, 60):
        ox.append(60.0)
        oy.append(i)
    for i in range(-10, 61):
        ox.append(i)
        oy.append(60.0)
    for i in range(-10, 61):
        ox.append(-10.0)
        oy.append(i)
    for i in range(-10, 40):
        ox.append(20.0)
        oy.append(i)
    for i in range(0, 40):
        ox.append(40.0)
        oy.append(60.0 - i)
    
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
    # Import Config Files
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
    Time = 50
    kT = round(Time/T)
    
    X = np.zeros((Adt.shape[0],kT))
    Xreal = np.zeros((Adt.shape[0]+Bdt.shape[1],kT))
    dX = np.zeros((Adt.shape[0]+Bdt.shape[1],1))
    
    Ref = np.zeros((4,kT))
    Y = np.zeros((4,kT))
    U = np.zeros((Bdt.shape[1],kT))
    
    t = np.arange(0.0,kT)*T
    
    # Initialise LQG Class
    LQG = lqg.LQG(Adt,Bdt,Cdt,Ddt,Kdt,Kidt,Ldt,U_e,kT,T)
    
    # Simulation
    for k in range(1,kT-1):
        
        Ref[0,int(10/T):] = 1
        Ref[1,int(15/T):] = 30*math.pi/180
        Ref[1,int(20/T):] = 0
        Ref[2,int(25/T):] = 30*math.pi/180
        Ref[2,int(30/T):] = 0
        Ref[3,int(40/T):] = 90*math.pi/180
        
        Y[:,k] = Xreal[[4,6,8,10],k]
        #Y[:,k] = X[[0,2,4,6],k]
        
        U[:,k] = LQG.calculate(U,Y,Ref,k)
        
        #dX = QD.Quad_Dynamics(Xreal[:,k],U[:,k]) # Forward Euler Integration Nonlinear Dynamics
        #Xreal[:,k+1] = Xreal[:,k] + T*dX.T
        
        #Xy = odeint(QD.Quad_Dynamics,k,Xreal[:,k],args=(U[:,k],))
        
        Xy = solve_ivp(QD.Quad_Dynamics,(0,k),Xreal[:,k],method='RK23',vectorized=True,args=(U[:,k],))
        #print(Xy.y)
        Xreal[:,k+1] = Xy.y[:,1]
        
        #Xreal[:,k+1] = rungeKutta(T,Xreal[:,k],U[:,k])
        
        #X[:,k+1] = Adt @ X[:,k] + Bdt @ U[:,k]  # Fully Linear Dynamics
        
    # Plots
    t = np.arange(0.0,kT)*T
    '''
    fig, ax = plt.subplots()
    ax.plot(t,X[0,:])
    ax.grid(True)
    ax.set_title('z')
    
    fig2, ax2 = plt.subplots()
    ax2.plot(t,X[2,:]*180/math.pi)
    ax2.grid(True)
    ax2.set_title('phi')
    
    fig3, ax3 = plt.subplots()
    ax3.plot(t,X[4,:]*180/math.pi)
    ax3.grid(True)
    ax3.set_title('theta')
    
    fig4, ax4 = plt.subplots()
    ax4.plot(t,X[6,:]*180/math.pi)
    ax4.grid(True)
    ax4.set_title('psi')
    
    fig5, ax5 = plt.subplots()
    circle = plt.Circle((5,5),radius=2,color="red")
    rectangle = plt.Rectangle((10,10),width=2,height=4,color="red")
    ax5.add_patch(circle)
    ax5.add_patch(rectangle)
    plt.axis("scaled")
    ax5.grid(True)
    ax5.set_title('Shapes')
    '''
    
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
    
    plt.show()
    
if __name__ == '__main__':
    main()
    
