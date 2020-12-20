import sys
import os
import time
import control.matlab as cont
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math
import Linear_Quadratic_Gaussian as lqg
import AStar as AS

def main():
    
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
    Ref = np.zeros((4,kT))

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
        
        X = LQG.calculate(Ref,k)
        
    # Plots
    t = np.arange(0.0,kT)*T

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

    plt.show()

if __name__ == '__main__':
    main()

'''
    #Path Planning
    
    show_animation = True
    
    print(__file__ + " start!!")

    # start and goal position
    sx = 10.0  # [m]
    sy = 10.0  # [m]
    gx = 50.0  # [m]
    gy = 5.0  # [m]
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
        plt.plot(ox, oy, ".k")
        plt.plot(sx, sy, "og")
        plt.plot(gx, gy, "xb")
        plt.grid(True)
        plt.axis("equal")

    a_star = AS.AStarPlanner(ox, oy, grid_size, robot_radius)
    rx, ry = a_star.planning(sx, sy, gx, gy)
    
    if show_animation:  # pragma: no cover
        plt.plot(rx, ry, "-r")
        plt.pause(0.001)
        plt.show()
'''
