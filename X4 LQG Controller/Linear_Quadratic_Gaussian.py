import numpy as np
import math
import Quad_Dynamics as QD

class LQG:

    def __init__(self,Adt,Bdt,Cdt,Ddt,Kdt,Kidt,Ldt,U_e,kT,T):

        self.Adt = Adt
        self.Bdt = Bdt
        self.Cdt = Cdt
        self.Ddt = Ddt
        self.Kdt = Kdt
        self.Kidt = Kidt
        self.Ldt = Ldt
        self.U_e = U_e
        self.kT = kT
        self.T = T
        
        self.Xest = np.zeros((self.Adt.shape[0],self.kT))
        self.Xreal = np.zeros((self.Adt.shape[0]+self.Bdt.shape[1],self.kT))
        self.dX = np.zeros((self.Adt.shape[0]+self.Bdt.shape[1],1))
        self.U = np.zeros((self.Bdt.shape[1],self.kT))
        self.X = np.zeros((self.Adt.shape[0],self.kT))
        self.Xe = np.zeros((self.Cdt.shape[0],self.kT))
        self.Y = np.zeros((self.Cdt.shape[0],self.kT))
        self.e = np.zeros((self.Cdt.shape[0],self.kT))
        
        
    def calculate(self,Ref,k):
        
        self.Xest[:,k] = self.Adt @ self.Xest[:,k-1] + self.Bdt @ self.U[:,k-1] # No KF Linear Prediction
        
        self.Y[:,k] = self.X[[0,2,4,6],k]
        self.e[:,k] = self.Y[:,k] - self.Xest[[0,2,4,6],k]
        self.Xest[:,k] = self.Xest[:,k] + self.Ldt @ self.e[:,k]
        
        self.Xe[:,k] = self.Xe[:,k-1] + (Ref[:,k] - self.Xest[[0,2,4,6],k]) # Integrator
        self.U[:,k] = self.U_e - (self.Kdt @ self.Xest[:,k]) - (self.Kidt @ self.Xe[:,k])
        
        self.X[:,k+1] = self.Adt @ self.X[:,k] + self.Bdt @ self.U[:,k]

        #self.dX = QD.Quad_Dynamics(self.Xreal[:,k],self.U[:,k]) # Forward Euler Integration Nonlinear Dynamics
        #self.Xreal[:,k+1] = self.Xreal[:,k] + T*self.dX.T

        return self.X

    
