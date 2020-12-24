import numpy as np
import math

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
        self.U = np.zeros((self.Bdt.shape[1],self.kT))
        self.Xe = np.zeros((self.Cdt.shape[0],self.kT))
        self.Y = np.zeros((self.Cdt.shape[0],self.kT))
        self.e = np.zeros((self.Cdt.shape[0],self.kT))
        
        
    def calculate(self,U,Y,Ref,k):
        
        self.Y[:,k] = Y[:,k]
        self.U[:,k-1] = U[:,k-1]
        
        self.Xest[:,k] = self.Adt @ self.Xest[:,k-1] + self.Bdt @ self.U[:,k-1] # KF Linear Prediction
        self.e[:,k] = self.Y[:,k] - self.Xest[[0,2,4,6],k]
        self.Xest[:,k] = self.Xest[:,k] + self.Ldt @ self.e[:,k]
        
        self.Xe[:,k] = self.Xe[:,k-1] + (Ref[:,k] - self.Xest[[0,2,4,6],k]) # Integrator
        self.U[:,k] = self.U_e - (self.Kdt @ self.Xest[:,k]) - (self.Kidt @ self.Xe[:,k])

        return self.U[:,k]

    
