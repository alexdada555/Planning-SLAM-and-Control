import numpy as np
import math

class LQG:

    def __init__(self,Adt,Bdt,Cdt,Ddt,Kdt,Kidt,Ldt,U_e):

        self.Adt = Adt
        self.Bdt = Bdt
        self.Cdt = Cdt
        self.Ddt = Ddt
        self.Kdt = Kdt
        self.Kidt = Kidt
        self.Ldt = Ldt
        self.U_e = U_e[:].reshape(self.Bdt.shape[1],1)
        
        self.U = np.zeros((self.Bdt.shape[1],1))
        self.prevU = np.zeros((self.Bdt.shape[1],1))

        self.Xest = np.zeros((self.Adt.shape[0],1))
        self.prevXest = np.zeros((self.Adt.shape[0],1))
        
        self.Xe = np.zeros((self.Cdt.shape[0],1))
        self.prevXe = np.zeros((self.Cdt.shape[0],1))

        self.Ref = np.zeros((self.Cdt.shape[0],1))
        self.Y = np.zeros((self.Cdt.shape[0],1))
        self.e = np.zeros((self.Cdt.shape[0],1))
        
    def calculate(self,prevU,Y,Ref,Linear):
        
        self.Y = Y[:]
        self.prevU = prevU[:]
        self.Ref = Ref[:]

        if Linear == True:
            self.Xest = self.Adt @ self.prevXest + self.Bdt @ (self.prevU) # KF Linear Prediction
            self.e = self.Y - self.Xest[[0,2,4,10]]
            self.Xest = self.Xest + self.Ldt @ self.e
        
            self.Xe = self.prevXe + (self.Ref - self.Xest[[0,2,4,10]]) # Integrator
            self.U = - (self.Kdt @ self.Xest) - (self.Kidt @ self.Xe)
        else:
            self.Xest = self.Adt @ self.prevXest + self.Bdt @ (self.prevU-self.U_e) # KF Linear Prediction
            self.e = self.Y - self.Xest[[0,2,4,10]]
            self.Xest = self.Xest + self.Ldt @ self.e
        
            self.Xe = self.prevXe + (self.Ref - self.Xest[[0,2,4,10]]) # Integrator
            self.U =  np.clip(self.U_e - (self.Kdt @ self.Xest) - (self.Kidt @ self.Xe),0,800)
            
        self.prevXe = self.Xe
        self.prevXest = self.Xest
        
        return self.U[:]
