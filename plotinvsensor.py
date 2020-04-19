# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 21:25:38 2020

@author: thima
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def plotinvsensor():
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
      
    # max range
    Rmax = 10
    
    # max openning
    Amax = 5.0
    
    # leitora do sensor
    d = 5
    
    mu = np.array([d, 0])
    
    sigma = np.array([[0.5  ,   0],[0   ,   0.05]])
#    sigma = np.array([[0.05 ,   0],[0   ,   0.005]])
#    sigma = np.array([[0.005  ,   0],[0   ,   0.0005]])
    
    invSigma = np.linalg.inv(sigma)
    Pmin = 0.1
    K = 0.02
#    K = 0.0002
#    K = 0.0000004
    
    x = np.linspace(0,Rmax,400)
    y = np.linspace(-Amax, Amax, 400)
    
    xx, yy = np.meshgrid(x, y)
    
    z = np.zeros((len(x), len(y)))
    
    max_value = 0
    
    for i in range(len(x)):
        
        for j in range(len(y)):
            
            if i < mu[0]:
                P = Pmin
                
            else:
                P = 0.5
                
            X = np.array([i, j])
            dummy = (X - mu)[np.newaxis]
            dummy2 = np.dot(invSigma , dummy.T)
            z[i, j] =  P + (K/(2*np.pi*sigma[0,0]*sigma[1 ,1]) + 0.5 - P)*np.exp(-.5 * np.dot((X - mu) , dummy2))
            
            if z[i , j] > max_value:
                max_value = z[i , j]
                           
    ax.plot_surface(xx, yy, z, cmap=cm.coolwarm,                       linewidth=0, antialiased=False)
   
    ax.set_xlabel('distancia')
    ax.set_ylabel('angulo')
    ax.set_zlabel('p(x|z,theta)')
    plt.show()
  

if __name__ == "__main__":
    
    plotinvsensor()
