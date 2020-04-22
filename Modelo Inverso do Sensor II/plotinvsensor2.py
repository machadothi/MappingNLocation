# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 19:40:52 2020

@author: thiago.cunha
"""
import numpy as np
import matplotlib.pyplot as plt

def plotinvsensor():

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # tamanho da celula
    tamCel = 50         # mm
    # numero de celulas
    numCelX = 200         # 10 x 10 m
    numCelY = 200
    # mapa
    mapa = np.zeros((numCelX, numCelY))

    # pose do robo
    px = 3000
    py = 1500
    pth = 15*np.pi/180

    # leitura do laser
    range_value = 2000
    fi = 20*np.pi/180
    mu = np.array([range_value, fi])

    # dados do laser
    precisao = 50    # mm
    passo = 1*np.pi/180     # rad

    sigma = np.array([[2, 0],[0, 0.2]])     # altera o shape da curva
# sigma = [0.5 0 ; 0 0.05];     # altera o shape da curva
# sigma = [0.05 0 ; 0 0.005];     # altera o shape da curva
    invSigma = np.linalg.inv(sigma)
    Pmin = 0.4

    K = 0.5     # altera o maximo da curva (0.7-0.8)
#    K = 0.02;     # altera o maximo da curva (0.7-0.8)
#    K = 0.0002;     # altera o maximo da curva

    x = np.linspace(1, numCelX, numCelX)    # X
    y = np.linspace(1, numCelY, numCelY)    # Y
    xx, yy = np.meshgrid(x, y)
    max_value = 0

    z = np.zeros((len(y), len(x)))
    for i in range(len(x)):
        # posicao da celula (i,j) no referencial global
        xi = (i-1)*tamCel + tamCel/2
        for j in range(len(y)):
            yj = (j-1)*tamCel + tamCel/2
            # raio em relacao ao robo
            r = np.sqrt((xi-px)**2 + (yj-py)**2)
            # angulo em relacao ao robo
            b = np.arctan2((yj-py), (xi-px)) - pth
            # diferença em relacao a leitura do sensor
            delta = (np.array([r, b]) - mu)
            if abs(delta[1]) > 2*passo:
#           if abs(delta(2)) > 2*passo || abs(delta(1)) > 2*precisao
                continue
                 
            if r < mu[0]:
                P = Pmin
            else:
                P = 0.5

            delta[0] /= 1000  # distancia em metros
            dummy = np.dot(delta,invSigma)
            dummy2 = np.dot(dummy , delta.T)
            z[i, j] =  P + (K/(2*np.pi*sigma[0,0]*sigma[1 ,1]) + 0.5 - P)*np.exp(-.5 *dummy2)
            mapa[i,j] = mapa[i,j] + np.log(z[i,j]/(1-z[i,j]))
            if z[i,j] > max_value:
                max_value = z[i,j]

    print(max_value)
    ax.plot_surface(xx, yy, z.T, cmap='viridis', edgecolor='none')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('p(x|z,theta)')
    plt.show()


if __name__ == "__main__":
    
    plotinvsensor()
