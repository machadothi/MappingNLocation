# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 19:40:52 2020

@author: thiago.cunha
"""
import numpy as np
import matplotlib.pyplot as plt
import restthru as rt

def plotinvsensor():

    #Restthru Setup
    host = 'http://localhost:4950'
    get_pose = '/motion/pose'
    get_range = '/perception/laser/1/distances'
    rt.http_init()
    

    #Plot Setup
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
    p, status = rt.http_get(host + get_pose)
   
    px = p['x']
    py = p['y']
    pth = p['th']*np.pi/180

    # leitura do laser
    range_value = [1580, 1581, 1583, 1582, 1583, 1587, 1588, 1594, 1596, 1597, 1605, 1607, 1614, 1619, 1628, 1632, 1644, 1649, 1661, 1668, 1678, 1693, 1703, 1713, 1729, 1740, 1755, 1772, 1784, 1804, 1821, 1838, 1859, 1881, 1903, 1922, 1945, 1970, 1997, 2028, 2054, 2086, 2120, 2152, 2189, 2224, 2265, 2304, 2320, 2315, 2318, 2327, 2320, 2338, 2331, 2351, 2349, 2367, 2367, 2382, 2385, 2408, 2408, 2431, 2434, 2459, 2475, 2485, 2510, 2510, 2494, 2481, 2464, 2451, 2436, 2425, 2415, 2405, 2394, 2389, 2377, 2373, 2366, 2361, 2354, 2350, 2345, 2343, 2340, 2340, 2338, 2339, 2339, 2341, 2345, 2346, 2349, 2356, 2357, 2365, 2372, 2381, 2385, 2395, 2407, 2415, 2426, 2440, 2452, 2408, 2289, 2180, 2085, 1998, 2005, 2018, 2034, 2050, 2071, 2091, 2109, 2131, 2155, 2175, 2199, 2228, 2256, 2282, 2315, 2344, 2375, 2414, 2423, 2380, 2333, 2294, 2252, 2213, 2179, 2145, 2113, 2080, 2054, 2022, 1996, 1972, 1948, 1923, 1902, 1880, 1862, 1842, 1824, 1807, 1793, 1776, 1763, 1750, 1737, 1723, 1712, 1701, 1690, 1680, 1671, 1662, 1655, 1650, 1641, 1633, 1629, 1623, 1619, 1615, 1613, 1606, 1607, 1602, 1599, 1602, 1598]

    #fi = 20*np.pi/180
    
    for r in range(len(range_value)):
        
        fi = r*np.pi/180
        mu = np.array([range_value[r], fi])
    
        # dados do laser
        precisao = 50    # mm
        passo = 1*np.pi/180     # rad
        
        #sigma[0 0] - Erro de distancia
        #sigma[1 1] - Erro de angulação
        # sigma = np.array([[2, 0],[0, 0.2]])     # altera o shape da curva (erro do sensor)
        # sigma = [0.5 0 ; 0 0.05];     # altera o shape da curva
        sigma = np.array( [[0.05, 0 ],[0, 0.005]]);     # altera o shape da curva
        invSigma = np.linalg.inv(sigma)
        Pmin = 0.4
    
        # K = 0.5     # altera o maximo da curva (0.7-0.8)
        # K = 0.02;     # altera o maximo da curva (0.7-0.8)
        K = 0.0002;     # altera o maximo da curva
    
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
                mapa[i,j] += np.log(z[i,j]/(1-z[i,j])) #odds
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
