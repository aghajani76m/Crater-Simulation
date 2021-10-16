import numpy as np
import json
import math
import random
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
from Functions import *

try:
    amounts = inputs()   
    floor_type = amounts[0]  #land==1 and sea==0
    L0 = amounts[1] #m
    v0 = amounts[2] #km/s 
    rowi = amounts[3] #kg/m^3
    theta = amounts[4]  #degree
    rowt = amounts[5] #kg/m^3
    dw = amounts[5]
    E = energy_func(v0,rowi,L0) 
    If = If_func(v0,rowi,L0,theta)
    z_star = z_star_func(If,v0,rowi)
    zb = zb_func(z_star,L0,theta,rowi)
    L = L_func(z_star,zb,rowi,L0,theta)
    v = v_func(v0,rowi,L0,theta,z_star,zb,If)
    if floor_type==0:
        v =  v_seafloor_func(v,dw,rowi,L,theta)
    Dtc = Dtc_func(v,rowi,rowt,L,theta,floor_type)
    Dfr,dfr = Dfr_func(Dtc)
    E = str(E).split(".")
    x = len(E[0])
    print("the Energy of this object is {}^({}) joule.".format(int(E[0][:2]),x-2))
    if Dfr>2000 and dfr>400:
        print(".::. it would be dangerous! it has side effects for humans, animals and plants \ne.g. Thermal Radiation, Seismic Effects, Ejecta Deposit and Air Blast.")
    N = 250
    
    Q,Z = specify_Q_Z(Dtc,dfr,L0,N)
    values = values_octaves(Q,N)   
    
    total = [0]*N
    for i in range(N):
        total[i] = [0]*N
        for j in range(N):
            total[i][j] = Z[i][j] + values[i][j]
    total = np.array(total)
    X = np.linspace(-34, 34, N)
    Y = np.linspace(-34, 34, N)
    X, Y = np.meshgrid(X, Y)
    fig = plt.figure(figsize = (10,6), dpi =  100, facecolor = 'w', edgecolor = 'k')
    ax = fig.gca(projection='3d')
    
    # Creating plot 
    if floor_type==1:
        colour = cm.afmhot_r
    else:
        colour = cm.RdBu
    surf = ax.plot_surface(X, Y, total,rstride=1, cstride=1, cmap=colour,  linewidth=0, antialiased=False)
    # Adjust the limits, ticks and view angle
    ax.set_zticks([0,0.05,0.1,.15,.2])
    ax.set_zticklabels(['','','','',''])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('Diameter of crater is: "{:,}" meter!'.format(int(Dfr)))
    ax.set_zlabel('Depth of crater is: "{:,}" meter!'.format(int(dfr)))
    ax.set_title('This is the form of crater!')
    ax.view_init(30, -30)
    
    plt.show()
except:
    print("something went wrong! please try again and pay attention to guidances.")