import json
import math
import random
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys

#function for getting data, guidance and error handling
def inputs():
    inputs = []
    
    #floor_type
    inputs.append(input("please specify the type of floor.\nsea = 0 and land = 1\n").replace(" ",""))
    try: 
        inputs[0] = int(inputs[0])
        if inputs[0] != 1 and inputs[0] != 0: 
            print("the floor type 0 specify the sea and 1 specify the land.\n please just write 0 or 1 to specify the type of target.")
            sys.exit(1)
        
        
    except:
        print("please just write a right \"number\", not any thing else!")
        sys.exit(1)
    #Diameter L0               
    inputs.append(input("specify the diameter of object when enters in atmosphere of Earth in \"meter\".\nsome suggestions:.:Asteroids:. Itokawa(500 m), Apophis(325m), Bennu(250m), \n.:Comets:. Tempel1(6300 m), Wild2(4200m), Hartley(1500m), \n.:Past events:. Chelyabinsk(20m), Tunguska(50m,stone), Barringer(50m,iron), Ries(1500m) Chicxulub(1400m) \n.:Other things:. Wembley Stadium(320m), Double-decker bus(5m) and... \n").replace(" ",""))
    try:
        inputs[1] = float(inputs[1])
        if inputs[1]<0: 
            print("its can not be zero or negative! please write a right value.")
            sys.exit(1)
        
        if inputs[1]<50: print("tiny! it would not form any crater! maybe some fragments and small holes.\n")
        elif inputs[1]<100: print("it would not form any crater! just some fragments and small holes.\n")
        elif inputs[1]>1000: print("what a huge meteor!! God bless us!\n")        
    except: 
        print("please just write a right \"number\", not any thing else!")
        sys.exit(1)
    
    #velocity
    inputs.append(input("specify the velocity of the projectile in \"km/s\".\nsome suggestions: This is the velocity of the projectile before it enters the atmosphere. \nThe minimum impact velocity on Earth is \"11 km/s\".\nTypical impact velocities are \"17 km/s\" for asteroids and \"51 km/s\" for comets. \nThe maximum Earth impact velocity for objects orbiting the sun is \"72 km/s\".\n").replace(" ",""))
    try:
        inputs[2] = float(inputs[2])
        if inputs[2]<0: print("its can not be zero or negative! please write a right value.")
        if inputs[2]>=300000: 
            print("wooow!! excuse me. this speed is not possible! please reduce the speed...")
            sys.exit(1)
        if inputs[2]<1: print("its tired?! it would spend a lot if time to impact the earth!\n")
        elif inputs[2]>72: print("woow!! its more than the maximum Earth impact velocity of objects orbiting the sun!\n")        
    except:
        print("please just write a right \"number\", not any thing else!")
        sys.exit(1)
    
    #impactor Density
    inputs.append(input("specify the density of object in \"kg/m^3\". \nsome suggestions: ice(920 kg/m^3), porous rock(1600 kg/m^3), gabbro rock(3500 kg/m^3), iron(8000 kg/m^3)\n").replace(" ",""))
    try:
        inputs[3] = float(inputs[3])
        if inputs[3]<0: 
            print("its can not be zero or negative! please write a right value.")
            sys.exit(1)
        if inputs[3]<1000:print("its density is lower than water?! ok!\n")
        elif inputs[3]>8000 and inputs[3]<22590:print("what a huge density! it would be very strong and with high resistant.\n")
        elif inputs[3]>22590:print("it is denser than any solid matter we know! is it a new matter entering the earth?!\n")        
    except: 
        print("please just write a right \"number\", not any thing else!")
        sys.exit(1)
    
    #angle of object from a plane tangent to the impact surface    
    inputs.append(input("specify the angle of objects in \"degree\" between 0 and 90.\nguidance: The impact angle is measured from a plane tangent to the impact surface. \nThis angle is 90 degrees for a vertical impact. The most probable angle of impact is 45 degrees.\n").replace(" ",""))
    try:
        inputs[4] = float(inputs[4])
        if inputs[4]<=0 and inputs[4]>90: print("the angle shuld be between 0 and 90 in degree.")
    except: 
        print("please just write a right \"number\", not any thing else!")
        sys.exit(1)
    
    #Depth of water or target Density
    if inputs[0]==0:
        inputs.append(input("specify the depth of water in \"meter\".\nguidance: the average depth of oceans is about 3000 meter. the deepest place in oceans is about 11000 meter!\n").replace(" ",""))
        try:
            inputs[5] = float(inputs[5])
            if inputs[5]<0: 
                print("its can not be zero or negative! please write a right value.")
                sys.exit(1)
            if inputs[5]<=20: print("are you sure its sea?! i think it is a pool!\n")
            elif inputs[5]>=11000: print("Woow! it is deeper than deepest place of oceans in the world!\n")            
        except: 
            print("please just write a right \"number\", not any thing else!")
            sys.exit(1)
        
    if inputs[0]==1:
        inputs.append(input("specify the target density. please write a number between 2500(sedimentary rock) and 2750(crystalline rock)in \"kg/m^3\". \n").replace(" ",""))
        try:
            inputs[5] = float(inputs[5]) 
            if inputs[5]<0: 
                print("its can not be zero or negative! please write a right value.")
                sys.exit(1)
            if inputs[5]<2500 or inputs[5]>2750: 
                print("please a number between 2500 and 2700")
                sys.exit(1)
        except: 
            print("please just write a right \"number\", not any thing else!") 
            sys.exit(1)
    return inputs

#computing the Kinetic Energy of object with special ralativity equation
def energy_func(v0,rowi,L0):
    c = 3*(10**5)
    m = np.pi *rowi*(L0**3)/6
    gamma = 1/(np.sqrt(1-(v0**2/c**2)))
    E = m * c**2 *(gamma-1)
    E *= 10**6
    return int(E)   #megajoule

#computing whether the impactor begins to break up well above the surface or not
def If_func(v0,rowi,L0,theta):
    cD = 2
    H = 8000 #m
    Yi = 10**(2.107+(0.0624*np.sqrt(rowi)))
    If = (4.07*cD*H*Yi)/(rowi*L0* v0**2 *np.sin(np.deg2rad(theta)))
    return If

#computing the height when impactor begins to Break up 
def z_star_func(If,v0,rowi):
    if If>1:
        z_star = -1
        return z_star
    cD = 2
    H = 8000 #m
    row0 = 1
    Yi = 10**(2.107+(0.0624*np.sqrt(rowi)))
    z_star = -H * ( np.log(Yi/(row0 * v0**2)) + 1.308 - (0.314*If) - (1.303*np.sqrt(1-If)))
    return z_star #meter

#computing the height when critical impactor Diameter occurs 
def zb_func(z_star, L0, theta,rowi):
    cD = 2
    H = 8000 #m
    row0 = 1
    fp = 7
    row_z_star = row0 * np.exp(-z_star/H)
    l = L0 * np.sin(np.deg2rad(theta)) * np.sqrt(rowi/(cD*row_z_star))   #meter 
    zb = z_star - 2*H * np.log(1+(l*np.sqrt(fp**2 - 1)))
    return zb  #meter

#computing the increased Diameter of impactor or fragments of that, when impact the earth    
def L_func(z_star,zb,rowi,L0,theta):
    cD = 2
    H = 8000 #m
    row0 = 1
    row_z_star = row0 * np.exp(-z_star/H)
    l = L0 * np.sin(np.deg2rad(theta)) * np.sqrt(rowi/(cD*row_z_star))
    if z_star<0:#not break up before impact
        L = L0
    elif z_star>0 and zb<0:#not arrivig at critical impactor Diameter before impact
        L = L0 * np.sqrt(1+((2*H/l)**2 * (np.exp(z_star/(2*H))-1)**2))/np.sin(np.deg2rad(theta))
    elif z_star>0 and zb>0:
        L = L0 * np.sqrt(1+((2*H/l)**2 * (np.exp(z_star/(2*H))-1)**2))
    return L  #meter

#computing the velocity of impactor when impact the earth
def v_func(v0,rowi,L0,theta,z_star,zb,If):
    cD = 2
    H = 8000 #m
    fp = 7
    row0 = 1
    row_z_star = row0 * np.exp(-z_star/H)
    l = L0 * np.sin(np.deg2rad(theta)) * np.sqrt(rowi/(cD*row_z_star)) 
    if If>1 or z_star<0:
        v = v0 * np.exp(-(3*row0*cD*H)/(4*rowi*L0*np.sin(np.deg2rad(theta))))
    elif If<1 and zb<0:
        v_z_star = v0 * np.exp(-(3*row_z_star*cD*H)/(4*rowi*L0*np.sin(np.deg2rad(theta))))
        L_z_star = (H**3 * L0**2 /(3 * l**2)) * ( 3*(4+((l/H)**2))*np.exp(z_star/H) + (6* np.exp(2*z_star/H)) - (16 * np.exp((3/2)*z_star/H)) - (3* (l/H)**2) -2)
        v = v_z_star * np.exp((-(3*row_z_star*cD)/(4*rowi*L0**3 *np.sin(np.deg2rad(theta)))) * L_z_star)
    elif If<1 and zb>0:
        v_z_star = v0 * np.exp(-(3*row_z_star*cD*H)/(4*rowi*L0*np.sin(np.deg2rad(theta))))
        alpha = np.sqrt(fp**2 -1)
        L_zb =  (l* L0**2 *alpha/24) * (8*(3 + alpha**2)+ 3*alpha*l*(2+alpha**2)/H)
        v = v_z_star * np.exp((-(3*row_z_star*cD)/(4*rowi*L0**3 *np.sin(np.deg2rad(theta))) ) * L_zb)
    return v #km/s

#computing the velocity of impactor when arrives at the seafloor
def v_seafloor_func(v,dw,rowi,L,theta):
    if dw>=500:
        cD = 1
    elif dw<500 and dw>200:
        cD = 1.5
    else: cD = 2
    roww = 1000
    v_floor = v * np.exp(-((3/2)*(roww*cD*dw))/(rowi*L*np.sin(np.deg2rad(theta))))
    return v_floor

#computing the transient Diameter of crater
def Dtc_func(v,rowi,rowt,L,theta,floor_type):
    gE = 9.8
    v *= 1000
    Dtc = 1.161 * (rowi/rowt)**(1/3) * L**(0.78) * v**(0.44) * gE**(-0.22) * np.sin(np.deg2rad(theta))**(1/3)
    if floor_type==0:#sea
        rowt = 2700
        Dtc = 1.365 * (rowi/rowt)**(1/3) * L**(0.78) * v**(0.44) * gE**(-0.22) * np.sin(np.deg2rad(theta))**(1/3)
    return Dtc #meter

#computing the final Diameter of crater
def Dfr_func(Dtc):
    #simple crater
    if Dtc<=2560:
        Dfr = 1.25 * Dtc
        vbr = 0.032* Dfr**3
        dtc = Dtc/(2*np.sqrt(2))
        hfr = (0.07 * (Dtc**4))/(Dfr**3)
        tbr = 2.8 * vbr * ((dtc+hfr)/(dtc * Dfr**2))
        dfr = dtc+hfr-tbr
    #complex crater
    elif Dtc>2560:
        Dc = 3.2
        Dtc /= 1000
        Dfr = (1.17 * Dtc**1.13)/Dc**0.13
        dfr =  0.294 * Dfr**(0.3)
        dfr *= 1000        
        Dfr *= 1000
    return Dfr,dfr #meter

#generating data to visualize the Simple crater 
def simple_crater(Sigma,c,N):
    X = np.linspace(-34, 34, N)
    Y = np.linspace(-34, 34, N)
    X, Y = np.meshgrid(X, Y)
    
    # Mean vector and covariance matrix
    mu = np.array([1, 1])
    # Pack X and Y into a single 3-dimensional array
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    
    def multivariate_gaussian(pos, mu, Sigma):
        """Return the multivariate Gaussian distribution on array pos.
    
        pos is an array constructed by packing the meshed arrays of variables
        x_1, x_2, x_3, ..., x_k into its _last_ dimension.
    
        """
    
        n = mu.shape[0]
        Sigma_det = np.linalg.det(Sigma)
        Sigma_inv = np.linalg.inv(Sigma)
        M = np.sqrt((2*np.pi)**n * Sigma_det)
        # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
        # way across all the input variables.        
        fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)
    
        return np.exp(-fac / 10) / M
    
    # The distribution on the variables X, Y packed into pos.
    
    Z = multivariate_gaussian(pos, mu, Sigma)
    Z *= c
    Z = Z.tolist()
    return Z

#generating data to visualize the Complex crater   
def complex_crater(Sigma1,Sigma2,c1,c2,N):
    X = np.linspace(-34, 34, N)
    Y = np.linspace(-34, 34, N)
    X, Y = np.meshgrid(X, Y)
    
    # Mean vector and covariance matrix
    mu = np.array([1, 1])
    
    # Pack X and Y into a single 3-dimensional array
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    
    def multivariate_gaussian(pos, mu, Sigma):
        """Return the multivariate Gaussian distribution on array pos.
    
        pos is an array constructed by packing the meshed arrays of variables
        x_1, x_2, x_3, ..., x_k into its _last_ dimension.
    
        """
    
        n = mu.shape[0]
        Sigma_det = np.linalg.det(Sigma)
        Sigma_inv = np.linalg.inv(Sigma)
        M = np.sqrt((2*np.pi)**n * Sigma_det)

        # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
        # way across all the input variables.
        fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)
    
        return np.exp(-fac / 10) / M
    
    # The distribution on the variables X, Y packed into pos.
    
    Z1 = multivariate_gaussian(pos, mu, Sigma1)
    Z1 *= c1
    Z2 = multivariate_gaussian(pos, mu, Sigma2)
    Z2 *= c2
    Z = Z1+Z2
    Z = Z.tolist()
    return Z

#specifying and visualizing the crater whether it is Simple or Complex        
def specify_Q_Z(Dtc,dfr,L0,N):
     
    if L0<100:
        c = 0
        Sigma = np.array([[ 45 , 2.5], [2.5,  50]])
        Q = 80
        Z = simple_crater(Sigma, c, N)
    elif dfr<=150:
        c = -20
        Sigma = np.array([[ 45 , 2.5], [2.5,  50]])
        Q = 120
        Z = simple_crater(Sigma, c, N)
    elif Dtc<=3200:
            if dfr<200:
                c = -20
                Sigma = np.array( [[ 45 , 2.5], [2.5,  50]])
                Q = 120  
                Z = simple_crater(Sigma, c, N)
            elif dfr<250:
                c = -40
                Sigma = np.array(  [[ 45 , 2.5], [2.5,  50]])
                Q = 80  
                Z = simple_crater(Sigma, c, N)
            elif dfr<300:
                c = -80
                Sigma = np.array(  [[ 45 , 2.5], [2.5,  50]])
                Q = 40 
                Z = simple_crater(Sigma, c, N)
            elif dfr<350:
                c = -100
                Sigma = np.array(  [[ 55 , .5], [.5,  60]])
                Q = 30 
                Z = simple_crater(Sigma, c, N)
            elif dfr<400:
                c = -150
                Sigma = np.array(  [[ 60 , .5], [.5,  75]])
                Q = 25  
                Z = simple_crater(Sigma, c, N)
            elif dfr<450:
                c = -200
                Sigma = np.array(  [[ 60 , .5], [.5,  75]])
                Q = 20    
                Z = simple_crater(Sigma, c, N)
            elif dfr<500:
                c = -250
                Sigma = np.array(  [[ 60 , .5], [.5,  75]])
                Q = 20   
                Z = simple_crater(Sigma, c, N)
            elif dfr<550:
                c = -300
                Sigma = np.array(  [[ 60 , .5], [.5,  75]])
                Q = 15   
                Z = simple_crater(Sigma, c, N)
            elif dfr<600:
                c = -350
                Sigma = np.array(  [[ 60 , .5], [.5,  75]])
                Q = 15 
                Z = simple_crater(Sigma, c, N)
            elif dfr<700:
                c = -100
                Sigma = np.array(  [[ 30 , 2.5], [2.5,  35]])
                Q = 22   
                Z = simple_crater(Sigma, c, N)
            else :
                c = -200
                Sigma = np.array(  [[ 30 , 2.5], [2.5,  30]])
                Q = 17 
                Z = simple_crater(Sigma, c, N)
  
    elif Dtc<=7000:
        if dfr<620:
             
            c1 = -40
            c2 =  2
            Sigma1 = np.array(  [[ 65 , 2.5], [2.5,  75]])
            Sigma2 = np.array(  [[ 6 , .05], [.05,  6]])
            Q = 130
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)
        elif dfr<=700:
             
            c1 = -80
            c2 = -5
            Sigma1 = np.array([[ 60 , 2.5], [2.5,  70]])
            Sigma2 = np.array([[ 8 , .05], [.05,  8]])
            Q = 100
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)
        elif dfr<= 800:
             
            c1 = -200
            c2 =  13
            Sigma1 = np.array(  [[ 75 , 2.5], [2.5,  80]])
            Sigma2 = np.array(  [[ 10 , 2.5], [2.5,  10]])
            Q = 65
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)
        else:
             
            c1 = -500
            c2 =  28
            Sigma1 = np.array([[ 150 , 2.5], [2.5,  150]])
            Sigma2 = np.array([[ 20 , .5], [.5,  20]])
            Q = 50            
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)                                           
    else:
        if dfr<1200:
             
            c1 = -40
            c2 =  1.5
            Sigma1 = np.array(  [[ 65 , 2.5], [2.5,  75]])
            Sigma2 = np.array(  [[ 6 , .05], [.05,  6]])
            Q = 130
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)
        else:
             
            c1 = -80
            c2 =  5
            Sigma1 = np.array([[ 60 , 2.5], [2.5,  70]])
            Sigma2 = np.array([[ 8 , .05], [.05,  8]])
            Q = 100
            Z = complex_crater(Sigma1,Sigma2,c1,c2,N)            
            
       
    return Q,Z                   

#generating noise for realize and idealize the visualization 
def values_octaves(Q,N):
    def multiple_octaves(octaves, start_amplitude,Q,N):
        WIDTH = N
        HEIGHT = N
        parameters = []
        for i in range(octaves):
            parameters.append({
                'offset': random.random() * 2 * math.pi,
                'frequency': 2**i,
                'amplitude': start_amplitude / float(i+1),
            })
    
        def noise(x, y):
            value = 0
            for p in parameters:
                x_part = math.sin(
                    (x / float(WIDTH))
                    * p['frequency']
                    * 2 * math.pi
                    + p['offset']
                )
                y_part = math.sin(
                    (y / float(HEIGHT))
                    * p['frequency']
                    * 2 * math.pi
                    + (x_part % (2 * math.pi))
                )
                value += y_part * p['amplitude']
    
            return value/Q
    
        values1 = []
        for x in range(WIDTH):
             for y in range(WIDTH):
                 values1.append(noise(x,y))
        
    
        return values1

    value = multiple_octaves(5, .9, Q, N)

    values2 = []
    i = 0
    j = N
    for c in range(N):
        values2.append(value[i:j])
        i+=N
        j+=N
       
    return values2
        
        
        
        
        
        
  
        
                
        
        
        
        
        
        
        
        
        
    
            
    
    
    
    
    
    
    
    
    
