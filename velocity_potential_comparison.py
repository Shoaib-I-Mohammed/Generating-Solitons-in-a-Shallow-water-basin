#This code visualises the velocity potential for the three soliton solution being
#generated in the rectangular basin. It looks at the aymptotic pertubation of phi
#up to epsilon**2. Produced by Shoaib Mohammed

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation
from math import floor, log10, isnan

def B(p,q,r): #n>=0
    B_pqr = (r-q)*np.power(p,2) + (p-r)*np.power(q,2) + (q-p)*np.power(r,2)
    return B_pqr

def KX(n,X,Y,tau): #generate will not be using this!
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + np.power(k6,2)*Y - np.power(k6,3)*tau
    
    if n%2==1:
        K_Xn = -(3*(2*Ahat)**(1/2)/2)**(n)*B_246*np.sinh((theta1+theta3+theta5-theta2-theta4-theta6)/2)\
             -((2*Ahat)**(1/2)/2)**(n)*B_145*np.sinh((theta1+theta4+theta5-theta2-theta3-theta6)/2)\
             +((2*Ahat)**(1/2)/2)**(n)*(B_235*B_146)**(1/2)*(-np.sinh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                      np.sinh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    
    else:
        K_Xn = (3*(2*Ahat)**(1/2)/2)**(n)*B_246*np.cosh((theta1+theta3+theta5-theta2-theta4-theta6)/2)\
             +((2*Ahat)**(1/2)/2)**(n)*B_145*np.cosh((theta1+theta4+theta5-theta2-theta3-theta6)/2)\
             +((2*Ahat)**(1/2)/2)**(n)*(B_235*B_146)**(1/2)*(np.cosh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                      np.cosh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    return K_Xn

def KY(m,X,Y,tau): #m>=1
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + np.power(k6,2)*Y - np.power(k6,3)*tau
    
    if m%2==1:     
        K_Yn = (k6**2-k5**2)**m*(B_235*B_146)**(1/2)*(np.sinh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                 np.sinh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    else:
        K_Yn = (k6**2-k5**2)**m*(B_235*B_146)**(1/2)*(np.cosh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                 np.cosh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    return K_Yn

def K_XY(n,m,X,Y,tau): #n>=0,m>=1
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + np.power(k6,2)*Y - np.power(k6,3)*tau    
    
    if (n+m)%2 == 1:     
        K_XYnm = ((2*Ahat)**(1/2)/2)**(n)*(k6**2 - k5**2)**m*(B_235*B_146)**(1/2)*((-1)**n*np.sinh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                 np.sinh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    else:
        K_XYnm = ((2*Ahat)**(1/2)/2)**(n)*(k6**2 - k5**2)**m*(B_235*B_146)**(1/2)*((-1)**n*np.cosh((theta1+theta3+theta6-theta2-theta4-theta5)/2)+\
                                                 np.cosh((theta1+theta4+theta6-theta2-theta3-theta5)/2))
    return K_XYnm

def phi(x,y,z,t):
    K=KX(0,x,y,t)
    
    psi = 6*(2**(1/2)/3)**(5/3)*KX(1,x,y,t)/K
    
    psixx = 2*2**(1/2)*(KX(3,x,y,t)/K - 3*KX(2,x,y,t)*KX(1,x,y,t)/(K)**(2) + 2*(KX(1,x,y,t)/K)**3)
    
    psixxxx = 6*(2**(1/2)/3)**(1/3)*(KX(5,x,y,t)/K\
                                     - 5*KX(4,x,y,t)*KX(1,x,y,t)/(K)**2 - 10*KX(3,x,y,t)*KX(2,x,y,t)/(K)**2\
                                     + 20*KX(3,x,y,t)*(KX(1,x,y,t)**2)/(K)**3 + 30*(KX(2,x,y,t)**2)*KX(1,x,y,t)/(K)**3\
                                     - 60*KX(2,x,y,t)*(KX(1,x,y,t)**3)/(K)**4 + 24*(KX(1,x,y,t)/K)**5)   
    
    psiyy = 6*(2**(1/2)/3)**(1/3)*eps*(K_XY(1,2,x,y,t)/K - 2*K_XY(1,1,x,y,t)*KY(1,x,y,t)/(K)**2\
                                     - KY(2,x,y,t)*KX(1,x,y,t)/(K)**2 + 2*KX(1,x,y,t)*(KY(1,x,y,t))**2/(K)**3)
    
    psiyyyy = 9*(2)**(1/2)*eps**2*(K_XY(1,4,x,y,t)/K - 4*K_XY(1,3,x,y,t)*KY(1,x,y,t)/(K)**2 - KY(4,x,y,t)*KX(1,x,y,t)/(K)**2\
                                  - 4*K_XY(1,1,x,y,t)*KY(3,x,y,t)/(K)**2 - 6*K_XY(1,2,x,y,t)*KY(2,x,y,t)/(K)**2 + 12*K_XY(1,2,x,y,t)*(KY(1,x,y,t)**2)/(K)**3\
                                  + 8*KY(3,x,y,t)*KY(1,x,y,t)*KX(1,x,y,t)/(K)**3 + 24*K_XY(1,1,x,y,t)*KY(2,x,y,t)*KY(1,x,y,t)/(K)**3\
                                  + 6*KY(2,x,y,t)**2*KX(1,x,y,t)/(K)**3 - 24*K_XY(1,1,x,y,t)*KY(1,x,y,t)**3/(K)**4 - 36*KY(2,x,y,t)*KY(1,x,y,t)**2*KX(1,x,y,t)/K**4\
                                  + 24*KY(1,x,y,t)**4*KX(1,x,y,t)/(K)**5)
    
    
    phi1=(c*atilde*l/h)*(psi-eps*z**2*(psixx+psiyy)/2+eps**2*z**4*(psixxxx+psiyyyy)/24)
    phi2=(c*atilde*l/h)*(psi-eps*z**2*(psixx+psiyy)/2)
    phi3=(c*atilde*l/h)*(psi)
    return phi1, phi2, phi3

def generate(X,Y,tau):
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + np.power(k6,2)*Y - np.power(k6,3)*tau

    K = A_135*np.exp(theta1+theta3+theta5)+\
        A_235*np.exp(theta2+theta3+theta5)+\
        A_136*np.exp(theta1+theta3+theta6)+\
        A_236*np.exp(theta2+theta3+theta6)+\
        A_145*np.exp(theta1+theta4+theta5)+\
        A_245*np.exp(theta2+theta4+theta5)+\
        A_146*np.exp(theta1+theta4+theta6)+\
        A_246*np.exp(theta2+theta4+theta6)
      
    K_X = A_135*(k1+k3+k5)*np.exp(theta1+theta3+theta5)+\
          A_235*(k2+k3+k5)*np.exp(theta2+theta3+theta5)+\
          A_136*(k1+k3+k6)*np.exp(theta1+theta3+theta6)+\
          A_236*(k2+k3+k6)*np.exp(theta2+theta3+theta6)+\
          A_145*(k1+k4+k5)*np.exp(theta1+theta4+theta5)+\
          A_245*(k2+k4+k5)*np.exp(theta2+theta4+theta5)+\
          A_146*(k1+k4+k6)*np.exp(theta1+theta4+theta6)+\
          A_246*(k2+k4+k6)*np.exp(theta2+theta4+theta6)
    
    K_XX= A_135*np.power(k1+k3+k5,2)*np.exp(theta1+theta3+theta5)+\
          A_235*np.power(k2+k3+k5,2)*np.exp(theta2+theta3+theta5)+\
          A_136*np.power(k1+k3+k6,2)*np.exp(theta1+theta3+theta6)+\
          A_236*np.power(k2+k3+k6,2)*np.exp(theta2+theta3+theta6)+\
          A_145*np.power(k1+k4+k5,2)*np.exp(theta1+theta4+theta5)+\
          A_245*np.power(k2+k4+k5,2)*np.exp(theta2+theta4+theta5)+\
          A_146*np.power(k1+k4+k6,2)*np.exp(theta1+theta4+theta6)+\
          A_246*np.power(k2+k4+k6,2)*np.exp(theta2+theta4+theta6)

    Z = 2 * ( (K_XX/K) - np.power(K_X/K,2) )
    return Z

def Psigen(X,Y,tau):
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + np.power(k6,2)*Y - np.power(k6,3)*tau

    K = A_135*np.exp(theta1+theta3+theta5)+\
        A_235*np.exp(theta2+theta3+theta5)+\
        A_136*np.exp(theta1+theta3+theta6)+\
        A_236*np.exp(theta2+theta3+theta6)+\
        A_145*np.exp(theta1+theta4+theta5)+\
        A_245*np.exp(theta2+theta4+theta5)+\
        A_146*np.exp(theta1+theta4+theta6)+\
        A_246*np.exp(theta2+theta4+theta6)
      
    K_X = A_135*(k1+k3+k5)*np.exp(theta1+theta3+theta5)+\
          A_235*(k2+k3+k5)*np.exp(theta2+theta3+theta5)+\
          A_136*(k1+k3+k6)*np.exp(theta1+theta3+theta6)+\
          A_236*(k2+k3+k6)*np.exp(theta2+theta3+theta6)+\
          A_145*(k1+k4+k5)*np.exp(theta1+theta4+theta5)+\
          A_245*(k2+k4+k5)*np.exp(theta2+theta4+theta5)+\
          A_146*(k1+k4+k6)*np.exp(theta1+theta4+theta6)+\
          A_246*(k2+k4+k6)*np.exp(theta2+theta4+theta6)
    
    
    Z = 2*K_X/K 
    return Z



'Step 1: Choosing parameters'
delt = 10**(-3) #alter this for more accuracy and phase shift, do not set to <=0
Ahat = 1 #alter this to alter k_1,k_2 etc, cannot be too big otherwise the plot wont work

choice=4 #keep choice 4 to see results in report, if desired the user can simulate the other values, but they will not be done 
         #in their respective basin
if choice==1:
    l=2 #lamda, i.e. wavelength of wave eta, should be larger than h in meters
    h=0.25 #depth of the water body in meters
    BasinL=30 #length of basin's wall in meters, ensure that BasinL>BasinW
    BasinW=25 #width of basin's wall in meters
elif choice==2:
    l=1.1
    h=np.sqrt(0.5) 
    BasinL=30 
    BasinW=25 
elif choice==3: #for the tabletop square basin cm
    l=13.5 
    h=7.5
    BasinL=250
    BasinW=250
elif choice==4: #rectangle basin, also the dimensions for shifting axis
    l=0.95
    h=0.35
    BasinL=30
    BasinW=27 
elif choice==5: #for circle basin
    BasinL=12.5 #length of basin's wall in meters, ensure that BasinL>BasinW
    BasinW=BasinL #width of basin's wall in meters    
    h=1
    l=1.4



#do not alter anything from below up until Step 2
k4 = np.sqrt(0.5*Ahat)  
k5 = np.sqrt(Ahat)*(np.sqrt(0.5)+delt)
k6 = np.sqrt(Ahat)*(np.sqrt(2)+np.sqrt(0.5)+delt)
k1 = -k6
k2 = -k5
k3 = -k4
a = ((k6*(k6**2-k4**2))/(k5*(k5**2-k4**2)))**(1/2)
b = 1
c = 1/a
Ystar = 0.5*np.log( a*k5*(k5**2-k4**2)/(c*k6*(k6**2-k4**2)) )/(k6**2-k5**2)
A_135 =       B(k1,k3,k5)
A_235 =     a*B(k2,k3,k5)
A_136 =     c*B(k1,k3,k6)
A_236 =   a*c*B(k2,k3,k6)
A_145 =     b*B(k1,k4,k5)
A_245 =   a*b*B(k2,k4,k5)
A_146 =   b*c*B(k1,k4,k6)
A_246 = a*b*c*B(k2,k4,k6)

B_135 =       B(k1,k3,k5)
B_235 =       B(k2,k3,k5)
B_136 =       B(k1,k3,k6)
B_236 =       B(k2,k3,k6)
B_145 =       B(k1,k4,k5)
B_245 =       B(k2,k4,k5)
B_146 =       B(k1,k4,k6)
B_246 =       B(k2,k4,k6)


eps=1/(l/h)**2 #epsilon
atilde=h*eps #dimensional amplitude of wave eta in meters
Abar=3*(2**(1/2)/3)**(4/3)*atilde*Ahat #amplitude of eta in meters
g=10 #gravitational acceleration, can change to 9.8 for more accurate results m/s^2
c=(g*h)**(1/2) #vertical velocity m/s
print('wavelength=%.2f'%l, 'amplitude=%.2f'%Abar, 'height=%.2f'%h,'epsilon=%.2f'%eps)
tau=-10 #start time (dimensionless)
dx = 0.1
dy = dx
dz=dx



'step 2: choose a value of tau'

rota=np.arctan(BasinL/BasinW)#to scale the x,y axis and rotate the basin
vart=np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))
taumeet1=-eps*((3*np.sqrt(8))**(1/3)*BasinL*(k6-k5)*(1/np.cos(vart))-4**(1/3)*l*np.log(B_246**2/(B_146*B_235)))/(2*l*(3**(1/3)*(k6-k5)+4**(1/3)*eps*(k6**3-k5**3)))
tmeet1=l*taumeet1/(2**(1/2)*eps*c) #meeting time for [1,2] and [3,4] soliton
taumeet2=(4**(1/3)*l*((k5-k6)**(-1)*(np.cos(vart)**2)*np.log(B_246**2/(B_235*B_146))-k3**(-1)*np.cos(2*vart)*np.log(B_246/B_236))+(3*8**(1/2))**(1/3)*BasinW*np.sin(vart))/(2*eps**(-1)*l*(4**(1/3)*eps*(k3**2*np.cos(2*vart)-(k6**2+k6*k5+k5**2)*np.cos(vart)**2)-3**(1/3)*np.sin(vart)**2))
tmeet2=l*taumeet2/(2**(1/2)*eps*c) #meeting time for [3,4] and [5,6] soliton

#choose tau below to view [3,4] soliton beginning to generate
tau=-(3*np.sqrt(8))**(1/3)*eps*(BasinW*np.sin(vart)+BasinL*np.cos(vart))/(2*l*(3**(1/3)+4**(1/3)*k3**2*eps))-4**(1/3)*eps*np.log(B_246/B_236)/(2*k3*(3**(1/3)+4**(1/3)*k3**2*eps))     

#choose tau below to view [1,2] soliton beginning to generate and meeting with [3,4] soliton
#tau=taumeet1

#choose tau below to view [5,6] soliton beginning to generate
#tau=eps*(3**(1/3)*np.sqrt(2)*(k6-k5)*(BasinL*np.cos(2*vart)-BasinW*np.sin(2*vart))+4**(1/3)*l*np.log(B_246**2/(B_235*B_146))*np.cos(vart))/(2*l*np.cos(vart)*(3**(1/3)*(k6-k5)+4**(1/3)*eps*(k6**3-k5**3)))

#choose tau below to view [3,4] and [5,6] soliton meeting at wall two
#tau=taumeet2

Ndtau=11
dtau=-tau/(Ndtau-3)

'step 3, choose a gen value from 1-5 more information below'
#1 to show a 3D contour plot of phi3, 2 to show a 4D scatter contour plot of phi2 from h<z<=h+gen(x,y,t), 
#3 to see same plot 4D plot of phi2 but from from 0<=z<=h, 4 to show a 4D plot of phi1 from h<z<=h+gen(x,y,t)
#and 5 to see the plot of phi3 from 0<=z<=h
gen=5
if gen==1: 
    fig = plt.figure(figsize = (8, 9))
    gs = fig.add_gridspec(nrows=6,ncols=20)
    ax = fig.add_subplot(gs[:-1,:-2])
    cbar_ax = fig.add_subplot(gs[1:-2,-1])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    time_text = ax.text(0.3, 1.05, '', transform=ax.transAxes, fontsize='xx-large')
    lw=1
    fac = 0.25
    levels = np.linspace(-2, 2, 17) # fac*[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
    print('levels',levels)
    ax.set_xlim(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1))
    ax.set_ylim(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1))
    x_bds = np.arange(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), dx)
    y_bds = np.arange(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1),dy)

    x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain

    wframe = None
    
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
            l34.remove()
            line12p.remove()
            line56p.remove()
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
        Z = phi((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l,1+eps*(4/3)**(1/3)*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l),eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l,tau), tau)[2] 
        #Z = (c*atilde*l/h)*3*(2**(1/2)/3)**(5/3)*Psigen((3/(2)**(1/2))**(1/3)*((x-c*t)), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y, tau)
        wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
            
        square = patches.Rectangle((-BasinL/2, -BasinW/2), BasinL, BasinW, edgecolor='red', facecolor='none')#dimensions of the basin        
        rotate = mpl.transforms.Affine2D().rotate_deg(-np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))*180/np.pi) + ax.transData #the appropriate rotation to obtain the ninefold amplitude in the basin
        square.set_transform(rotate)
        ax.add_patch(square)
        time_text.set_text(r'$\phi$(x,y,z,t=%.2fs)' %t)


        l34=ax.axvline(x=(1+np.sqrt(2)*(np.sqrt(2)/3)**(1/3)*eps*k3**2)*c*t+l*np.sqrt(2)**(1/3)*np.log(A_246/A_236)/(24**(1/3)*k3), color='black', linestyle='--', linewidth=lw)
        yp12p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)-(k6**3-k5**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_146))/((k6**2-k5**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line12p, = ax.plot(x_bds, yp12p, '--k', linewidth=lw)
        yp56p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)+(k5**3-k6**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_245))/((k5**2-k6**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line56p, = ax.plot(x_bds, yp56p, '--k', linewidth=lw)

        plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\phi$  (m$^2$/s)')
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau

        plt.show()
        
elif gen==2:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    gs = fig.add_gridspec(nrows=6,ncols=20)
    cbar_ax = fig.add_subplot(gs[1:-2,-1])
    time_text = ax.text(0.3, 1.05, 1, '', transform=ax.transAxes, fontsize='xx-large')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')

    x = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 500)
    y = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 500)
   
    X, Y = np.meshgrid(x, y)
    wframe = None
    bar=None
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
 
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
        Z=h+atilde*(4/3)**(1/3)*generate(X,Y,t)
        A = phi((3/(2)**(1/2))**(1/3)*((X-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,1+eps*(4/3)**(1/3)*generate((3/(2)**(1/2))**(1/3)*((X-c*t)/l),eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,tau), tau)[0] 
        x_flat = X.flatten()
        y_flat = Y.flatten()
        z_flat = Z.flatten()
        phi_flat = A.flatten()
        
        
        wframe = ax.scatter(x_flat, y_flat, z_flat,vmin=-1.5 , vmax=1.5, c=phi_flat, cmap='viridis',alpha=0.5)
        ax.set_title(r'$\phi$(x,y,z,t=%.2fs)' %t)
        ax.title.set_size(20)
        bar=plt.colorbar(wframe, cax=cbar_ax, orientation='vertical', label='$\phi$ (m$^2$/s)')
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau
        plt.show()
        
        
        
elif gen==3: #phi1 aith z=h+atilde*(4/3)**(1/3)*generate(X,Y,tau)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    gs = fig.add_gridspec(nrows=6,ncols=20)
    cbar_ax = fig.add_subplot(gs[1:-2,-1])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    time_text = ax.text(0.3, 1.05, 1, '', transform=ax.transAxes, fontsize='xx-large')
    x = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 50)
    y = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 50)
    z = np.linspace(0,h,50)
    X, Y, Z = np.meshgrid(x, y, z)
    wframe = None
    bar=None
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
 
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
       
        A = phi((3/(2)**(1/2))**(1/3)*((X-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,z/h, tau)[0] 
        x_flat = X.flatten()
        y_flat = Y.flatten()
        z_flat = Z.flatten()
        phi_flat = A.flatten()
        

        wframe = ax.scatter(x_flat, y_flat, z_flat,vmin=-1.5 , vmax=1.5, c=phi_flat, cmap='viridis',alpha=0.5)
        ax.set_title(r'$\phi$(x,y,z,t=%.2fs)' %t)
        ax.title.set_size(20)
        bar=plt.colorbar(wframe, cax=cbar_ax, orientation='vertical', label='$\phi$ (m$^2$/s)')
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau
        plt.show()
elif gen==4: #phi1 aith z=h+atilde*(4/3)**(1/3)*generate(X,Y,tau)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    gs = fig.add_gridspec(nrows=6,ncols=20)
    cbar_ax = fig.add_subplot(gs[1:-2,-1])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    time_text = ax.text(0.3, 1.05, 1, '', transform=ax.transAxes, fontsize='xx-large')
    x = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 500)
    y = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 500)
   
    X, Y = np.meshgrid(x, y)
    wframe = None
    bar=None
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
 
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
        Z=h+atilde*(4/3)**(1/3)*generate(X,Y,t)
        A = phi((3/(2)**(1/2))**(1/3)*((X-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,1+eps*(4/3)**(1/3)*generate((3/(2)**(1/2))**(1/3)*((X-c*t)/l),eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,tau), tau)[1] 
        x_flat = X.flatten()
        y_flat = Y.flatten()
        z_flat = Z.flatten()
        phi_flat = A.flatten()
        
        
        wframe = ax.scatter(x_flat, y_flat, z_flat,vmin=-1.5 , vmax=1.5, c=phi_flat, cmap='viridis',alpha=0.5)
        ax.set_title(r'$\phi$(x,y,z,t=%.2fs)' %t)
        ax.title.set_size(20)
        bar=plt.colorbar(wframe, cax=cbar_ax, orientation='vertical', label='$\phi$ (m$^2$/s)')
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau
        plt.show()
        
        
elif gen==5:     
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    gs = fig.add_gridspec(nrows=6,ncols=20)
    cbar_ax = fig.add_subplot(gs[1:-2,-1])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    time_text = ax.text(0.3, 1.05, 1, '', transform=ax.transAxes, fontsize='xx-large')
    x = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 50)
    y = np.linspace(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), 50)
    z = np.linspace(0,h,50)
    X, Y, Z = np.meshgrid(x, y, z)
    wframe = None
    bar=None
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
 
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
       
        A = phi((3/(2)**(1/2))**(1/3)*((X-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*Y/l,z/h, tau)[1] 
        x_flat = X.flatten()
        y_flat = Y.flatten()
        z_flat = Z.flatten()
        phi_flat = A.flatten()
        
        
        wframe = ax.scatter(x_flat, y_flat, z_flat,vmin=-1.5 , vmax=1.5, c=phi_flat, cmap='viridis',alpha=0.5)
        ax.set_title(r'$\phi$(x,y,z,t=%.2fs)' %t)
        ax.title.set_size(20)
        bar=plt.colorbar(wframe, cax=cbar_ax, orientation='vertical', label='$\phi$ (m$^2$/s)')
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau
        plt.show()    
    














