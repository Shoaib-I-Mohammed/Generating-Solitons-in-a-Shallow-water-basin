# This code visualises the exact solution of the KPE that describes the interaction 
# of three solitary waves in 3D, provided by Prof. Onno Bokhove and modified by Shoaib Mohammed

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

def B(p,q,r):
    B_pqr = (r-q)*np.power(p,2) + (p-r)*np.power(q,2) + (q-p)*np.power(r,2)
    return B_pqr

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

def generatePsi(X,Y,tau,facc):
    theta1 = k1*X + facc*np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + facc*np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + facc*np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + facc*np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + facc*np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + facc*np.power(k6,2)*Y - np.power(k6,3)*tau
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
    
    Psi = 2 * ( (K_X/K) )
    return Psi

def soliton(ki,kj):
    amp = 0.5*np.power(kj-ki,2)
    speed = ki*ki+ki*kj+kj*kj
    angle = np.arctan(-(ki+kj))*180/np.pi # degree
    return amp,speed,angle

fig = plt.figure(figsize = (8, 9))
gs = fig.add_gridspec(nrows=6,ncols=20)
ax = fig.add_subplot(gs[:-1,:-2])
cbar_ax = fig.add_subplot(gs[1:-2,-1])

delt = 10**(-6)
Atilde = 0.5
k4 = np.sqrt(0.5*Atilde) # 0.5
k5 = np.sqrt(Atilde)*(np.sqrt(0.5))+delt # 0.50000071
k6 = np.sqrt(Atilde)*(np.sqrt(2)+np.sqrt(0.5))+delt # 1.50000071
k1 = -k6
k2 = -k5
k3 = -k4
a = ((k6*(k6**2-k4**2))/(k5*(k5**2-k4**2)))**(1/2)
b = 1
c = 1/a
Ystar = 0.5*np.log( a*k5*(k5**2-k4**2)/(c*k6*(k6**2-k4**2)) )/(k6**2-k5**2)
tau= -30 # start time
dt = 5
Ndt = 11


A_135 =       B(k1,k3,k5)
A_235 =     a*B(k2,k3,k5)
A_136 =     c*B(k1,k3,k6)
A_236 =   a*c*B(k2,k3,k6)
A_145 =     b*B(k1,k4,k5)
A_245 =   a*b*B(k2,k4,k5)
A_146 =   b*c*B(k1,k4,k6)
A_246 = a*b*c*B(k2,k4,k6)
print('As Ystar',A_135,A_235,A_136,A_236,A_145,A_245,A_146,A_246,Ystar)
print('k4,k5,k6',k4,k5,k6)

L = 50
dx = 0.1
dy = dx

ax.set_xlim(-L, L)
ax.set_ylim(Ystar-L,Ystar+ L)

ax.set_xlabel('X')
ax.set_ylabel('Y')

lw = 1 # line width


# Make data.
x = np.arange(-L, L, dx)
xp= np.arange(0, L, dx)
xn= np.arange(Ystar-L,0, dx)
y = np.arange(Ystar-L,Ystar+L, dy)
X, Y = np.meshgrid(x, y)

time_text = ax.text(0.3, 1.05, '', transform=ax.transAxes, fontsize='xx-large')

fac = 0.25
levels = np.linspace(0.0, 9*Atilde, 17) # fac*[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
print('levels',levels)

wframe = None

for frame in range(Ndt): # range(n): number of frames
    if wframe:
        wframe.remove()
        cbar_ax.cla()
        
    
    # Plot the new wireframe and pause briefly before continuing.
    Z = generate(X, Y, tau)
    wframe = ax.contourf(X, Y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
    
    time_text.set_text(r'$u(X,Y,\tau=%.2f)$' %tau)
    plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$u$')
    
    plt.pause(1)
    tau=tau+dt

plt.show()
