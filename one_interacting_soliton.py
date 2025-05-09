
#This code visualises the exact solution of the KPE and shows the one soliton solution in 3D while also providing its 
#properties, which include Amplitude, Speed and Angle, created by Shoaib Mohammed

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

def A(p,q,r):
    A_pqr = (r-q)*np.power(p,2) + (p-r)*np.power(q,2) + (q-p)*np.power(r,2)
    return A_pqr

def generate(X,Y,tau):
    theta1 = k1*X + np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + np.power(k2,2)*Y - np.power(k2,3)*tau
    K = np.exp(theta1)+\
        a*np.exp(theta2)
	
    K_X = (k1)*np.exp(theta1)+\
        (k2)*a*np.exp(theta2)
                      
    K_XX = (k1)**2*np.exp(theta1)+\
        (k2)**2*a*np.exp(theta2)
            
    Z = 2 * ( (K_XX/K) - np.power(K_X/K,2) )
    return Z

def soliton(ki,kj):
    amp = 0.5*np.power(kj-ki,2)
    speed = ki*ki+ki*kj+kj*kj
    angle = np.arctan(-(ki+kj))*180/np.pi # degree
    return amp,speed,angle

fig = plt.figure(figsize = (8, 9))
gs = fig.add_gridspec(nrows=6,ncols=20)
ax = fig.add_subplot(gs[:-1,:-2])
cbar_ax = fig.add_subplot(gs[1:-2,-1])


k1 = -0.25 #alter this and k2 to observe the solitons at different angles, amplitudes and speeds
k2 = 0.75
print(soliton(k1,k2))
a = 1 #alter this for phase shift, doesn't really do much otherwise
        
tau= -4 # start time
dt = 1 #change in time
Ndt = 10 #number of changes in time, includes line 74 as first

L = 20 #the (X,Y) domain
dx = 0.1
dy = dx
Ystar = 0.0
X1 = -L
X2 = L

ax.set_xlim(X1, X2)
ax.set_ylim(Ystar-L,Ystar+ L)

ax.set_xlabel('X')
ax.set_ylabel('Y')

lw = 1 # line width
lw2 = 4

# Make data.
x = np.arange(-L, L, dx)
xp= np.arange(0, L, dx)
xn= np.arange(Ystar-L,0, dx)
y = np.arange(Ystar-L,Ystar+L, dy)
X, Y = np.meshgrid(x, y)


time_text = ax.text(0.3, 1.05, '', transform=ax.transAxes, fontsize='xx-large')

fac = 0.25
levels = np.linspace(0, 2.0, 21) # fac*[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
print('levels',levels)

wframe = None

for frame in range(Ndt): # range(n): number of frames
    
    if wframe:
        wframe.remove()
        cbar_ax.cla()
        line12p.remove()
        p12p.remove()
    
    # Plot the new wireframe and pause briefly before continuing.
    Z = generate(X, Y, tau)
    wframe = ax.contourf(X, Y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
    
    time_text.set_text(r'$u(X,Y,\tau=%.2f)$' %tau)
    
   
    # Settings for Y or lambda SP1 solution
    yp12 = facc*(-(k1-k2)*x+(k1**3-k2**3)*tau+np.log(a))/(k1**2-k2**2)
    line12p, = ax.plot(x, yp12, '--k', linewidth=lw)
    p12p, = ax.plot((-(k1**2-k2**2)*L+(k1**3-k2**3)*tau-np.log(a))/(k1-k2), facc*L, 'xk', linewidth=lw2)
    
    plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$u$')
    
    plt.pause(1)
    tau=tau+dt

plt.show()
