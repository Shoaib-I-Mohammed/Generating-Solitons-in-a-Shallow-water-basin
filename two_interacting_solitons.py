#
#
#This code was initially provided by Prof. Onno Bokhove and altered by Shoaib Mohammed. It visualises the
#exact solution of the KPE, which shows the two-soliton solution in 3D

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation

def A(p,q,r):
    A_pqr = (r-q)*np.power(p,2) + (p-r)*np.power(q,2) + (q-p)*np.power(r,2)
    return A_pqr

def generate(X,Y,tau,facc):
    theta1 = k1*X + facc*np.power(k1,2)*Y - np.power(k1,3)*tau
    theta2 = k2*X + facc*np.power(k2,2)*Y - np.power(k2,3)*tau
    theta3 = k3*X + facc*np.power(k3,2)*Y - np.power(k3,3)*tau
    theta4 = k4*X + facc*np.power(k4,2)*Y - np.power(k4,3)*tau
    theta5 = k5*X + facc*np.power(k5,2)*Y - np.power(k5,3)*tau
    theta6 = k6*X + facc*np.power(k6,2)*Y - np.power(k6,3)*tau
    K = (k3-k1)*np.exp(theta1+theta3)+\
        a*(k3-k2)*np.exp(theta2+theta3)+\
	b*(k4-k1)*np.exp(theta1+theta4)+\
        a*b*(k4-k2)*np.exp(theta2+theta4)

    K_X = (k1+k3)*(k3-k1)*np.exp(theta1+theta3)+\
        (k2+k3)*a*(k3-k2)*np.exp(theta2+theta3)+\
        (k1+k4)*b*(k4-k1)*np.exp(theta1+theta4)+\
        (k2+k4)*a*b*(k4-k2)*np.exp(theta2+theta4)
    
    K_XX = (k1+k3)**2*(k3-k1)*np.exp(theta1+theta3)+\
        (k2+k3)**2*a*(k3-k2)*np.exp(theta2+theta3)+\
        (k1+k4)**2*b*(k4-k1)*np.exp(theta1+theta4)+\
        (k2+k4)**2*a*b*(k4-k2)*np.exp(theta2+theta4)
    
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

facc = 1
delt = 10**(-2) #alter size for bigger amplitude and bigger phaseshift
k1 = -1
k2 = -delt
k3 = -k2
k4 = -k1
k5 = 1
k6 = 1
a= (((k3-k1)*(k4-k1))/((k3-k2)*(k4-k2)))**(1/2) #so fourfold amplitude is located at the centre
b= (((k3-k1)*(k3-k2))/((k4-k1)*(k4-k2)))**(1/2) #so fourfold amplitude is located at the centre
#a=1 #for the Y-wave
#b=1 #for the Y-wave

    
    

tau= -5 # start time
dt = 1
Ndt = 11
print('k1,k2,k3,k4',k1,k2,k3,k4,a,b)


L = 20
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
lw2 = 2

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
        line12n.remove()
        line34p.remove()
        line34n.remove()
        #p12p.remove()
        #p12n.remove()
        #p34p.remove()
        #p34n.remove()
    # Plot the new wireframe and pause briefly before continuing
    Z = generate(X, Y, tau,facc)
    wframe = ax.contourf(X, Y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)

    
    time_text.set_text(r'$u(X,Y,\tau=%.2f)$' %tau)
    
    
    #Centerlines for the line soliton solutions, with both phase shifts
    yp12p = ((-k1+k2)*x+(k1**3-k2**3)*tau+np.log(a*(k4-k2)/(k4-k1)))/(k1**2-k2**2)
    line12p, = ax.plot(x, yp12p, '--k', linewidth=lw)
    yp12n = facc*((-k1+k2)*x+(k1**3-k2**3)*tau+np.log(a*(k3-k2)/(k3-k1)))/(k1**2-k2**2)
    line12n, = ax.plot(x, yp12n, '--k', linewidth=lw)
    yp34p = facc*((-k3+k4)*x+(k3**3-k4**3)*tau-np.log(b*(k4-k1)/(k3-k1)))/(k3**2-k4**2)
    line34p, = ax.plot(x, yp34p, ':k', linewidth=lw)
    yp34n = facc*((-k3+k4)*x+(k3**3-k4**3)*tau-np.log(b*(k4-k2)/(k3-k2)))/(k3**2-k4**2)
    line34n, = ax.plot(x, yp34n, ':k', linewidth=lw)

    
    #p12p, = ax.plot(((-k1+k2)*L+(k1**3-k2**3)*tau+np.log(a*(k4-k2)/(k4-k1)))/(k1**2-k2**2), facc*L, 'xk', linewidth=lw2)
    #p12n,= ax.plot(((-k1+k2)*L+(k1**3-k2**3)*tau+np.log(a*(k3-k2)/(k3-k1)))/(k1**2-k2**2), facc*L, 'xk', linewidth=lw2)
    #p34p, = ax.plot((-k4**2*L+k4**3*tau-np.log(2.0))/k4, facc*L, 'ok', linewidth=lw2)
    #p34n, = ax.plot((-k4**2*L+k4**3*tau-np.log(2.0))/k4, facc*L, 'ok', linewidth=lw2)
    
    plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$u$')
    
    plt.pause(1)
    tau=tau+dt
    
plt.show()
