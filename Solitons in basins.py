#This code visualises a simulation of the three-soliton solution with dimensions 
#generated in a shallow water basin, which is either a square, rectangle or circle. 
#It also shows the shifting x domain and a fixed x domain. Created by Shoaib Mohammed
 
#Code that should not be edited, go to step 1
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation
from math import floor, log10

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



fig = plt.figure(figsize = (8, 9))
gs = fig.add_gridspec(nrows=6,ncols=20)
ax = fig.add_subplot(gs[:-1,:-2])
cbar_ax = fig.add_subplot(gs[1:-2,-1])


'Step 1: Choosing parameters'
delt = 10**(-3) #alter this for more accuracy and phase shift, do not set to <=0
Ahat = 1 #alter this to alter k_1,k_2 etc, cannot be too big otherwise the plot wont work

choice=5 #select a choice for different x,y axis, wavelenghts, depths and amplitudes for eta
if choice==1:
    l=np.sqrt(2) #lamda, i.e. wavelength of wave eta, should be larger than h in meters
    h=np.sqrt(0.25) #depth of the water body in meters
    BasinL=30 #length of basin's wall in meters, ensure that BasinL>BasinW
    BasinW=25 #width of basin's wall in meters
elif choice==2:
    l=1 
    h=np.sqrt(0.5) 
    BasinL=30 
    BasinW=25 
elif choice==3: #used for the tabletop square basin in cm, when computing time divide by 10
    l=13.5
    h=7.5
    BasinL=250
    BasinW=BasinL
elif choice==4: #used for rectangle basin in m, also the dimensions for shifting axis
    l=0.95
    h=0.35
    BasinL=30
    BasinW=27 
elif choice==5: #used for circle basin in m
    BasinL=12.5 
    BasinW=BasinL     
    h=1
    l=1.8



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


'Step 2 choosing time and domain bounds, this is if gen=1 or gen=2 (skip to Step 3 otherwise)'

L =50 #X_0, i.e the dimensionless X domain: [-X_0,X_0]
W=L   #Y_0, i.e the dimensionless Y domain: [-Y_0,Y_0]
tau=-50 #start time (dimensionless)
Ndtau = 11 #no. of changes in time (dimensionless) (this includes tau_0)
dtau = 10 #change in time (dimensionless)

       
#only alter above, go to step 3
tau_0=tau #start time (dimensionless)
tau_1=tau+dtau*(Ndtau-1) #end time (dimensionless)
dx = 0.1
dy = dx

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

lw = 1 # line width

time_text = ax.text(0.3, 1.05, '', transform=ax.transAxes, fontsize='xx-large')

fac = 0.25
levels = np.linspace(0, 9*Abar, 15) 


'Step 3 - Generation choices, choose a gen number from 1-4 and then go to if section, then manually choose which tau you would like to observe'

gen=4 #set to 1 for moving x axis, 2 for fixed x axis, 3 for a close up of a square/rectangular basin, 4 for a circular basin
wframe = None
if gen==1:
    for frame in range(Ndtau): # range(n): number of frames
        t=l*tau/(2**(1/2)*eps*c)# time in seconds
        ax.set_xlim(-(2**(1/2)/3)**(1/3)*l*L+l*tau/(2**(1/2)*eps), ((2**(1/2)/3)**(1/3)*l*L+l*tau/(2**(1/2)*eps)))
        ax.set_ylim(-(eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W), (eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W))
        x_bds = np.arange(-(2**(1/2)/3)**(1/3)*l*L+l*tau/(2**(1/2)*eps), (2**(1/2)/3)**(1/3)*l*L+l*tau/(2**(1/2)*eps), dx)
        y_bds = np.arange((-eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W), eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W,dy)
        x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain

        if wframe:
            wframe.remove()
            cbar_ax.cla()

    
        # Plot the new wireframe and pause briefly before continuing.
        Z = (4/3)**(1/3)*atilde*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l, tau) #eta with shifting domain
        
        wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
    
        time_text.set_text(r'$\eta$(x,y,t=%.2fs)' %t)
    
        plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\eta$  (m)')
    
        plt.pause(1)
        plt.pause(1)
        tau=tau+dtau

    plt.show()
elif gen==2:
    ax.set_xlim(-(2**(1/2)/3)**(1/3)*l*L+l*tau_0/(2**(1/2)*eps), ((2**(1/2)/3)**(1/3)*l*L+l*tau_1/(2**(1/2)*eps)))
    ax.set_ylim(-(eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W), (eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W))
    x_bds = np.arange(-(2**(1/2)/3)**(1/3)*l*L+l*tau_0/(2**(1/2)*eps), (2**(1/2)/3)**(1/3)*l*L+l*tau_1/(2**(1/2)*eps), dx)
    y_bds = np.arange((-eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W), eps**(-1/2)*(2**(1/2)/3)**(2/3)*l*W,dy)
    x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain
    
    for frame in range(Ndtau): # range(n): number of frames
        t=l*tau/(2**(1/2)*eps*c) #time in seconds    
        if wframe:
            wframe.remove()
            cbar_ax.cla()
            l34.remove()
            line12p.remove()
            line56p.remove()
        # Plot the new wireframe and pause briefly before continuing.
        Z = (4/3)**(1/3)*atilde*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l, tau) #eta with fixed domain
        
        wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
    
        #square = patches.Rectangle((-BasinL/2, -BasinW/2), BasinL, BasinW, edgecolor='red', facecolor='none') #dimensions of the basin        
        #rotate = mpl.transforms.Affine2D().rotate_deg(-np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))*180/np.pi) + ax.transData #the appropriate rotation to obtain the ninefold amplitude in the basin
        #square.set_transform(rotate)
        #ax.add_patch(square)
        time_text.set_text(r'$\eta$(x,y,t=%.2fs)' %t)
        l34=ax.axvline(x=(1+np.sqrt(2)*(np.sqrt(2)/3)**(1/3)*eps*k3**2)*c*t)
        yp12p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)-(k6**3-k5**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_146))/((k6**2-k5**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line12p, = ax.plot(x_bds, yp12p, '--k', linewidth=lw)
        yp56p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)+(k5**3-k6**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_245))/((k5**2-k6**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line56p, = ax.plot(x_bds, yp56p, '--k', linewidth=lw)
    
        plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\eta$  (m)')
    
        plt.pause(1)
        plt.pause(1)
        tau=tau+dtau

    plt.show()
elif gen==3:
    if BasinL==BasinW:
        ax.set_xlabel('x (cm)')
        ax.set_ylabel('y (cm)')

        vart=np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))
        taumeet=eps*((k6-k5)**(-1)*np.log(B_246**2/(B_235*B_146))-(6**(1/3)*l**(-1)*BasinL*np.tan(vart)-k3**(-1)*(1-np.tan(vart))*np.log(B_246/B_236)))/(2*eps*(k3**2*(np.tan(vart)-1)+k6**2+k6*k5+k5**2)+6**(1/3)*np.tan(vart))
        tmeet=l*taumeet/(2**(1/2)*eps*c)
        #choose tau below to view [3,4] soliton beginning to generate
        #tau=-3**(1/3)*eps*BasinL/(l*(3**(1/3)+2**(2/3)*k3**2*eps))-4**(1/3)*eps*np.log(B_246/B_236)/(2*k3*(3**(1/3)+4**(1/3)*k3**2*eps)) 
        
        #choose tau below to view [1,2] soliton being generated/ stop generating at the top corner
        tau=eps*((9*2**(1/2))**(1/3)*eps**(1/2)*BasinL*(k5**2-k6**2)+l*np.log(B_246**2/(B_235*B_146)))/(l*(6**(1/3)*(k6-k5)+2*eps*(k6**3-k5**3))) #time for [1,2] and [5,6] soliton to be generated
        
        #choose tau below to view [1,2] and [3,4] soliton meeting at the wall
        #tau=taumeet
        
        Ndtau=11
        dtau=-tau/(Ndtau-3)
        ax.set_xlim(-np.sqrt(2)*BasinL/2, np.sqrt(2)*BasinL/2)
        ax.set_ylim(0,np.sqrt(2)*BasinL)
        x_bds = np.arange(-np.sqrt(2)*BasinL/2, np.sqrt(2)*BasinL/2, dx)
        y_bds = np.arange(0, np.sqrt(2)*BasinL,dy)
        x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain
        
        xmeet=(1+(4/3)**(1/3)*eps*k3**2)*c*tmeet+np.sqrt(2)**(1/3)*l*np.log(B_246/B_236)/(24**(1/3)*k3)
        ymeet= xmeet+ BasinL/np.sqrt(2)
        print('xmeet=%.2f'%xmeet,'ymeet=%.2f'%ymeet,'L_[1,2],[5,6]=%.2f'%abs(np.sqrt(2)*xmeet), 'L_[3,4]=%.2f'%abs(np.sqrt(2)*ymeet))
        
        for frame in range(Ndtau): # range(n): number of frames
            if wframe:
                wframe.remove()
                cbar_ax.cla()
                l34.remove()
                line12p.remove()
                line56p.remove()
            t=l*tau/(2**(1/2)*eps*c)
            tactual=t/10
            # Plot the new wireframe and pause briefly before continuing.
            Z = (4/3)**(1/3)*atilde*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l, tau) 
                
            wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
                
            square = patches.Rectangle((-BasinL/2, -BasinW/2), BasinL, BasinW, edgecolor='red', facecolor='none')#dimensions of the basin        
            rotate = mpl.transforms.Affine2D().rotate_deg(45) + ax.transData #the appropriate rotation to obtain the ninefold amplitude in the basin
            square.set_transform(rotate)
            ax.add_patch(square)
            time_text.set_text(r'$\eta$(x,y,t=%.2f s)' %tactual)
        
        
            l34=ax.axvline(x=(1+np.sqrt(2)*(np.sqrt(2)/3)**(1/3)*eps*k3**2)*c*t+l*np.sqrt(2)**(1/3)*np.log(A_246/A_236)/(24**(1/3)*k3), color='black', linestyle='--', linewidth=lw)
            yp12p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)-(k6**3-k5**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_146))/((k6**2-k5**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
            line12p, = ax.plot(x_bds, yp12p, '--k', linewidth=lw)
            yp56p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)+(k5**3-k6**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_245))/((k5**2-k6**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
            line56p, = ax.plot(x_bds, yp56p, '--k', linewidth=lw)
        
            plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\eta$  (cm)')
        
            plt.pause(2)
            plt.pause(1)
            tau=tau+dtau

        plt.show()
    else:
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
        ax.set_xlim(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1))
        ax.set_ylim(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1))
        x_bds = np.arange(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1), dx)
        y_bds = np.arange(round(-(np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2-1), round((np.sin(vart)*BasinL+np.cos(vart)*BasinW)/2+1),dy)
        x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain
        
        xmeet1=(1+(4/3)**(1/3)*eps*k3**2)*c*tmeet1+np.sqrt(2)**(1/3)*l*np.log(B_246/B_236)/(24**(1/3)*k3)
        ymeet1= xmeet1/np.tan(vart)+ BasinL/(2*np.sin(vart))
        xmeet2=(1+(4/3)**(1/3)*eps*k3**2)*c*tmeet2+np.sqrt(2)**(1/3)*l*np.log(B_246/B_236)/(24**(1/3)*k3)
        ymeet2= -np.tan(vart)*xmeet2-BasinW/(2*np.cos(vart))
        
        W_12=abs((2*xmeet1+BasinL*np.cos(vart)-BasinW*np.sin(vart))/(2*np.sin(vart)))
        W_34=abs((2*xmeet1+BasinL*np.cos(vart)+BasinW*np.sin(vart))/(2*np.sin(vart)))
        L_56=abs((2*xmeet2-BasinL*np.cos(vart)+BasinW*np.sin(vart))/(2*np.cos(vart)))
        S_34=abs((2*(xmeet1*np.cos(vart)+xmeet2*np.sin(vart))+(BasinL*np.cos(vart)+BasinW*np.sin(vart))*(np.cos(vart)+np.sin(vart)))/(np.sin(2*vart)))
        print('xmeet1=%.2f'%xmeet1,'ymeet1=%.2f'%ymeet1, 'xmeet2=%.2f'%xmeet2,'ymeet2=%.2f'%ymeet2,'W_[1,2]=%.2f'%W_12, 'L_[5,6]=%.2f'%L_56, 'S_[3,4]=%.2f'%S_34)
        
        for frame in range(Ndtau): # range(n): number of frames
            if wframe:
                wframe.remove()
                cbar_ax.cla()
                l34.remove()
                line12p.remove()
                line56p.remove()
            t=l*tau/(2**(1/2)*eps*c)
            # Plot the new wireframe and pause briefly before continuing.
            Z = (4/3)**(1/3)*atilde*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l, tau) 
                
            wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
                
            square = patches.Rectangle((-BasinL/2, -BasinW/2), BasinL, BasinW, edgecolor='red', facecolor='none')#dimensions of the basin        
            rotate = mpl.transforms.Affine2D().rotate_deg(-np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))*180/np.pi) + ax.transData #the appropriate rotation to obtain the ninefold amplitude in the basin
            square.set_transform(rotate)
            ax.add_patch(square)
            time_text.set_text(r'$\eta$(x,y,t=%.2fs)' %t)
        
            l34=ax.axvline(x=(1+np.sqrt(2)*(np.sqrt(2)/3)**(1/3)*eps*k3**2)*c*t+l*np.sqrt(2)**(1/3)*np.log(A_246/A_236)/(24**(1/3)*k3), color='black', linestyle='--', linewidth=lw)
            yp12p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)-(k6**3-k5**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_146))/((k6**2-k5**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
            line12p, = ax.plot(x_bds, yp12p, '--k', linewidth=lw)
            yp56p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)+(k5**3-k6**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_245))/((k5**2-k6**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
            line56p, = ax.plot(x_bds, yp56p, '--k', linewidth=lw)
        
            plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\eta$  (m)')
        
            plt.pause(2)
            plt.pause(1)
            tau=tau+dtau

        plt.show()
else: #To represent circle basins, e.g.FloWave
    R=BasinL
    vart=np.arctan((3/np.sqrt(2))**(1/3)*eps**(1/2)*(k5+k6))
    p=4**(1/3)*eps*(k3**2-k6**2-k6*k5-k5**2)/(np.tan(vart)*(3**(1/3)+4**(1/3)*eps*k3**2))
    q=(np.sqrt(2)**(1/3)*l/(3**(1/3)*(k6-k5)*np.tan(vart)))*(np.log(B_246/np.sqrt(B_235*B_146))+(3**(1/3)*(k6-k5)+4**(1/3)*eps*(k6**3-k5**3))*np.log(B_246/B_236)/(2*k3*(3**(1/3)+4**(1/3)*eps*k3**2)))
    omega=-np.arctan(p)-np.arcsin(abs(q)*np.cos(np.arctan(p))/R)
    taumeet=-(3*np.sqrt(8))**(1/3)*eps*R*np.cos(omega)/(l*(3**(1/3)+4**(1/3)*eps*k3**2))-4**(1/3)*eps*np.log(B_246/B_236)/(2*k3*(3**(1/3)+4**(1/3)*eps*k3**2))
    
    #select tau below to view starting time for [3,4] soliton
    tau=-3**(1/3)*np.sqrt(2)*eps*R/(l*(3**(1/3)+4**(1/3)*k3**2*eps))-4**(1/3)*eps*np.log(B_246/B_236)/(2*k3*(3**(1/3)+4**(1/3)*k3**2*eps))
    
    #select tau below to view starting time for [1,2] and [5,6] soliton
    #tau=-(eps*((3*np.sqrt(8))**(1/3)*R*(k6-k5)-4**(1/3)*l*np.cos(vart)*np.log(B_246/(B_235*B_146)**(1/2))))/(l*np.cos(vart)*(3**(1/3)*(k6-k5)+4**(1/3)*eps*(k6**3-k5**3)))
    
    #select tau below to view time [1,2] and [3,4] soliton meet
    #tau=taumeet
    
    #select tau below to view time [1,2] and [5,6] soliton stop generating, warning! plot will go in reverse with this value
    #tau=4**(1/3)*eps*np.log(B_246/np.sqrt(B_235*B_146))/(3**(1/3)*(k6-k5)+4**(1/3)*eps*(k6**3-k5**3))
    
    Ndtau=1
    dtau=-tau/(Ndtau-3)
    ax.set_xlim(-BasinL-1, BasinL+1)
    ax.set_ylim(-BasinL-1, BasinL+1)
    x_bds = np.arange(-BasinL-1, BasinL+1, dx)
    y_bds = np.arange(-BasinL-1, BasinL+1,dy)
    x, y = np.meshgrid(x_bds, y_bds) #dimensioned x and y domain
    
    xmeet=-R*np.cos(omega),
    ymeet=R*np.sin(omega)
    L1256=R*(np.pi/2-omega+vart)
    L34=2*R*omega
    Lpist=L34+2*L1256
    print('xmeet1=%.2f'%xmeet,'ymeet1=%.2f'%ymeet, 'L_[1,2],[5,6]=%.2f'%L1256, 'L_[3,4]=%.2f'%L34, 'L_piston=%.2f'%Lpist)
    
    for frame in range(Ndtau): # range(n): number of frames
        if wframe:
            wframe.remove()
            cbar_ax.cla()
            l34.remove()
            line12p.remove()
            line56p.remove()
        t=l*tau/(2**(1/2)*eps*c)
        # Plot the new wireframe and pause briefly before continuing.
        Z = (4/3)**(1/3)*atilde*generate((3/(2)**(1/2))**(1/3)*((x-c*t)/l), eps**(1/2)*(3/(2)**(1/2))**(2/3)*y/l, tau) 
            
        wframe = ax.contourf(x, y, Z, levels, cmap=cm.twilight, linewidth=0, antialiased=False, alpha=0.8)
            
        circle = patches.Circle((0,0), R, edgecolor='red', facecolor='none')#dimensions of the basin        
        ax.add_patch(circle)
        time_text.set_text(r'$\eta$(x,y,t=%.2fs)' %t)
    
        l34=ax.axvline(x=(1+np.sqrt(2)*(np.sqrt(2)/3)**(1/3)*eps*k3**2)*c*t+l*np.sqrt(2)**(1/3)*np.log(A_246/A_236)/(24**(1/3)*k3), color='black', linestyle='--', linewidth=lw)
        yp12p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)-(k6**3-k5**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_146))/((k6**2-k5**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line12p, = ax.plot(x_bds, yp12p, '--k', linewidth=lw)
        yp56p = ((k6-k5)*(3/np.sqrt(2))**(1/3)*(x_bds-c*t)+(k5**3-k6**3)*eps*c*np.sqrt(2)*t+l*np.log(A_246/A_245))/((k5**2-k6**2)*eps**(1/2)*(3/np.sqrt(2))**(2/3))
        line56p, = ax.plot(x_bds, yp56p, '--k', linewidth=lw)
    
        plt.colorbar(wframe,cax=cbar_ax,orientation='vertical',label='$\eta$  (m)')
    
        plt.pause(2)
        plt.pause(1)
        tau=tau+dtau        
        
    
    
    


















