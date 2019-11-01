# -*- coding: utf-8 -*-
"""
Created on Sun May 26 18:37:06 2019

@author: Neelotpal
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:39:36 2019

@author: Neelotpal
"""

import numpy as np
import matplotlib.pyplot as plt



k=np.arange(0)
vval=[]

L=3
dx=0.05
N=(int(L/dx)+1)
gamma=1.4
R=8.314


denp=np.zeros(N)
pp=np.zeros(N)
vp=np.zeros(N)
Tp=np.zeros(N)
ep=Tp[:]

U1=np.zeros(N)
U1_temp=np.zeros(N)
U2=np.zeros(N)
U2_temp=np.zeros(N)
U3=np.zeros(N)
U3_temp=np.zeros(N)
F1=np.zeros(N)
F2=np.zeros(N)
F3=np.zeros(N)
J2=x=np.zeros(N)
diff_n_U1=np.zeros(N)
diff_np1_U1=np.zeros(N)
diff_n_U2=np.zeros(N)
diff_np1_U2=np.zeros(N)
diff_n_U3=np.zeros(N)
diff_np1_U3=np.zeros(N)
Ap=np.zeros(N)
diff_Apf=np.zeros(N)
diff_Apb=np.zeros(N)

for l in [1]:   
    for i in range(N):
        if (i*dx<=0.5):
            denp[i]=1.
            Tp[i]=1.
        elif (i*dx<=1.5):
           denp[i]=1.0-0.366*(i*dx-0.5)#1#
           Tp[i]=1.0-0.167*(i*dx-0.5)#1#
        elif (i*dx<=2.1):
            denp[i]=0.634-0.702*(i*dx-1.5)#1
            Tp[i]=0.833-0.4908*(i*dx-1.5)#1#
        else:
            denp[i]=0.5892+(0.10228)*(i*dx-2.1)#0.6+l+(0.020028)*((N-1)*dx-2.1)#
            Tp[i]=0.93968+0.0622*(i*dx-2.1)
        
    for i in range(N):
        Ap[i]=1+2.2*((i*dx-1.5)**2)
         
    for i in range(N):
        diff_Apf[i]=((1+2.2*(((i+1)*dx-1.5)**2))-(Ap[i]))/dx
    
    for i in range(N):
        diff_Apb[i]=((Ap[i])-(1+2.2*(((i-1)*dx-1.5)**2)))/dx
        
        
    pend=denp[N-1]*Tp[N-1]   
        
    vp=0.59/(denp*Ap)
    U1=U1_temp=denp*Ap
    U2=U2_temp=U1*vp
    U3=U1_temp=denp*((Tp/(gamma-1))+(gamma/2)*vp*vp)*Ap
    
    
    for time in range(20000):
        F1=U2
        F2=((U2*U2/U1)+((gamma-1)/gamma)*(U3-(gamma/2)*(U2*U2/U1)))
        F3=(gamma*U2*U3/U1)-((gamma*(gamma-1)/2)*(U2*U2*U2/(U1*U1)))
        J2=((gamma-1)/gamma)*(U3-(gamma/2)*(U2*U2/U1))*diff_Apf/Ap
        
        for i in range(1,N-1):
            diff_n_U1[i]=-(F1[i+1]-F1[i])/dx
            diff_n_U2[i]=(-(F2[i+1]-F2[i])/dx)+J2[i]
            diff_n_U3[i]=-(F3[i+1]-F3[i])/dx
        dt_all=np.zeros(N-2)
        for i in range(1,N-1):
            dt_all[i-1]=0.5*dx/(vp[i]+(Tp[i]**0.5))
            dt=min(dt_all)
                
        U1_temp=U1+dt*diff_n_U1
        U2_temp=U2+dt*diff_n_U2
        U3_temp=U3+dt*diff_n_U3
    
        F1=U2_temp
        F2=((U2_temp*U2_temp/U1_temp)+((gamma-1)/gamma)*(U3_temp-(gamma/2)*(U2_temp*U2_temp/U1_temp)))
        F3=(gamma*U2_temp*U3_temp/U1_temp)-((gamma*(gamma-1)/2)*(U2_temp*U2_temp*U2_temp/(U1_temp*U1_temp)))
        J2=((gamma-1)/gamma)*(U3_temp-(gamma/2)*(U2_temp*U2_temp/U1_temp))*diff_Apb/Ap
    
        for i in range(1,N-1):
            diff_np1_U1[i]=-(F1[i]-F1[i-1])/dx
            diff_np1_U2[i]=(-(F2[i]-F2[i-1])/dx)+J2[i]
            diff_np1_U3[i]=-(F3[i]-F3[i-1])/dx
            
        for i in range(1,N-1):
            U1[i]=U1[i]+0.5*dt*(diff_n_U1[i]+diff_np1_U1[i])
            U2[i]=U2[i]+0.5*dt*(diff_n_U2[i]+diff_np1_U2[i])
            U3[i]=U3[i]+0.5*dt*(diff_n_U3[i]+diff_np1_U3[i])
        
        U1[N-1]=2*U1[N-2]-U1[N-3]
        U2[0]=2*U2[1]-U2[2]
        U2[N-1]=2*U2[N-2]-U2[N-3]
        vp=U2/U1
        U3[0]=U1[0]*((1/(gamma-1))+(gamma/2)*vp[0]*vp[0])
        
        U3[N-1]=(pend*Ap[N-1]/(gamma-1))+((gamma/2)*U2[N-1]*vp[N-1])
        Tp=(gamma-1)*((U3/U1)-(gamma/2)*vp*vp)
        
   
denp=U1/Ap
pp=denp*Tp
vect=np.arange(0,N*dx,dx)
plt.figure(1,dpi=200)
plt.plot(vect,pp)   
plt.xlabel("x'")
plt.ylabel("p'")

plt.title("Pressure variation without artificial viscosity, Time step=20000, Exit p'= "+str(pend))

plt.show()