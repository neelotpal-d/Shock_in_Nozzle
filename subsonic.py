# -*- coding: utf-8 -*-
"""
Created on Tue May 28 00:22:30 2019

@author: Neelotpal
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 27 23:14:50 2019

@author: Neelotpal
"""

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



k=np.arange(0.1,0.2,0.05)
vval=[]

L=3
dx=0.05
N=(int(L/dx)+1)
gamma=1.4
R=8.314



denp=np.zeros(N)
denp_temp=np.zeros(N)
pp=np.zeros(N)
pp_temp=np.zeros(N)
vp=np.zeros(N)
Tp=np.zeros(N)
vp_temp=np.zeros(N)
Tp_temp=np.zeros(N)
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
S1=np.zeros(N)
S2=np.zeros(N)
S3=np.zeros(N)

diff_n_U1=np.zeros(N)
diff_np1_U1=np.zeros(N)
diff_n_U2=np.zeros(N)
diff_np1_U2=np.zeros(N)
diff_n_U3=np.zeros(N)
diff_np1_U3=np.zeros(N)
Ap=np.zeros(N)
diff_Apf=np.zeros(N)
diff_Apb=np.zeros(N)
   
for i in range(N):
    denp[i]=1#-0.0003*i*dx
    Tp[i]=1#-0.00009333*i*dx
    vp[i]=0.05#+0.11*i*dx
    
for i in range(N):
    Ap[i]=1+2.2*((i*dx-1.5)**2)
     
for i in range(N):
    diff_Apf[i]=((1+2.2*(((i+1)*dx-1.5)**2))-(Ap[i]))/dx

for i in range(N):
    diff_Apb[i]=((Ap[i])-(1+2.2*(((i-1)*dx-1.5)**2)))/dx
    
    
pend=0.95#denp[N-1]*Tp[N-1]   
    

U1=U1_temp=denp*Ap
U2=U2_temp=U1*vp
U3=U1_temp=denp*((Tp/(gamma-1))+(gamma/2)*vp*vp)*Ap

denp=U1/Ap
pp=denp*Tp

for time in range(200000):
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
    
    for i in range(1,N-1):
        S1[i]=0.2*(abs(pp[i-1]-2*pp[i]+pp[i+1])/(pp[i-1]+2*pp[i]+pp[i+1]))*(U1[i+1]-2*U1[i]+U1[i-1])
        S2[i]=0.2*(abs(pp[i-1]-2*pp[i]+pp[i+1])/(pp[i-1]+2*pp[i]+pp[i+1]))*(U2[i+1]-2*U2[i]+U2[i-1])
        S3[i]=0.2*(abs(pp[i-1]-2*pp[i]+pp[i+1])/(pp[i-1]+2*pp[i]+pp[i+1]))*(U3[i+1]-2*U3[i]+U3[i-1])
    U1_temp=U1+dt*diff_n_U1+S1
    U2_temp=U2+dt*diff_n_U2+S2
    U3_temp=U3+dt*diff_n_U3+S3
     
    U1_temp[N-1]=2*U1_temp[N-2]-U1_temp[N-3]
    U2_temp[0]=2*U2[1]-U2_temp[2]
    U2_temp[N-1]=2*U2_temp[N-2]-U2_temp[N-3]
    vp_temp=U2_temp/U1_temp
    U3_temp[0]=U1_temp[0]*((1/(gamma-1))+(gamma/2)*vp_temp[0]*vp_temp[0]) 
    U3_temp[N-1]=(pend*Ap[N-1]/(gamma-1))+((gamma/2)*U2_temp[N-1]*vp_temp[N-1])    

    Tp_temp=(gamma-1)*((U3_temp/U1_temp)-(gamma/2)*vp_temp*vp_temp)
    denp_temp=U1_temp/Ap
    pp_temp=denp_temp*Tp
    
    F1=U2_temp
    F2=((U2_temp*U2_temp/U1_temp)+((gamma-1)/gamma)*(U3_temp-(gamma/2)*(U2_temp*U2_temp/U1_temp)))
    F3=(gamma*U2_temp*U3_temp/U1_temp)-((gamma*(gamma-1)/2)*(U2_temp*U2_temp*U2_temp/(U1_temp*U1_temp)))
    J2=((gamma-1)/gamma)*(U3_temp-(gamma/2)*(U2_temp*U2_temp/U1_temp))*diff_Apb/Ap

    for i in range(1,N-1):
        diff_np1_U1[i]=-(F1[i]-F1[i-1])/dx
        diff_np1_U2[i]=(-(F2[i]-F2[i-1])/dx)+J2[i]
        diff_np1_U3[i]=-(F3[i]-F3[i-1])/dx
    
    for i in range(1,N-1):
        S1[i]=0.2*(abs(pp_temp[i-1]-2*pp_temp[i]+pp_temp[i+1])/(pp_temp[i-1]+2*pp_temp[i]+pp_temp[i+1]))*(U1_temp[i+1]-2*U1_temp[i]+U1_temp[i-1])
        S2[i]=0.2*(abs(pp_temp[i-1]-2*pp_temp[i]+pp_temp[i+1])/(pp_temp[i-1]+2*pp_temp[i]+pp_temp[i+1]))*(U2_temp[i+1]-2*U2_temp[i]+U2_temp[i-1])
        S3[i]=0.2*(abs(pp_temp[i-1]-2*pp_temp[i]+pp_temp[i+1])/(pp_temp[i-1]+2*pp[i]+pp_temp[i+1]))*(U3_temp[i+1]-2*U3_temp[i]+U3_temp[i-1])    
    
    for i in range(1,N-1):
        U1[i]=U1[i]+0.5*dt*(diff_n_U1[i]+diff_np1_U1[i])+S1[i]
        U2[i]=U2[i]+0.5*dt*(diff_n_U2[i]+diff_np1_U2[i])+S2[i]
        U3[i]=U3[i]+0.5*dt*(diff_n_U3[i]+diff_np1_U3[i])+S3[i]
    
    U1[N-1]=2*U1[N-2]-U1[N-3]
    U2[0]=2*U2[1]-U2[2]
    U2[N-1]=2*U2[N-2]-U2[N-3]
    vp=U2/U1
    U3[0]=U1[0]*((1/(gamma-1))+(gamma/2)*vp[0]*vp[0])
    U3[N-1]=(pend*Ap[N-1]/(gamma-1))+((gamma/2)*U2[N-1]*vp[N-1])
    Tp=(gamma-1)*((U3/U1)-(gamma/2)*vp*vp)
    
#print(denp)
denp=U1/Ap
pp=denp*Tp
#plt.figure(1,figsize=(40,20),dpi=70)
plt.plot(pp)   


plt.show()