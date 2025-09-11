import sys
import numpy as np
import matplotlib.pyplot as plt 


from params import *
from equations import *

def getOxygen(rho0,hyprho0,ts,tstart,tend):
    o2 = rho0 + 0*ts
    o2[ts>tstart] = hyprho0
    o2[ts>tend] = 1*rho0
    return o2

T = 5000.0
ts, dt = np.linspace(0,T,int(T/dt)+1,retstep=True)

go = 1.0
rho0 = 1.0 # normal oxygen level
deltaRho0 = 0.4
hyprho0 = rho0 - deltaRho0
gstd = 0.0 # noise 
rstd = 0.0

tstart = 1000
duration = 500


gi = 0*ts; gi[0] = 0.1
Lam = 0*ts; Lam[0] = 0.5
a = 0*ts; a[0] = 1.0
E = 0*ts; E[0] = 1.0
M = 0*ts; M[0] = 1.0
G = 0*ts; G[0] = 1.0
H = 0*ts; H[0] = 1.0

rho = getOxygen(rho0,hyprho0,ts,tstart,tstart+duration)

for i in range(1,len(ts)):
    gop = go*(1 + gstd*np.random.normal(1))

    gi[i] = G[i-1]*gop # glucose in cell 

    a[i] = a[i-1] + dt* dadt(a[i-1],gi[i-1],Lam[i-1],rho[i-1]*(1 + max(-1,rstd*np.random.normal(1))),gc,rhoc,af,bf,n) # ATP production

    E[i] = E[i-1] + dt*dEdt(E[i-1],Lam[i-1],gi[i-1],af,bf,gc,n) # AGE. age formation by ldha via ros production

    M[i] = M[i-1] + dt*dMdt(M[i-1],E[i-1],aw,Mprod,bw,Ec,n) # ecm

    G[i] = G[i-1] + dt*dGdt(G[i-1],a[i-1],aw,bw) # glut1

    Lam[i] = Lam[i-1] + dt*dLamdt(H[i-1],Lam[i-1],aw,bw,hc,n) # LDHA

    H[i] = H[i-1] + dt*dHdt(a[i-1],H[i-1],aw,bw,ac,n) # hif

    print(i,len(ts),rho0, gi[i],a[i],E[i],M[i],G[i],Lam[i],H[i])


data2save = np.zeros((len(ts),9))
data2save[:,0] = ts
data2save[:,1] = rho
data2save[:,2] = gi
data2save[:,3] = a
data2save[:,4] = E
data2save[:,5] = M
data2save[:,6] = G
data2save[:,7] = Lam
data2save[:,8] = H

np.savetxt('data_'+str(gstd)+'_'+str(rstd)+'.dat',data2save)

f, ax = plt.subplots(5,1)
ax[0].plot(ts,rho)
ax[1].plot(ts,a)
ax[2].plot(ts,E)
ax[3].plot(ts,M)
ax[4].plot(ts,H)
plt.show()

