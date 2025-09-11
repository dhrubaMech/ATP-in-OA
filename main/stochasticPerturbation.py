import numpy as np
import matplotlib.pyplot as plt 

from params import *
from equations import *

def getOxygen(rho0,hyprho0):
    o2 = rho0 + 0*ts
    o2[ts>5000] = hyprho0
#    o2[ts>7000] = rho0
    o2[ts>5000] = 0.5*(rho0+hyprho0)+0.5*(rho0-hyprho0)*np.cos(2*np.pi*(ts[ts>5000]-5000)/1000)
#    o2[ts>1000] = rho0 + 0.1*np.sin(2*np.pi*ts[ts>1000]/5)
    return o2

T = 5000.0
ts, dt = np.linspace(0,T,int(T/dt)+1,retstep=True)

#hyprho0 = 2.0
rho = getOxygen(rho0,rho0)

gi = 0*ts; gi[0] = 0.1
Lam = 0*ts; Lam[0] = 0.5
a = 0*ts; a[0] = 1.0
E = 0*ts; E[0] = 1.0
M = 0*ts; M[0] = 1.0
G = 0*ts; G[0] = 1.0
H = 0*ts; H[0] = 1.0

gnois = 0.3
rhoNois = 0.3

for i in range(1,len(ts)):
    gi[i] = G[i-1]*go + gnois*np.random.normal(1) # glucose in cell 
    rho[i-1] += rhoNois*np.random.normal(1)

    a[i] = a[i-1] + dt* dadt(a[i-1],gi[i-1],Lam[i-1],rho[i-1],gc,rhoc,af,bf,n) # ATP production

    E[i] = E[i-1] + dt*dEdt(E[i-1],Lam[i-1],gi[i-1],af,bf,gc,n) # AGE. age formation by ldha via ros production

    M[i] = M[i-1] + dt*dMdt(M[i-1],E[i-1],aw,Mprod,bw,Ec,n) # ecm

    G[i] = G[i-1] + dt*dGdt(G[i-1],a[i-1],aw,bw) # glut1

    Lam[i] = Lam[i-1] + dt*dLamdt(H[i-1],Lam[i-1],aw,bw,hc,n) # LDHA

    H[i] = H[i-1] + dt*dHdt(a[i-1],H[i-1],aw,bw,ac,n) # hif

print(rho0, gi[-1],a[-1],E[-1],M[-1],G[-1],Lam[-1],H[-1])
  
f, ax = plt.subplots(8,1,figsize=(10,8))
ax[0].plot(ts,gi);ax[0].set_ylabel('g')
ax[1].plot(ts,a); ax[1].set_ylabel('a')
ax[2].plot(ts,E); ax[2].set_ylabel('E')
ax[3].plot(ts,M); ax[3].set_ylabel('M')
ax[4].plot(ts,G); ax[4].set_ylabel('G')
ax[5].plot(ts,Lam);ax[5].set_ylabel('Lam')
ax[6].plot(ts,H);  ax[6].set_ylabel('H')
ax[7].plot(ts,rho);  ax[7].set_ylabel('O2')
plt.show()
