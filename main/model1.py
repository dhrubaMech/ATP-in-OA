import numpy as np
import matplotlib.pyplot as plt 

af = 1.0
bf = 1.0

aw = 0.1
bw = 0.1

gc = 1.0


T = 2000.0
dt = 0.05
ts, dt = np.linspace(0,T,int(T/dt)+1,retstep=True)

go = 1.0 # outside glucose
rho0 = 1.0 # O2
hyprho0 = 1e-2  
gc = 1.0
rhoc = 1.0
Ec = 1.0
ac = 0.9
hc = 0.9
n = 4

gi = 0*ts; gi[0] = 0.1
Lam = 0*ts; Lam[0] = 0.5
a = 0*ts; a[0] = 1.0
E = 0*ts; E[0] = 1.0
M = 0*ts; M[0] = 1.0
G = 0*ts; G[0] = 1.0
H = 0*ts; H[0] = 1.0

def getOxygen():
    o2 = rho0 + 0*ts
    o2[ts>1000] = hyprho0
    o2[ts>1500] = rho0
#    o2[ts>1000] = rho0 + 0.1*np.sin(2*np.pi*ts[ts>1000]/5)

    return o2

def hill1(x,k,n):
    return (x**n)/(k**n + x**n)

def dgidt(Gn,gin):
    return af*Gn*go - bf*gin

def dadt(an,gin,Lamn,rhon):
    glyco = hill1(gin,gc,1)
    oxphos = (1-Lamn)*hill1(rhon,rhoc,n)*hill1(gin,gc,1)
    return af*glyco + 20*af*oxphos - bf*an

def dEdt(En,Lamn,gin):
    return af*Lamn*hill1(gin,gc,1) - bf*En

def dMdt(Mn,En):
    return aw - bw*Mn*hill1(En,Ec,1)

def dGdt(Gn,an):
    return aw*(1-hill1(an,ac,n)) - bw*Gn

def dLamdt(Hn,Lamn):
    return aw*hill1(Hn,hc,1) - bw*Lamn

def dHdt(an,Hn):
#    print(aw, bw, Hn, an, bw*Hn*(1 - hill1(an,0*ac,n)))
#    return aw - bw*Hn*(1 - hill1(an,0*ac,n))
    return aw - bw*Hn*(ac**n)/(an**n)

rho = getOxygen()

for i in range(1,len(ts)):
    gi[i] = gi[i-1] + dt*dgidt(G[i-1],gi[i-1])

    a[i] = a[i-1] + dt*dadt(a[i-1],gi[i-1],Lam[i-1],rho[i-1])

    E[i] = E[i-1] + dt*dEdt(E[i-1],Lam[i-1],gi[i-1])

    M[i] = M[i-1] + dt*dMdt(M[i-1],E[i-1])

    G[i] = G[i-1] + dt*dGdt(G[i-1],a[i-1])

    Lam[i] = Lam[i-1] + dt*dLamdt(H[i-1],Lam[i-1])

    H[i] = H[i-1] + dt*dHdt(a[i-1],H[i-1])

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
