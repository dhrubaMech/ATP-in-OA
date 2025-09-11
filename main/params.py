import numpy as np
import matplotlib.pyplot as plt 

af = 1.0
bf = 1.0

aw = 0.1
bw = 0.1


gc = 1.0

dt = 0.05

T = 15000.0
ts, dt = np.linspace(0,T,int(T/dt)+1,retstep=True)

go = 1.0 # outside glucose
rho0 = 1.0 # O2

gc = 1.0
rhoc = 1.0
Ec = 0.4
ac = 1.0
hc = 1.0
Mprod = 0.1   
n = 3
