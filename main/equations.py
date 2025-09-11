import numpy as np

def hill1(x,k,n):
    return (x**n)/(k**n + x**n)

def dgidt(Gn,gin):
    return af*Gn*go - bf*gin

def dadt(an,gin,Lamn,rhon,gc,rhoc,af,bf,n):
    glyco = hill1(gin,gc,1)
    oxphos = (1-Lamn)*hill1(rhon,rhoc,n)*hill1(gin,gc,n)
#    oxphos = hill1(rhon,rhoc,n)*hill1(gin,gc,n)*(1 - hill1(Lamn,1,n))
    return 0.5*af*glyco + 20*af*oxphos - bf*an

def dEdt(En,Lamn,gin,af,bf,gc,n):
    return af*Lamn*hill1(gin,gc,n) - bf*En

def dMdt(Mn,En,aw,Mprod,bw,Ec,n):
    return aw*Mprod - bw*Mn*hill1(En,Ec,n)

def dGdt(Gn,an,aw,bw):
    return 1*aw*(1.5-an) - 0.01*bw*Gn

def dLamdt(Hn,Lamn,aw,bw,hc,n):
    return aw*hill1(Hn,hc,n) - bw*Lamn

def dHdt(an,Hn,aw,bw,ac,n):
    return aw - bw*Hn*hill1(an,ac,n)
