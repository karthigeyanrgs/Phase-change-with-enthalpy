import fenics
import matplotlib.pyplot as plt
import numpy as np
import fenics as fs
def semiphasefieldwithenthalpy(N):
    phi=[]
    h=np.empty(N)
    phi=np.linspace(0,1,30)
    phi[0:N]=0
    h[0:N/3]=np.linspace(-10,0,N/3)
    phi[(2*N)/3:N]=1
    h[(2*N)/3:N]=np.linspace(6,16,N/3)
    phi[N/3:(2*N)/3]=np.linspace(0.01,0.99,N/3)
    h[N/3:(2*N)/3]=np.linspace(0,6,N/3)
    plt.plot(h,phi)
    plt.ylim([0,2])
    plt.show()
    return [phi,h]
