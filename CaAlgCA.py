import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

os.chdir(sys.path[0])

K1 = 10**-2.88 # From DOI: 10.1016/j.ejpb.2013.12.017
K2 = 10**-4.36
K3 = 10**-5.69

Ca3 = 1880 
Ca2 = 67

Kalg = 10**-3.5
KalgCa = 4500 # From doi:10.1016 / S0168-3659(03)00098-1

pH = np.linspace(1, 10, 500)
h = 10**-pH 

alphaCA = h**3 / (h**3 + h**2*K1 + h*K1*K2 + K1*K2*K3)
alpha1 = h**2*K1 / (h**3 + h**2*K1 + h*K1*K2 + K1*K2*K3)
alpha2 = h**1*K1*K2 / (h**3 + h**2*K1 + h*K1*K2 + K1*K2*K3)
alpha3 = K1*K2*K3 / (h**3 + h**2*K1 + h*K1*K2 + K1*K2*K3)

alphaAlg = 1- (h / (h+ Kalg))

plt.figure(figsize=(6,5))

plt.plot(pH, alphaCA*100,':', label="H$_3$Cit")
plt.plot(pH, alpha1*100,'-.', label="H$_2$Cit$^{-1}$")
plt.plot(pH, alpha2*100,'--', label="HCit$^{-2}$")
plt.plot(pH, alpha3*100,'-', label="Cit$^{-3}$")
plt.plot(pH, alphaAlg*100,'-', label="Alg$^{-}$")


# plt.yscale('log')
# plt.ylim(0.01, 100)
plt.xlabel('pH')
plt.ylabel('Species Distribution (%)')
plt.legend(frameon=False)
plt.tight_layout()
plt.show()


tCit = 0.1 # mol/L
tCa = 0.0125
tAlg =  0.05

Fcit2 = Ca2*alpha2*tCit / (1 + Ca2*alpha2*tCit + Ca3*alpha3*tCit + KalgCa*(alphaAlg*tAlg)*0.5)
Fcit3 = Ca3*alpha3*tCit / (1 + Ca2*alpha2*tCit + Ca3*alpha3*tCit + KalgCa*(alphaAlg*tAlg)*0.5)
FAlg2 = KalgCa*(alphaAlg*tAlg*0.5) / (1 + Ca2*alpha2*tCit + Ca3*alpha3*tCit  + KalgCa*(alphaAlg*tAlg*0.5))

plt.figure(figsize=(6,5))

plt.plot(pH, Fcit2*100, '-', label="Ca.Cit$^{-2}$")
plt.plot(pH, Fcit3*100, '-', label="Ca.Cit$^{-3}$")
plt.plot(pH, FAlg2*100, '-', label="Ca.Alg$^{-2}$")


# plt.yscale('log')
plt.ylim(0.01, 100)
plt.xlabel('pH')
plt.ylabel('Species Distribution (%)')
plt.legend(frameon=False)
plt.tight_layout()
plt.show()