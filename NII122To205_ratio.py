import matplotlib.pyplot as plt
import os
from astropy.io import ascii
import numpy as np
import sys

def NII122To205_ratio(ne,te):
    # example NII122To205_ratio(1E3, 1E8)
    te    =   te/1e4 
# ------collisional strength
    omega10 = 0.431 * te ** (0.099 + 0.014 * np.log(te)) #cm^3 s-1 
    omega20 = 0.273 * te ** (0.166 + 0.030 * np.log(te)) # 
    omega21 = 1.15  * te ** (0.137 + 0.024 * np.log(te)) # Appendix Table F.2 in Draine 2005 book 

# ----- partition factors 
    g0  =  1.
    g1  =  2*1. + 1.
    g2  =  2*2. + 1.

# ----- collisional rate coefficients 
    cons   =  8.629e-8                # EQ 2.27 in Draine 2005 book  still do not know where it is from
    k10    =  omega10 * cons / te**0.5 / g1  # EQ 2.27 in Draine 2005 book  
    k20    =  omega20 * cons / te**0.5 / g2  # EQ 2.27 in Draine 2005 book  
    k21    =  omega21 * cons / te**0.5 / g2  # EQ 2.27 in Draine 2005 book  

# ----- Upper Energy levels  in K
    E10eV  =  0.00604            #eV from http://physics.nist.gov/PhysRefData/ASD/lines_form.html 
    E20eV  =  0.01622
    E21eV  =  E20eV-E10eV         
    eV2K   =  1.1605E4
    E10K   =  E10eV * eV2K    #  upper energy level in K          
    E21K   =  E21eV * eV2K    
    E20K   =  E20eV * eV2K    

# ----- Einstein A coefficient 
    A21    =  7.46e-6            #s-1 ;form nist.gov 
    A10    =  2.08e-6            #http://physics.nist.gov/PhysRefData/ASD/lines_form.html 
    A20    =  1.12e-12

    
# ----- excitation rates R_if
    C10 = k10 
    C20 = k20
    C21 = k21
    
    n_gamma_10 = 0
    n_gamma_20 = 0
    n_gamma_21 = 0
    
    R10 = ne  * C10 + A10 * (1 + n_gamma_10)   
    R20 = ne  * C20 + A20 * (1 + n_gamma_20)
    R21 = ne  * C21 + A21 * (1 + n_gamma_21)

    R01 = g1 / g0 * ne  * (C10 * np.exp(-E10K /(te *1e4))   +  A10 * n_gamma_10)    
    R02 = g2 / g0 * ne  * (C20 * np.exp(-E20K /(te *1e4))   +  A20 * n_gamma_20)
    R12 = g2 / g1 * ne  * (C21 * np.exp(-E21K /(te *1e4))   +  A21 * n_gamma_21)

    n0     =  R10 * R20 + R10 * R21 + R21 * R20
    n1Ton0 =  R01 * R20 + R01 * R21 + R21 * R02 
    n2Ton0 =  R02 * R10 + R02 * R12 + R12 * R01 

    f122To205 = n2Ton0 * A21 * E21K / (n1Ton0 * A10 * E10K)
#   print(f122To205)
#   n1toNtot  = n1Ton0 / (n0 + n1Ton0 + n2Ton0)
    return f122To205

#if(n_elements(CIItoNII205) ne 0) then begin
#   NtoH  =   np.array([-4.11,-4.10]) # abundance of N to H  Oberst et al. 2006, 2011 
#   CtoH  =   np.array([-3.4, -3.86]) #  ??? 
#    ;NtoH=7.8e-5 ;;N/H=7.8e-5
#    ;CtoH=1.4e-4 ;;C/H=1.4e-4
#    n1NII=n1Ton0/(n0+n1Ton0+n2Ton0)*10^NtoH[abundance-1]     ;the absolute value is not important
#    A10C=2.29e-6
#    E10KC=91.
#    omega10C=(1.55+1.25*te)/(1.+0.35*te^1.25)
#    ;omega10C=2.0756
#    g0c=2.
#    g1c=4.
#    k10c=omega10c*cons/te^0.5/g1c
#    k01c=k10c*g1c/g0c*exp(-E10KC/te/1e4)
#    n1Ton0C=nee*k01c/(nee*k10c+A10C)
#    n1CII=n1Ton0C/(1.+n1Ton0C)*10^CtoH[abundance-1]
#    CIItoNII205=n1CII*A10C*E10KC/(n1NII*A10*E10K)
#endif




Te = 8E4 # electron temperature in K 
fig, ax = plt.subplots()
x=[]
y=[]
for ii in range(int(1E6)): 
#     ii=ii.
      x.append(ii) 
      y.append(NII122To205_ratio(ii,Te))

ax.plot(x,y, '-.', color='blue', label='$T_e$= 8000 K' )


Te = 1E5
x=[]
y=[]
for ii in range(int(1E6)): 
#     ii=ii.
      x.append(ii) 
      y.append(NII122To205_ratio(ii,Te))


ax.plot(x,y, '-', color='black', label='$T_e$= 10000 K' )

Te = 1.2E5
x=[]
y=[]
for ii in range(int(1E6)): 
#     ii=ii.
      x.append(ii) 
      y.append(NII122To205_ratio(ii,Te))

ax.plot(x,y, '--', color='green', label='$T_e$= 12000 K')



ax.set_xscale('log')
#ax.set_yscale('log')

plt.xlabel('Electron Density ($n_e$) [cm$^{-3}$]')
plt.ylabel('Line Ratio [NII] 122/250')
plt.grid(True)
plt.legend(loc=2, borderaxespad=0.,prop={'size':8})




plt.savefig('NII122_205.pdf')
os.system("open NII122_205.pdf")


