import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

folder=str("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

n_points=40
start_P=40
stop_P=60
P_ideal=np.linspace(start_P,stop_P,n_points)


def sigm(x,d0,n,m):
    return (1-(x/d0)**n)/(1-(x/d0)**m)

def lin(x,a,b):
    return a*x+b

def gauss(x,d0,n,m):
    return (1-(x/d0)**n)/(1-(x/d0)**m)


# Recupere les donnees
P,T,O2 = np.loadtxt(str(folder+"2-CoordinancesO_map.dat"), unpack=True )
O2_3K=O2[T>2500]
P3K=P[T>2500]
T3K=T[T>2500]
O2_2K=O2[T<2500]
P2K=P[T<2500]
T2K=T[T<2500]
O2_25K=O2[T>2000]
P25K=P[T>2000]
T25K=T[T>2000]
O2_25K=O2_25K[T25K<3000]
P25K=P25K[T25K<3000]
T25K=T25K[T25K<3000]

plt.plot(P2K,O2_2K,"g.")
plt.plot(P25K,O2_25K,"b.")
plt.plot(P3K,O2_3K,"r.")
plt.show()

# Recupere les donnees
P,T,C3 = np.loadtxt(str(folder+"3-CoordinancesC_map.dat"), unpack=True )
n_points=50

C3_3K=C3[T>2500]
P3K=P[T>2500]
T3K=T[T>2500]
C3_2K=C3[T<2500]
P2K=P[T<2500]
T2K=T[T<2500]
C3_25K=C3[T>2000]
P25K=P[T>2000]
T25K=T[T>2000]
C3_25K=C3_25K[T25K<3000]
P25K=P25K[T25K<3000]
T25K=T25K[T25K<3000]


plt.plot(P2K,C3_2K,"g.")
plt.plot(P25K,C3_25K,"b.")
plt.plot(P3K,C3_3K,"r.")
plt.show()

# Recupere les donnees
P,T,C4 = np.loadtxt(str(folder+"4-CoordinancesC_map.dat"), unpack=True )
n_points=25

C4_3K=C4[T>2500]
P3K=P[T>2500]
T3K=T[T>2500]
C4_2K=C4[T<2500]
P2K=P[T<2500]
T2K=T[T<2500]
C4_25K=C4[T>2000]
P25K=P[T>2000]
T25K=T[T>2000]
C4_25K=C4_25K[T25K<3000]
P25K=P25K[T25K<3000]
T25K=T25K[T25K<3000]


plt.plot(P2K,C4_2K,"g.")
plt.plot(P25K,C4_25K,"b.")
plt.plot(P3K,C4_3K,"r.")
plt.show()

for i in range(len(P2K)):
    for j in range(i+1,len(P2K)):
        if P2K[i] < P2K[j]:
            stock=P2K[i] 
            P2K[i]=P2K[j]
            P2K[j]=stock
            stock=C4_2K[i]
            C4_2K[i]=C4_2K[j]
            C4_2K[j]=stock
for i in range(len(P25K)):
    for j in range(i+1,len(P25K)):
        if P25K[i] < P25K[j]:
            stock=P25K[i] 
            P25K[i]=P25K[j]
            P25K[j]=stock
            stock=C4_25K[i]
            C4_25K[i]=C4_25K[j]
            C4_25K[j]=stock
for i in range(len(P25K)):
    for j in range(i+1,len(P25K)):
        if P25K[i] < P25K[j]:
            stock=P25K[i] 
            P25K[i]=P25K[j]
            P25K[j]=stock
            stock=C4_25K[i]
            C4_25K[i]=C4_25K[j]
            C4_25K[j]=stock


cut_off=[48,49,52]
a_2K,b_2k = curve_fit(lin,P2K[P2K>cut_off[0]],C4_2K[P2K>cut_off[0]])[0]
a_25K,b_25k = curve_fit(lin,P25K[P25K>cut_off[1]],C4_25K[P25K>cut_off[1]])[0]
a_3K,b_3k = curve_fit(lin,P3K[P3K>cut_off[2]],C4_3K[P3K>cut_off[2]])[0]

plt.plot(P2K,C4_2K,"g.")
plt.plot(P2K[P2K>cut_off[0]],lin(P2K[P2K>cut_off[0]],a_2K,b_2k),"g-")
plt.plot(P25K,C4_25K,"b.")
plt.plot(P25K[P25K>cut_off[1]],lin(P25K[P25K>cut_off[1]],a_25K,b_25k),"b-")
plt.plot(P3K,C4_3K,"r.")
plt.plot(P3K[P3K>cut_off[2]],lin(P3K[P3K>cut_off[2]],a_3K,b_3k),"r-")
plt.xlabel("P (GPa)")
plt.ylabel("T (K)")
plt.show()

file = open(str(folder+"4-CoordinancesC_map_smoothed.dat"),"w") 
for i in range(len(P_ideal)):
    if P_ideal[i] < cut_off[0] :
        file.write("2000 "+str(P_ideal[i])+" "+str(0)+"\n")
    else :
        if lin(P_ideal[i],a_2K,b_2k) > 1:
            file.write("2000 "+str(P_ideal[i])+" "+str(1)+"\n")
        else:
            file.write("2000 "+str(P_ideal[i])+" "+str(lin(P_ideal[i],a_2K,b_2k))+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    if P_ideal[i] < cut_off[1]:
        file.write("2500 "+str(P_ideal[i])+" "+str(0)+"\n")
    else:
        if lin(P_ideal[i],a_25K,b_25k) > 1:
            file.write("2500 "+str(P_ideal[i])+" "+str(1)+"\n")
        else:
            file.write("2500 "+str(P_ideal[i])+" "+str(lin(P_ideal[i],a_25K,b_25k))+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    if P_ideal[i] < cut_off[2]:
        file.write("3000 "+str(P_ideal[i])+" "+str(0)+"\n")
    else:
        if lin(P_ideal[i],a_3K,b_3k) > 1:
            file.write("3000 "+str(P_ideal[i])+" "+str(1)+"\n")
        else:
            file.write("3000 "+str(P_ideal[i])+" "+str(lin(P_ideal[i],a_3K,b_3k))+"\n")
file.write("\n")
file.close()

# Recupere les donnees
T,P,C2 = np.loadtxt(str(folder+"2-CoordinancesC_map.dat"), unpack=True )
n_points=50

C2_2K=C2[T<2500]
P2K=P[T<2500]
T2K=T[T<2500]
C2_25K=C2[T>2000]
P25K=P[T>2000]
T25K=T[T>2000]
C2_25K=C2_25K[T25K<3000]
P25K=P25K[T25K<3000]
C2_3K=C2[T>2500]
P3K=P[T>2500]
T3K=T[T>2500]

for i in range(len(P2K)):
    for j in range(i+1,len(P2K)):
        if P2K[i] < P2K[j]:
            stock=P2K[i] 
            P2K[i]=P2K[j]
            P2K[j]=stock
            stock=C2_2K[i]
            C2_2K[i]=C2_2K[j]
            C2_2K[j]=stock
for i in range(len(P25K)):
    for j in range(i+1,len(P25K)):
        if P25K[i] < P25K[j]:
            stock=P25K[i] 
            P25K[i]=P25K[j]
            P25K[j]=stock
            stock=C2_25K[i]
            C2_25K[i]=C2_25K[j]
            C2_25K[j]=stock
for i in range(len(P3K)):
    for j in range(i+1,len(P3K)):
        if P3K[i] < P3K[j]:
            stock=P3K[i] 
            P3K[i]=P3K[j]
            P3K[j]=stock
            stock=C2_3K[i]
            C2_3K[i]=C2_3K[j]
            C2_3K[j]=stock


d0_2K,n_2K,m_2K = curve_fit(sigm,P2K,C2_2K)[0]
d0_25K,n_25K,m_25K = curve_fit(sigm,P25K,C2_25K,[48,50,60])[0]
d0_3K,n_3K,m_3K = curve_fit(sigm,P3K,C2_3K)[0]

plt.plot(P2K,C2_2K,"g.",P_ideal,sigm(P_ideal,d0_2K,n_2K,m_2K),"g-")
plt.plot(P25K,C2_25K,"r.",P_ideal,sigm(P_ideal,48,50,60),"r-")
plt.plot(P3K,C2_3K,"c.",P_ideal,sigm(P_ideal,d0_3K,n_3K,m_3K),"c-")
plt.xlabel("P (GPa)")
plt.ylabel("T (K)")
plt.show()

file = open(str(folder+"2-CoordinancesC_map_smoothed.dat"),"w") 

for i in range(len(P_ideal)):
    file.write("2000 "+str(P_ideal[i])+" "+str(sigm(P_ideal[i],d0_2K,n_2K,m_2K))+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    file.write("2500 "+str(P_ideal[i])+" "+str(sigm(P_ideal[i],48,50,60))+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    file.write("3000 "+str(P_ideal[i])+" "+str(sigm(P_ideal[i],d0_3K,n_2K,m_3K))+"\n")
file.write("\n")
 
file.close() 

file = open(str(folder+"3-CoordinancesC_map_smoothed.dat"),"w") 
for i in range(len(P_ideal)):
    value=0
    if P_ideal[i] < cut_off[0] :
        value=0
    else :
        if lin(P_ideal[i],a_2K,b_2k) > 1:
            value=1
        else:
            value=lin(P_ideal[i],a_2K,b_2k)
    value=1-value-sigm(P_ideal[i],d0_2K,n_2K,m_2K)
    if value < 0:
        value=0
    file.write("2000 "+str(P_ideal[i])+" "+str(value)+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    value=0
    if P_ideal[i] < cut_off[1] :
        value=0
    else :
        if lin(P_ideal[i],a_25K,b_25k) > 1:
            value=1
        else:
            value=lin(P_ideal[i],a_25K,b_25k)
    value=1-value-sigm(P_ideal[i],48,50,60)
    if value < 0:
        value=0
    file.write("2500 "+str(P_ideal[i])+" "+str(value)+"\n")
file.write("\n")
for i in range(len(P_ideal)):
    value=0
    if P_ideal[i] < cut_off[1] :
        value=0
    else :
        if lin(P_ideal[i],a_3K,b_3k) > 1:
            value=1
        else:
            value=lin(P_ideal[i],a_3K,b_3k)
    value=1-value-sigm(P_ideal[i],d0_3K,n_3K,m_3K)
    if value < 0:
        value=0
    file.write("3000 "+str(P_ideal[i])+" "+str(value)+"\n")
file.close() 

