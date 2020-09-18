import numpy as np
import matplotlib.pyplot as plt
M3,X3,Y3,Z3,VX3,VY3,VZ3=np.loadtxt("cestica_komplet,m,x,v.txt",unpack = True)

T2,M2,X2,Y2,Z2,VX2,VY2,VZ2=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/25-Cestica_BEZ_TRENJA,t,m,x,v.txt",unpack = True)

T1,M1,X1,Y1,Z1,VX1,VY1,VZ1=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/25-(put)Osnovna_svojstva_cestica-analiticki_sa_trenjem,t,m,x,v.txt",unpack = True)
Tana=np.loadtxt("VREME_CESTICE_25---11.txt",unpack = True)
N1 =len(M1)
#N2 =len(M2)
r1 = []
r2 = []
r3 = []

t1 = []
t2 =[]
t3 = []

tmax = N1
tmin = 0
dt = 0.05
for i in range(len(Tana)):
    r1.append((X1[i]**2+Y1[i]**2+Z1[i]**2)**0.5)    
    t1.append(Tana[i])
    tmin+=dt
tmin = 0  
for i in range(len(Tana)):
    r2.append((X2[i]**2 + Y2[i]**2 + Z2[i]**2)**0.5)    
    t2.append(Tana[i])
    
for i in range(len(X3)):
    r3.append((X3[i]**2 + Y3[i]**2 + Z3[i]**2)**0.5)    
    t3.append(Tana[i]*2)

plt.plot(X1,Y1,c = 'b')
plt.plot(X3,Y3,c = 'r')
    
#    tmin+=dt
#plt.axis([0, 1.8e7, 1e-22, 1e-20])    

#plt.plot(t1,r1,c = 'b')
#plt.plot(t2,r2,c = 'g')
#plt.plot(t3,r3,c = 'r')

#plt.title('Grafik putanje čestice')
#plt.text(75, 200, 'crvena - numerički model')
#plt.text(75, 150, 'plava - analitički potencijal')
plt.xlabel('T [Gyr]')
plt.ylabel('R [kpc]')
plt.show()
