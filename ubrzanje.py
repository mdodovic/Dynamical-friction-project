import numpy as np
import matplotlib.pyplot as plt

M,X,Y,Z,VX,VY,VZ=np.loadtxt("cestica_komplet,m,x,v.txt",unpack = True)

#T0,M0,X0,Y0,Z0,VX0,VY0,VZ0= np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/25-Cestica_BEZ_TRENJA,t,m,x,v.txt", unpack=True)#ovde ubacujes halo

t,m,x,y,z,vx,vy,vz=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/analiticki_sa_trenjem-CE-VIDIMO,t,m,x,v.txt",unpack = True)

#T,M,X,Y,Z,VX,VY,VZ= np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/25-(put)Osnovna_svojstva_cestica-analiticki_sa_trenjem,t,m,x,v.txt", unpack=True)#ovde ubacujes halo

n = len(m)
i = 1
aap=[]
tap=[]
dt=0.05
dT=0.05
kmtokpc=3.2408e-17
while i<n-1:
    dV= 1/kmtokpc*((vx[i-1]-vx[i])**2 + (vy[i-1]-vy[i])**2 + (vz[i-1]-vz[i])**2)**0.5
    aap.append(dV/dt)
    tap.append(dt)
    dt += dT
    i += 1




#N = len(M0)
i = 1
a = []
t = []
dt = 0.05
dT = 0.05
kmtokpc=3.2408e-17
#while i < N:
#    V1 = 1/kmtokpc*((VX0[i-1])**2 + (VY0[i-1])**2 + (VZ0[i-1])**2)**0.5
#    V2 = 1/kmtokpc*((VX0[i])**2 + (VY0[i])**2 + (VZ0[i])**2)**0.5    
#    a.append(((V2-V1)/dT))
#    t.append(dt)
#    dt += dT
#    i += 1
N = len(VY)
i = 1
af=[]
tf=[]
dt=0.05
dT=0.05
kmtokpc=3.2408e-17
while i<N-1:
    dV= 1/kmtokpc*((VX[i-1]-VX[i])**2 + (VY[i-1]-VY[i])**2 + (VZ[i-1]-VZ[i])**2)**0.5
    af.append(dV/dt)
    tf.append(dt)
    dt += dT
    i += 1

#file = open("C:/Users/matij/Desktop/Novi_momenti/Ubrzanje_cestice-#analiticki_bez_trenja,t,a.txt","w")
#for i in range(len(a)):
    #file.write(str(t[i]) + " " + str(a[i]) + "\n")
plt.plot(tap,aap)
plt.show()
