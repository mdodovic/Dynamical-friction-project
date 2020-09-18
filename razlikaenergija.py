import numpy as np
import matplotlib.pyplot as plt
#T3,EK3=np.loadtxt("PUNAEnergija,t,ek.txt",unpack = True)
#N3 =len(T3)
#raz_adt_ndt = []
#t = []
#T,enum,eana,er=np.loadtxt("theenersvekomp.txt",unpack = True)
T1,Ek1=np.loadtxt("25-(ene)Energija-analiticki_sa_trenjem,t,ek.txt",unpack = True)
T2,Ek2=np.loadtxt("25-Energija_Cestica_BEZ_TRENJA,t,ek.txt",unpack = True)
T3,Ek3=np.loadtxt("Energija,ek.txt",unpack = True)
Tana=np.loadtxt("VREME_CESTICE_25---11.txt",unpack = True)
t,ek=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/analiticki_sa_trenjem,ENERGIJA,t,ek.txt",unpack = True)

t1=[]
t2=[]
tr1=[]
t3=[]
tr2=[]
e1=[]
e2=[]
er1=[]
e3=[]
er2=[]
T1 = []
E1 = []
N=len(Tana)
i = 0
energija1 = []
while i < N:
    energija1.append(10*Ek1[i])
    i = i+1

for i in range(len(Tana)):
    e1.append(Ek1[i])    
    t1.append(Tana[i])  
for i in range(len(Tana)):
    e2.append(Ek2[i]/5)    
    t2.append(Tana[i])
for i in range(len(Tana)):
    er1.append(Ek2[i]-Ek1[i])    
    tr1.append(Tana[i])
for i in range(len(Tana)):
    e3.append(Ek3[i])    
    t3.append(Tana[i])
for i in range(len(Tana)):
    er2.append(Ek3[i]-Ek2[i])    
    tr2.append(Tana[i])
for i in range(len(Tana)):
    E1.append(ek[i])
    T1.append(Tana[i])
#e,t =np.loadtxt("treceener.txt",unpack = True) 

#for i in range(N3):
#    if T3[i] <=11:
        #raz_adt_ndt.append(EK3[i]-EK2[i])
        #.append(T3[i])
#for i in range(len(T)):
 #   er[i]*=1e14
#plt.plot(T,er,c = 'b')

#plt.plot(T,eana,c = 'r')
#plt.plot(T1,eana1,c = 'r')
#plt.plot(t,e,c = 'r')

#plt.plot(tr2,er2,c='purple')
#plt.plot(t2,e2,c='g')

#plt.plot(t1,e1,c='b') #analiticki-trenje
#plt.plot(t2,e2,c='g') #numericki
plt.plot(T1,E1)
#plt.title('Grafik razlike kinetičkih energija čestice')
plt.xlabel('T [Gyr]')
plt.ylabel('Ek [10^11 Msol*km^2/s^2]')
plt.show()
"""
plt.show()
file = open("ENERGIJArazlika.txt","w")
for i in range(len(t)):
    file.write(str(t[i]) + " " + str(raz_adt_ndt[i]) + "\n")
"""    