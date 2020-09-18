import numpy as np
import matplotlib.pyplot as plt


T1,a1=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/Ubrzanje_cestice-analiticki_sa_trenjem,t,a.txt",unpack = True)

T2,a2=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/Ubrzanje_cestice-analiticki_bez_trenja,t,a.txt",unpack = True)

a3=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/Ubrzanje,a.txt",unpack = True)


Tana=np.loadtxt("C:/Users/matij/Desktop/Novi_momenti/VREME_CESTICE_25---11__ubr.txt",unpack = True)

#plt.plot(T,ar,c = 'b')
#plt.plot(T,anum,c = 'b')
#plt.plot(tana,aana,c = 'r')
#plt.plot(T1,anum1,c = 'r')
#plt.title('Grafik razlike ubrznja čestice')
N = len(Tana)
i = 0
ubr = []
while i < N:
    ubr.append(abs(a1[i]/100))
    i = i+1

i=0
ubr1 = []
while i < N:
    ubr1.append(abs(a3[N-i-1]/100))
    i = i+1


j = 0
ubr2 = []
while j < N:
    ubr2.append(a2[j]/25)
    j = j+1


plt.plot(Tana,ubr,c = 'b')
#plt.plot(Tana,ubr1,c = 'r')
#plt.plot(Tana,ubr1,c = 'g')



plt.xlabel('T [Gyr]')
plt.ylabel('a [10^3 km/s/Gyr]')
plt.show()
#file = open("UBRYANJErazlika.txt","w")
#for i in range(len(t)):
#    file.write(str(t[i]) + " " + str(raz_adt_ndt[i]) + "\n")
