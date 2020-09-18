import math as math
import numpy as np
import matplotlib.pyplot as plt

#M, X, Y, Z, VX, VY, VZ = np.loadtxt('galaxy.txt', unpack = True)
#analiticki potencijal
gustina = 277.72/(0.71*0.71)  #Msol/kpc^2
scale_length = 8.18    #kpc
povrsinska_gustina = 27e4 #bezdimenzioni parametar

N=1

def masa(r):
    M = 4*math.pi*gustina*povrsinska_gustina*(scale_length**3)*(math.log((r+scale_length)/scale_length) - (r/(r+scale_length)))
    return M
def potencijal(r):
    Gamma = 4.515353945657138e-39 #kpc^3/(Msol*s^2)
    FI = -4*math.pi*Gamma*povrsinska_gustina*gustina*(scale_length**2)*(scale_length/r)*(math.log((r+scale_length)/scale_length))
    return FI
    
class telo:
    x = 0
    y = 0
    z = 0
    vx = 0
    vxh = 0
    vy = 0
    vyh = 0
    vz = 0
    vzh = 0
    ax = 0
    ay = 0
    az = 0
    m = 0
    
    def __init__(self,x,y,z,vx,vxh,vy,vyh,vz,vzh,ax,ay,az,m):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vxh = vxh
        self.vy = vy
        self.vyh = vyh
        self.vz = vz
        self.vzh = vzh
        self.ax = ax
        self.ay = ay
        self.az = az
        self.m = m

G = 4.515353945657138e-39 #kpc^3/(Msol*s^2)
Map = 884603000000.0 #Msol
R = 200 #kpc
vk = (G*Map/R)**(0.5) #kpc/s
#print(vk)
#vk = 24.6798*3.24078e-17 #kpc/s

tela = []        
tela.append(telo(R,0,0,0,0,vk,0,0,0,0,0,0,884603000))
x = []
y = []
z = []
vxp = []
vyp = []
vzp = []
u = []
ek = []
ep = []
Uzaplot = []
Tzaplot = []

dt = 86400*365*1e6#milion godina
dt = dt * 50 #za 50miliona godina 
t=0
dn = 0
tmax=86400*365*1e9#milijardu godina
tmax = tmax * 112.5 # za 11 milijardi godina
#V_krit = (G*1.989e42/(400*3.08*1e22))**0.5
dist_ap_x = 0
dist_ap_y = 0
dist_ap_z = 0

for i in range(N):
    tela[i].vxh = tela[i].vx + tela[i].ax * dt / 2
    tela[i].vyh = tela[i].vy + tela[i].ay * dt / 2
    tela[i].vzh = tela[i].vz + tela[i].az * dt / 2

for i in range(N):
    ek.append(1/2*tela[i].m*(tela[i].vx**2 + tela[i].vy**2 + tela[i].vz**2))
    tela[i].vx += -(tela[i].ax * dt)
    tela[i].vy += -(tela[i].ay * dt)
    tela[i].vz += -(tela[i].az * dt)
for i in range(N):
    rap = ((tela[i].x - dist_ap_x)**2 + (tela[i].y - dist_ap_y)**2 + (tela[i].z - dist_ap_z)**2)**(0.5)
    ep.append(potencijal(rap) * tela[i].m)
for i in range(N):
    tela[i].x += tela[i].vxh * dt
    tela[i].y += tela[i].vyh * dt
    tela[i].z += tela[i].vzh * dt
for i in range(N):
    pom = ek[i] #+ ep[i]
    u.append(pom)
Poc = sum(u)
Uzaplot.append(sum(u))
#Uzaplot.append((sum(u)-Poc)/Poc)
Tzaplot.append(t)
ep[:]=[]
ek[:]=[]
u[:]=[]
VREME=[]
while t<=tmax:
    #print(t/tmax*100)
    VREME.append(dn)
    dn+=0.050
    for j in range(N):
        Fxu = 0
        Fyu = 0
        Fzu = 0
        for k in range(N):
            if(k != j):
                r = ((tela[j].x - tela[k].x)**2 + (tela[j].y - tela[k].y)**2 + (tela[j].z - tela[k].z)**2)**(0.5)
                Fxu += - G*tela[j].m * tela[k].m *(tela[j].x - tela[k].x) / (r**3)
                Fyu += - G*tela[j].m * tela[k].m *(tela[j].y - tela[k].y) / (r**3) 
                Fzu += - G*tela[j].m * tela[k].m *(tela[j].z - tela[k].z) / (r**3)
        rap = ((tela[j].x - dist_ap_x)**2 + (tela[j].y - dist_ap_y)**2 + (tela[j].z - dist_ap_z)**2)**(0.5)                
        Fxap = - G * tela[j].m * masa(rap) * (tela[j].x - dist_ap_x)/((rap)**3)
        Fyap = - G * tela[j].m * masa(rap) * (tela[j].y - dist_ap_y)/((rap)**3)
        Fzap = - G * tela[j].m * masa(rap) * (tela[j].z - dist_ap_z)/((rap)**3) #masa po xyz ili kao radijus
        Fx = Fxu + Fxap
        Fy = Fyu + Fyap
        Fz = Fzu + Fzap                
        tela[j].ax = Fx / tela[j].m
        tela[j].ay = Fy / tela[j].m
        tela[j].az = Fz / tela[j].m
        
    for i in range(N):
        tela[i].vxh += tela[i].ax*dt
        tela[i].vyh += tela[i].ay*dt
        tela[i].vzh += tela[i].az*dt
    
    for i in range(N):
        ek.append(tela[i].m*(tela[i].vx**2 + tela[i].vy**2 + tela[i].vz**2)/2)
    for i in range(N):
        rap = ((tela[i].x - dist_ap_x)**2 + (tela[i].y - dist_ap_y)**2 + (tela[i].z - dist_ap_z)**2)**(0.5)
        ep.append(potencijal(rap) * tela[i].m)
    for i in range(N):
        tela[i].vx = -(tela[i].vxh + tela[i].ax * dt / 2)
        vxp.append(tela[i].vx)
        tela[i].vy = -(tela[i].vyh + tela[i].ay * dt / 2)
        vyp.append(tela[i].vy)        
        tela[i].vz = -(tela[i].vzh + tela[i].az * dt / 2)
        vzp.append(tela[i].vz)
    for i in range(N):
        tela[i].x += tela[i].vxh * dt; 
        x.append(tela[i].x);
        tela[i].y += tela[i].vyh * dt; 
        y.append(tela[i].y);
        tela[i].z += tela[i].vzh * dt; 
        z.append(tela[i].z);
    for i in range(N):
        pom = ek[i] #+ ep[i]
        u.append(pom)

    Uzaplot.append(sum(u))
    #Uzaplot.append((sum(u)-Poc)/Poc)
    Tzaplot.append(t)
    ep[:]=[]
    ek[:]=[]
    u[:]=[]
    t+=dt;
#
#plt.axis([0, 1.8e7, 1e-22, 1e-20])    
file = open("C:/Users/matij/Desktop/Novi_momenti/25-Energija_Cestica_BEZ_TRENJA,t,ek.txt","w")
for i in range(len(VREME)):
    file.write(str(VREME[i]) + " " + str(Uzaplot[i]) + "\n")

file = open("C:/Users/matij/Desktop/Novi_momenti/25-Cestica_BEZ_TRENJA,t,m,x,v","w")
for i in range(len(x)):
    file.write(str(VREME[i]) + " " + str(884603000) + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + " " + str(vxp[i]) + " " + str(vyp[i]) + " " + str(vzp[i]) + "\n")    

plt.title('Grafik putanje čestice - analitički potencijal')
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.plot(x,y, c= 'g')
plt.show()
