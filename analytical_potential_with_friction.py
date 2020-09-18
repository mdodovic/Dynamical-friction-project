import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#H0 = 2.1923877e-18 #1/s
G = 4.515353945657138e-39 #kpc^3/(Msol*s^2)
Mass = 884603000000.0 #Msol
gustina = 277.72/(0.71*0.71)  #Msol/kpc^2
scale_length = 8.18    #kpc
scale_length_3 = scale_length**3
povrsinska_gustina = 27e4 #bezdimenzioni parametar
#Ro_critical = 3*(H0**2)/(8*math.pi*G)
Rvir = 262.5
log_Rvir_sl= math.log((1+Rvir)/scale_length) - (Rvir / (Rvir + scale_length))
Half_mass_radius = 18.49558899734247 #masa na kojoj je masa haloa pola
CHANDRA_const = -4*math.pi*G*G
pi_sqrt = math.sqrt(math.pi)
Pot_const = -4*math.pi*G*povrsinska_gustina*gustina*(scale_length**3)
def potencijal(r):
    FI = (Pot_const/r)*(math.log((r+scale_length)/scale_length))
    return FI

def gustina_nfw(r):
    Ro = povrsinska_gustina*gustina/((r/scale_length)*((1 + r/scale_length)**2))
    return Ro

def masa(r):
    M =4.0*math.pi*gustina*povrsinska_gustina*scale_length_3*(math.log((r+scale_length)/ scale_length) - (r/(r+scale_length)))
    return M

def f(r):
    fh = -1.0 * ((G * masa(Rvir))/(r**3)) * (math.log((1+r)/scale_length) - (r/(r + scale_length))) / log_Rvir_sl*(r - 0)
    return fh

def sigma(rast):
    donja_granica = rast #polozaj cestice 
    gornja_granica = 10000. #beskonacno velika donja granica integrala
    #r = np.linspace(donja_granica, gornja_granica, gornja_granica, True) #pavljenje niza rastojanja od polozaja cestice do velike vrednosti 1Mpc
    korak = 1. #1kpc    
    x = []
    y = []
    r_prim = donja_granica
    while r_prim <= gornja_granica:
        f1 = gustina_nfw(r_prim)
        f2 = f(r_prim)
        x.append(r_prim)
        y.append(f1*f2)
        r_prim += korak
    k = -sp.integrate.simps(y, x)
    y[:] = []
    x[:] = []
        
        
    s = (1/gustina_nfw(rast)) * k
    return s
def xhi(rast, vx, vy, vz):
    v = (vx**2 + vy**2 + vz**2)**0.5
    X = v/(math.sqrt(2)*sigma(rast))
    return X    
    
def half_mass_radius():
    Uk_mass = Mass
    Red_mass = Uk_mass/2.    
    rad = 0
    Krit_rad = -1    
    while rad <= 1000 :
        if masa(rad) >= Red_mass:
            Krit_rad = rad
            break;
        rad += 0.0000001
    return Krit_rad
        
def log_Culon(r):
    bmax = r
    rh = Half_mass_radius
    vtyp = sigma(r) 
    M = Mass
    L = math.log(bmax/max(rh,G*M/(vtyp**2)))
    return L

def chandra(m, rast, vx, vy, vz, prosledjena_brzina_za_koordinatu):
    v = (vx**2 + vy**2 + vz**2)**0.5
    v_3 = v**3
    X = xhi(rast,vx,vy,vz)
    a = ((CHANDRA_const * m * gustina_nfw(rast) * log_Culon(rast))/v_3)*(sp.special.erf(X) - (2*X*(math.e**(-X*X)))/pi_sqrt) * prosledjena_brzina_za_koordinatu
    return a

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

N=1
G = 4.515353945657138e-39 #kpc^3/(Msol*s^2)
Map = 884603000000.0 #Msol
R = 200 #kpc
#vk = (G*Map/R)**(0.5) #kpc/s
#print(vk)
vk = 2.46798*3.24078e-17 #kpc/s
tela = []        
tela.append(telo(R,0,0,0,0,vk,0,0,0,0,0,0,88460300))
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
Kzaplot = []
Pzaplot = []

dt = 86400*365*1e6#milion godina
dt = dt * 50 #za 50miliona godina 
t=0
dn=0
tmax=86400*365*1e9#milijardu godina
DUR = 25 # Koliko milijardi godina ce trajati simulacija
tmax = tmax * DUR # za DUR milijardi godina
#V_krit = (G*1.989e42/(400*3.08*1e22))**0.5
dist_ap_x = 0
dist_ap_y = 0
dist_ap_z = 0

for i in range(N):
    tela[i].vxh = -(tela[i].vx + tela[i].ax * dt / 2)
    tela[i].vyh = -(tela[i].vy + tela[i].ay * dt / 2)
    tela[i].vzh = -(tela[i].vz + tela[i].az * dt / 2)

for i in range(N):
    ek.append(1/2*tela[i].m*(tela[i].vx**2 + tela[i].vy**2 + tela[i].vz**2))
    tela[i].vx += -(tela[i].ax * dt)
    tela[i].vy += -(tela[i].ay * dt)
    tela[i].vz += -(tela[i].az * dt)
#for i in range(N):
#    rap = ((tela[i].x - dist_ap_x)**2 + (tela[i].y - dist_ap_y)**2 + (tela[i].z - dist_ap_z)**2)**(0.5)
#    ep.append(potencijal(rap) * tela[i].m)
for i in range(N):
    tela[i].x += tela[i].vxh * dt
    tela[i].y += tela[i].vyh * dt
    tela[i].z += tela[i].vzh * dt
for i in range(N):
    pom = ek[i] #+ ep[i]
    u.append(pom)
Poc = sum(u)
Uzaplot.append(sum(u))
#Pzaplot.append(sum(ep))
#Kzaplot.append(sum(ek))
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
        #print(rap)                
        Fxap = - G * tela[j].m * masa(rap) * (tela[j].x - dist_ap_x)/((rap)**3)
        Fyap = - G * tela[j].m * masa(rap) * (tela[j].y - dist_ap_y)/(
        (rap)**3)
        Fzap = - G * tela[j].m * masa(rap) * (tela[j].z - dist_ap_z)/((rap)**3) #masa po xyz ili kao radijus
        Fx = Fxu + Fxap
        Fy = Fyu + Fyap
        Fz = Fzu + Fzap                
        tela[j].ax = Fx / tela[j].m + chandra(tela[j].m, rap, tela[j].vx, tela[j].vy, tela[j].vz, tela[j].vx)
        tela[j].ay = Fy / tela[j].m + chandra(tela[j].m, rap, tela[j].vx, tela[j].vy, tela[j].vz, tela[j].vy)
        tela[j].az = Fz / tela[j].m + chandra(tela[j].m, rap, tela[j].vx, tela[j].vy, tela[j].vz, tela[j].vz)
        
    for i in range(N):
        tela[i].vxh += tela[i].ax*dt
        tela[i].vyh += tela[i].ay*dt
        tela[i].vzh += tela[i].az*dt
    
    for i in range(N):
        ek.append(tela[i].m*(tela[i].vx**2 + tela[i].vy**2 + tela[i].vz**2)/2)
    
#    for i in range(N):
#        rap = ((tela[i].x - dist_ap_x)**2 + (tela[i].y - dist_ap_y)**2 + (tela[i].z - dist_ap_z)**2)**(0.5)
#        ep.append(potencijal(rap)* tela[i].m)
    
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
        pom = ek[i]# + ep[i]
        u.append(pom)
        
    #ENERGIJA JE KINETICKA        
        
    #Pzaplot.append(sum(ep))
    #Kzaplot.append(sum(ek))
    Uzaplot.append(sum(u))
    #Uzaplot.append((sum(u)-Poc)/Poc)
    Tzaplot.append(t)
    ep[:]=[]
    ek[:]=[]
    u[:]=[]
    
    t+=dt;
    
 
print(VREME)
print(len(VREME))
print(len(x))
   
file = open("C:/Users/matij/Desktop/Novi_momenti/nemamin-Energija-analiticki_sa_trenjem,t,ek.txt","w")
for i in range(len(VREME)):
    file.write(str(VREME[i]) + " " + str(Uzaplot[i]) + "\n")
file.close
file = open("C:/Users/matij/Desktop/Novi_momenti/nemamin-Osnovna_svojstva_cestica-analiticki_sa_trenjem,t,m,x,v.txt","w")
for i in range(len(x)):
    file.write(str(VREME[i]) + " " + str(884603000) + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + " " + str(vxp[i]) + " " + str(vyp[i]) + " " + str(vzp[i]) + "\n")    
file.close
plt.title('Grafik putanje čestice - analitički potencijal')
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')

plt.plot(x,y)
plt.show()
