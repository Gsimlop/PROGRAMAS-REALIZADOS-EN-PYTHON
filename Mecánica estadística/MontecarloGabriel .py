import numpy as np
from matplotlib import pyplot as plt
import random

#Numero de pasos totales de Monte Carlo

istepmax = int(input('Numero de pasos de Monte Carlo = '))


#Diferencia de Energía

de = float(input('Diferencia de energia entre niveles = '))


#Introducimos una semilla para los numeros random

aseed = 30

random.seed(30)

#Nombramos listas para las variables estadisticas
ntemp = []
n2temp = []
tempval = []
nexact = []
n2exact = []

#Start the metropolis loop
T = 0.
dT = 0.1

#Primer apartado de 0 a 20 grados
for itemp in range(1,200):
    nval = np.zeros(istepmax)
    T += dT
    n = 1
    nsum = 0
    nsuccess = 0
    for istep in range(istepmax):
        if n==1: 
            nn=0
        if n==0:
            u = random.random()
            if u <= np.exp(-de/T):
                nn = 1
            else:
                nn = 0
        nval[istep] = nn
        n = nn 
        nsum += n
    nave = nsum/istepmax
    exact = 1/(np.exp(de/T) + 1)
    nexact.append(exact)

    ntemp.append(nave)
    tempval.append(T)
    
X0 = np.ones(itemp) * exact

#Gráfica primer apartado(media)
plt.axes(xlim=(0,20),ylim=(-0.01,0.6)) 
plt.plot(tempval,ntemp,'o')
plt.plot(tempval,X0,'--')
plt.plot(tempval,nexact,'-')
plt.axhline(exact, linestyle='--', color = 'r', label = r'$\overline{n_{teórico}} $ = '+str(np.round(exact,3)))
plt.xlabel('T (kT)', fontsize=10)
plt.ylabel('Promedio de n',fontsize=10)
N = str(istepmax)
plt.title("Sistema de dos niveles: " + 'N = ' + N + '; $\epsilon$ = ' + str(de))
plt.legend(loc = 'best')
plt.show()

#Segundo y tercer apartado de 0 a 200 grados

ntemp = []
n2temp = []
tempval = []
nexact = []
n2exact = []
disperteo = []
disperexp = []

T = 0.
dT = 0.1

for i in range(1,2000):  #La temperatura
    nval = np.zeros(istepmax)
    T += dT
    n = 1
    nsum = 0
    nsuccess = 0
    for istep in range(istepmax):
        if n==1: 
            nn=0
        if n==0:
            u = random.random()
            if u <= np.exp(-de/T):
                nn = 1
            else:
                nn = 0
        nval[istep] = nn
        n = nn 
        nsum += n
    nave = nsum/istepmax
    exact = 1/(np.exp(de/T) + 1)
    
    dispersionteo = exact * (1 - exact)
    dispersionexp = nave * (1- nave)
    
    nexact.append(exact)
    disperteo.append(dispersionteo)
    disperexp.append(dispersionexp)
    
    ntemp.append(nave)
    tempval.append(T)


#Definimos una matriz de 1 para graficas el limite de 1/2

X=np.ones(i)*1/2

#Definimos una matriz de 1 para graficar el limite 1/4

X2 = np.ones(i)*1/4

########################### Graficas ########################################

#Gráfica segundo apartado(media)
plt.axes(xlim=(0,200),ylim=(-0.01,0.6)) 
plt.plot(tempval,ntemp,'o')
plt.plot(tempval,X,'--')
plt.plot(tempval,nexact,'-')
plt.xlabel('T (kT)', fontsize=10)
plt.ylabel('Promedio de n',fontsize=10)
plt.axhline(0.5, linestyle='--', color = 'r', label = r'$\overline{n} $ = '+ str(0.5))
plt.title("Sistema de dos niveles")
plt.legend(loc = 'best')
plt.show()
#Gráfica tercer apartado(dispersión)
plt.axes(xlim=(0,200),ylim=(-0.01,0.6)) 
plt.plot(tempval,disperexp,'o')
plt.plot(tempval,X2,'--')
plt.plot(tempval,disperteo,'-')
plt.xlabel('T (kT)', fontsize=10)
plt.ylabel('Dispersion de n',fontsize=10)
plt.axhline(0.25,linestyle='--', color = 'r' ,label = r'$Dispersion$ = '+ str(0.25))
plt.title("Sistema de dos niveles")
plt.legend(loc = 'best')
plt.show()

#plt.errorbar(tempval,ntemp,yerr=n3temp,fmt='o')
#plt.show()
