########### PRÁCTICA 1. POZO INFINITO CON PARED MÓVIL ###########

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numba import jit
from mpl_toolkits.mplot3d import Axes3D

#Definimos las constantes

n=50
dz=1/n
dt=0.05
tmax=50
m = 2 

#Definimos las variables, que están relacionadas con la pared

a = 1    #Longitud pared instante inicial (t=0)
b = 2    # Velocidad pared
w = 2    # Frecuencia angular pared

#Definimos vector espacial y temporal

z = np.arange(0, 1+dz, dz)      #espacial
t = np.arange(0, tmax+dt, dt)   #temporal 

#Creamos una funcion que describe el movimiento de la pared.

@jit
def l(t):
    return a+b*t # Lineal 
    #return a+b*np.sin(w*t) # Periódico

#Y su derivada.
    
@jit
def lprima(t):
    return b # Lineal 
    #return w*b*np.cos(w*t) # Periódico

#Creamos el entorno para representar los casos

Z,T=np.meshgrid(z,t)

# Deshacemos el cambio de variable con X = Z*l
X = np.zeros((len(t), len(z)))
for i in range(len(t)):
    for j in range(len(z)):
        X[i,j] = Z[i,j]*l(t[i])


######## 1.1 APROXIMACIÓN ADIABÁTICA ########

#suponemos que lprima << l
aproxad = np.zeros((len(t),len(z))) #array de ceros que llenaremos con valores de espacio y tiempo

for i in range(len(t)):
    for j in range(len(z)):
        aproxad[i, j] = 2/l(t[i])*(np.sin(m*np.pi*z[j])**2)  

#REPRESENTACIÓN GRÁFICA
    
fig, adb = plt.subplots(1,1)
adb = plt.axes(projection='3d')
adb.set_xlabel('Espacio')
adb.set_ylabel('Tiempo')
adb.set_zlabel('Densidad de probabilidad')
adb.set_title('Densidad de probabilidad (Aproximación adiabática)')
adb.plot_surface(X, T,aproxad, rstride=1, cstride=1, edgecolor='none')
adb.set_ylim((t[-1],t[0]))
adb.set_zlim((0,6))

####### 1.2 RESOLUCIÓN NUMÉRICA #######

@jit
def phi0(z):
    return np.sqrt(2/l(0))*np.sin(m*np.pi*z)

phi = np.zeros((len(z)), dtype = np.complex64) #añadimos el tipo para tener complejos
phi[:] = phi0(z) #llenamos el array con la función inicial

@jit
def funcion(t,phi):
    phi[0] = 0.+0j 
    phi[-1] = 0.+0j
    retf = [np.complex64(x) for x in range(0)]
    for k in range(len(z)):
        if k == 0 or k == len(z)-1:
            retf += [0.+0j]
        else:
            retf += [(1j/(dz**2 * l(t)**2)) * (phi[k+1] - 2*phi[k] + phi[k-1])+ (1/(2*dz))*z[k]*(lprima(t)/l(t)) * (phi[k+1] - phi[k-1])]
    return np.array(retf, dtype = np.complex64)

# Resolucion por Runge-Kutta

solu = solve_ivp(fun = funcion, t_span = (t[0], t[-1]), t_eval = t, y0 = phi,method = 'RK45', rtol = 1e-5, atol= 1e-5)

        
#REPRESENTACIÓN GRÁFICA
fig, num = plt.subplots(1,1)
num = plt.axes(projection='3d')
num.set_xlabel('Espacio')
num.set_ylabel('Tiempo')
num.set_zlabel('Densidad de probabilidad')
num.set_title('Densidad de probabilidad (Runge-Kutta)')
num.plot_surface(X, T, np.abs(np.transpose(solu.y))**2, rstride=1,cstride=1,edgecolor='none')
num.set_ylim((t[-1],t[0]))
num.set_zlim((0,6))

