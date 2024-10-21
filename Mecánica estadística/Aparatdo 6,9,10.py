import numpy as np
import matplotlib.pyplot as plt
#Pedimos numero de iteraciones
I=int(input('Número de iteraciones: '))
#Creamos una lista para el recorrido
recorrido = []
#DEfinimos la Gaussiana
def gaussiana(x,mean,std):
    return (1/(std*np.sqrt(2*np.pi)))*np.exp((-(x-mean)**2)/(2*std**2))
#Definimos el camino          
def CCamino(N):
    salto = np.random.choice(4,N,p=[0.5,0,0.5,0])
    camino = [np.array([0,0])]
    for i in range(N):
        if salto[i] == 1:
            camino.append(np.array([1,0])+camino[i])
        elif salto[i] == 2:
            camino.append(np.array([-1,0])+camino[i])
        elif salto[i] == 3:
            camino.append(np.array([0,1])+camino[i])
        else:
            camino.append(np.array([0,-1])+camino[i])
    return np.array(camino)

#Definimos las iteraciones
def Iteraciones(I,Narray):
    posfx = []
    posfy = []
    for n in Narray:
        recorrido = []
        for i in range(I):
            recorrido.append(CCamino(int(n)))
        temp = np.array(recorrido)
        posfx.append(temp[:,int(n),0])
        posfy.append(temp[:,int(n),1])
    return np.array(posfx), np.array(posfy)
#Y defnimos los valores estadísticos que vamos a usar
def Valores(Posx,Posy):
    mediax = []
    mediay = []
    stdx = []
    stdy = []
    for i in range(len(Posfx)):
        mediax.append(np.mean(Posx[i]))
        mediay.append(np.mean(Posy[i]))
        stdx.append(np.std(Posx[i]))
        stdy.append(np.std(Posy[i]))
    return np.array(mediax),np.array(mediay),np.array(stdx),np.array(stdy)

def f(N):
    return np.sqrt(N/2)

def media(N,p,q):
    return N*(p-q)     

Narray = np.linspace(2,1000,15)

Posfx,Posfy = Iteraciones(I,Narray)

mediax, mediay, stdx, stdy = Valores(Posfx,Posfy)

array = np.linspace(1,1000,500)
#Graficamos 
plt.plot(array,f(array),'b',label='$\sqrt{N/2}$')
plt.plot(Narray,stdx,'p',label='$\sigma_x$')
plt.plot(Narray,stdy,'d',label='$\sigma_y$')
plt.legend()
plt.show()

plt.plot()
plt.plot(array,media(array,0,0.5),'b',label='$\{N(Pder-Pizq)}$')
plt.plot(Narray,mediax,'p',label='$\mu_x$')
plt.legend()
plt.plot()
plt.show()

plt.plot(array,media(array,0,0.5),'b',label='$\{N(Pder-Pizq)}$')
plt.plot(Narray,mediay,'p',label='$\mu_y$')
plt.legend()
plt.plot()
plt.show()










