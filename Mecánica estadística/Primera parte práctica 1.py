import numpy as np
import matplotlib.pyplot as plt
#Pedimos I y N
N = int(input('Numero de pasos:'))
I = int(input('Numero de iteraciones:'))
#Creamos dos matrices de ceros
x = np.zeros(N+1)
y = np.zeros(N+1)
#Definimos la trayectoria aleatoria de las partículas 
def camino(N):
    for i in range(N):
        salto = np.random.choice(['N','S','E','O'], p = [0.5,0,0.5,0])
        if salto == 'N':
            x[i+1] = x[i] + 1
            y[i+1] = y[i] 
        if salto == 'S':
            x[i+1] = x[i] - 1
            y[i+1] = y[i]
        if salto == 'E':
            x[i+1]= x[i]
            y[i+1] = y[i] + 1
        if salto == 'O':
            x[i+1] = x[i]
            y[i+1] = y[i] - 1    
    return (x,y)

x2 = np.zeros(I)
y2 = np.zeros(I)

for j in range(I):
    x, y = camino(N)
    x2[j] = x[N]
    y2[j] = y[N]
    plt.plot(x,y)
plt.xlabel('Direcion Oeste-Este')
plt.ylabel('Direcion Sur-Norte')
plt.title('Trayectorias' )
plt.show()
#Definimos la Gaussiana, que vendrá determinada por la desviación y la media
def gaussiana(media,desv,x):
    return 1./(np.sqrt(2*np.pi)*desv)*np.exp(-(x-media)**2./(2.*desv**2))
def media(N,p,q):
    return N*(p-q) 

#Desviacion y media en x
mediax = np.mean(x2)
desvx = np.std(x2)
#Desviacion y media en y
mediay = np.mean(y2)
desvy = np.std(y2)
#Creamos dos arrays para graficar
xm = np.linspace(-2.5*desvx+mediax, 2.5*desvx+mediax)
ym = np.linspace(-2.5*desvy+mediay, 2.5*desvy+mediay)
#Graficamos 
plt.hist(x2,bins = 'auto',edgecolor = 'black', linewidth = 1, rwidth = 0.5, normed = True )
plt.plot(xm,gaussiana(mediax,desvx,xm))
plt.text(max(x2)-desvx/2,(1/(desvx*np.sqrt(2*np.pi)))/1.5,'$\mu={0}$ \n $\sigma={1}$'.format(mediax,desvx.round(4)))
plt.show()

plt.hist(y2,bins = 'auto',edgecolor = 'black', linewidth = 1, rwidth = 0.5, normed = True )
plt.plot(ym,gaussiana(mediay,desvy,ym))
plt.text(max(y2)-desvy/2,(1/(desvy*np.sqrt(2*np.pi)))/1.5,'$\mu={0}$ \n $\sigma={1}$'.format(mediay,desvy.round(4)))
plt.show()

plt.figure()
plt.plot








