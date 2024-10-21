#Metropolis Monte Carlo para el modelo de Ising en 2D
#MJC
#Datos de entrada:
# (1) Longitud de la caja en un direccion. Calculo para una caja cubica.
# (2) Temperatura de la simulacion
# (3) Numero total de pasos de Monte Carlo

import numpy as np
from matplotlib import pyplot as plt
import random

################################Apartado 4,5 y 6##############################
#En estos apartados tenemos que responder una preguntas, introduciendo los 
# datos suministrados.

#Leemos los datos de entrada por pantalla.
#Longitud de la caja en una direccion.

length = int(input('Longitud de la caja = ')) 

#Temperatura

T = float(input('Temperatura = '))

#Numero total de pasos de Monte Carlo

istepmax = int(input('Numero de pasos de Monte Carlo = '))

#Inicializamos el generador de numeros aleatorios

aseed = 30
random.seed(30)

#Generamos un vector que guarda las coordenadas de cada uno de los espines, en 2 dimensiones
#Inicializamos este vector.

spin = np.zeros((length,length))

#Generamos un vector que guarda la energia por espin, localizando el espin por sus coordenadas en 2D.
#Inicializamos este vector

eatom = np.zeros((length,length))

#Inicializamos variables para valores promedio de todos los pasos de Monte Carlo
#vmag0 = magnetizacion total
#ntotal = numero total de espines

vmag0 = 0
ntotal = 0

#Initialize listas para la energia total (energy), magnetizacion total (magnetization) y 
#numero de eventos aceptados (success) que se guardan en cada paso de Monte Carlo


energy = []
magnetization = []
success = []

#Generamos la configuracion inicial de espines, asignando aleatoriamente un valor de -1 o de 1
#en una red cuadrada de longitud en cada lado la longitud inicial tomada como parametro de entrada.
#En cada punto de esta red cuadrada tendremos un valor de espin que puede ser 1 o -1.


for i in range(length):
    for j in range(length):
        if random.random() > 0.5:
            spin[i,j] = 1
        else:
            spin[i,j] = -1

#Calculamos la magnetizacion total, como la suma de todos los espines.

        vmag0 += spin[i,j]

#Contamos el numero total de espines generados.

        ntotal+= 1

print('Numero total de espines = ',ntotal)
print('Magnetizacion = ',vmag0)

#Representamos los espines en un grafico X,Y con dos colores distintos
plt.imshow(spin, aspect='auto', interpolation ='none', extent = [0,length ,0,length]);
plt.axis('equal')
plt.title("Situación Inicial; " + 'N = ' + str(istepmax) + '; T = ' + str(int(T)) + '; L = ' + str(length), weight = 'bold')
plt.xlabel('X', fontsize=10)
plt.ylabel('Y',fontsize=10)
plt.ylim((0,length))
plt.xlim((0,length))
plt.show()

#Calculamos la energia total del sistema.

totalenergy = 0.0

#Buscamos los vecinos de cada espin, considerando solo primeros vecinos y
#con condiciones periodicas.
#Este calculo se puede hacer de muchas maneras. La siguiente es solo una de ellas.
#Vecinos son los espines que estan arriba, abajo, hacia la derecha o hacia la izquierda.
#Utilizamos condiciones periodicas: un espin en la posicion length-1 tiene como vecino el del otro extremo (posicion 0)

for i in range(length):
    for j in range(length):
        if i == 0:
            down = length -1
        else:
            down = i - 1
        if i == length-1:
            up = 0
        else:
            up = i + 1
        if j == 0:
            left = length - 1
        else:
            left = i - 1
        if i == length-1:
            right = 0
        else:
            right = i + 1

        eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])

        totalenergy += 0.5*eatom[i,j]

print('Energía inicial total = ',totalenergy)

#Comenzamos el bucle de Monte Carlo. 
#Es lo mismo que hemos hecho hasta ahora pero moviendo uno de los espines cada vez de -1 a 1 o al reves.
#Inicializamos variables
#istep = numero total de pasos
#nsuccess = numero de pasos de Monte Carlo aceptados.
#nfail = numero de pasos de Monte Carlo no aceptados.
#isteplastsucess = 
#Eave = energia total promediada
#vmag = magnetizacion promediada

istep = 0
nsuccess = 0
nfail = 0
isteplastsuccess = 0
Eave = 0
vmag = vmag0

#Comenzamos el bucle de Monte Carlo para los pasos de simulacion que se han leido como dato de entrada.

for istep in range(istepmax):

#Seleccionamos uno de los espines de forma aleatoria, seleccionamos una posicion 
      
    ipick = random.randint(0,length-1)
    jpick = random.randint(0,length-1)

#Calculamos el cambio en energia del sistema debido a este cambio en la energia del espin.
#De esta forma no tenemos que calcular la energia de todo el sistema cada vez que cambiamos un espin.
#Para hacer este cambio tenemos que utilizar de nuevo condiciones periodicas, en el caso de elegir un espin en el borde.

    if ipick==0:
        down = length-1
    else:
        down = ipick-1
    if ipick==length-1:
        up = 0
    else:
        up = ipick+1
    if jpick==0:
        left = length-1
    else:
        left = jpick-1
    if jpick==length-1:
        right = 0
    else:
        right = jpick+1

    deltaenergy = 2.*spin[ipick,jpick]*(spin[up,jpick]+spin[down,jpick] + spin[ipick,left] + spin[ipick,right]) 

#Aceptamos o rechazamos el cambio en el espin considerando el metodo de Metropolis Monte Carlo
#Si el cambio en energia es negativo (la nueva configuracion tiene una energia mas baja) lo aceptamos.
#Si el cambio en energia es positivo, elegimos un numero aleatorio entre 0 y 1 y comparamos con el factor de Boltzmann.
#Si el numero es menor que ese factor aceptamos la configuracion.

    accept = 0
    if deltaenergy<=0:
        accept = 1
    elif random.random()<np.exp(-(deltaenergy/T)):
        accept = 1

    Etotold = totalenergy

    if accept==1:

#Se se acepta la configuracion, se le da la vuelta al espin
#Se calcula la nueva magnetizacion y se actualiza el numero de sucesos aceptados.
 
        vmag = vmag - spin[ipick,jpick]
        spin[ipick,jpick] = -spin[ipick,jpick]
        vmag = vmag + spin[ipick,jpick]
        isteplastsuccess=istep
        nsuccess+=1    

#Se calcula la energia total del sistema

        totalenergy = 0.0
        for i in range(length):
            for j in range(length):
                if i==1:
                    down = length-1
                else:
                    down = i-1
                if i==length-1:
                    up = 1
                else:
                    up = i+1
                if j==1:
                    left = length-1
                else:
                    left = j-1
                if j==length-1:
                    right = 1
                else:
                    right = j+1

                eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])
                totalenergy = totalenergy + 0.5*eatom[i,j]

        Eave += totalenergy
    else:
#Si no se acepta la configuracion nos quedamos con la anterior, actualizando los contadores.
        Eave += Etotold
        nfail += 1
        vmag = vmag

    magnetization.append(vmag) 
    energy.append(totalenergy) 
    success.append(nsuccess)

#    print('Numero total de pasos = ',istep,'Pasos aceptados = ',nsuccess,'Energia total =',totalenergy)


    
#Representamos la energia total en funcion del número de pasos, la magnetización y la configuracion final.
#Todas las representaciones asi como los bucles son muy mejorables en este programa.


plt.plot(energy,'r--')
plt.xlabel('Numero de pasos', fontsize=10, weight = 'bold')
plt.ylabel('Energia total',fontsize=10, weight = 'bold')
plt.title('Energía Total; '  + 'N = ' + str(istepmax) + '; T = ' + str(int(T)) + '; L = ' + str(length), weight = 'bold')
plt.grid()
plt.show()

plt.plot(magnetization,'r--')
plt.xlabel('Numero de pasos', fontsize=10, weight = 'bold')
plt.ylabel('Magnetizacion',fontsize=10, weight = 'bold')
plt.title('Magnetización; ' + 'N = ' + str(istepmax) + '; T = ' + str(int(T)) + '; L = ' + str(length), weight = 'bold')
plt.grid()
plt.show()

plt.plot(success,'r--')
plt.xlabel('Numero de pasos', fontsize=10, weight = 'bold')
plt.ylabel('Numero de pasos aceptados',fontsize=10, weight = 'bold')
plt.title("Modelo de Ising 2D", weight = 'bold')
plt.grid()
plt.show()

plt.axes(xlim=(0,length),ylim=(0,length))
plt.imshow(spin, aspect='auto', interpolation ='none',extent = [0,length,0,length]);
plt.axis('equal')
plt.xlabel('X', fontsize=10)
plt.ylabel('Y',fontsize=10)
plt.ylim((0,length))
plt.xlim((0,length))
plt.title("Situación Final; " + 'N = ' + str(istepmax) + '; T = ' + str(int(T)) + '; L = ' + str(length), weight = 'bold')
plt.show()

print('Pregunta 4: Observamos que a partir de T=1 los espines se alinean completamente, de forma decreciente')
print('Pregunta 5: Observamos que obtenemos los resultados deseados a partir de N = 20000.')


####################################Apartado 7#################################

Temp = np.linspace(0,5,50)
Mmeanlista = []
Emeanlista = []
for T in Temp:
    print(round(T,4))
    istepmax = 10000
    length = 8
    spin = np.zeros((length,length))
    eatom = np.zeros((length,length))
    vmag0 = 0
    ntotal = 0
    energy = []
    magnetization = []
    success = []
    for i in range(length):
        for j in range(length):
            if random.random() > 0.5:
                spin[i,j] = 1
            else:
                spin[i,j] = -1
    
    #Calculamos la magnetizacion total, como la suma de todos los espines.
    
            vmag0 += spin[i,j]
    
    #Contamos el numero total de espines generados.
    
            ntotal+= 1
    #Calculamos la energia total del sistema.
    
    totalenergy = 0.0
    
    #Buscamos los vecinos de cada espin, considerando solo primeros vecinos y
    #con condiciones periodicas.
    #Este calculo se puede hacer de muchas maneras. La siguiente es solo una de ellas.
    #Vecinos son los espines que estan arriba, abajo, hacia la derecha o hacia la izquierda.
    #Utilizamos condiciones periodicas: un espin en la posicion length-1 tiene como vecino el del otro extremo (posicion 0)
    
    for i in range(length):
        for j in range(length):
            if i == 0:
                down = length -1
            else:
                down = i - 1
            if i == length-1:
                up = 0
            else:
                up = i + 1
            if j == 0:
                left = length - 1
            else:
                left = i - 1
            if i == length-1:
                right = 0
            else:
                right = i + 1
    
            eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])
    
            totalenergy += 0.5*eatom[i,j]
    
    
    #Comenzamos el bucle de Monte Carlo. 
    #Es lo mismo que hemos hecho hasta ahora pero moviendo uno de los espines cada vez de -1 a 1 o al reves.
    #Inicializamos variables
    #istep = numero total de pasos
    #nsuccess = numero de pasos de Monte Carlo aceptados.
    #nfail = numero de pasos de Monte Carlo no aceptados.
    #isteplastsucess = 
    #Eave = energia total promediada
    #vmag = magnetizacion promediada
    
    istep = 0
    nsuccess = 0
    nfail = 0
    isteplastsuccess = 0
    Eave = 0
    vmag = vmag0
    
    #Comenzamos el bucle de Monte Carlo para los pasos de simulacion que se han leido como dato de entrada.
    
    for istep in range(istepmax):
    
    #Seleccionamos uno de los espines de forma aleatoria, seleccionamos una posicion 
          
        ipick = random.randint(0,length-1)
        jpick = random.randint(0,length-1)
    
    #Calculamos el cambio en energia del sistema debido a este cambio en la energia del espin.
    #De esta forma no tenemos que calcular la energia de todo el sistema cada vez que cambiamos un espin.
    #Para hacer este cambio tenemos que utilizar de nuevo condiciones periodicas, en el caso de elegir un espin en el borde.
    
        if ipick==0:
            down = length-1
        else:
            down = ipick-1
        if ipick==length-1:
            up = 0
        else:
            up = ipick+1
        if jpick==0:
            left = length-1
        else:
            left = jpick-1
        if jpick==length-1:
            right = 0
        else:
            right = jpick+1
    
        deltaenergy = 2.*spin[ipick,jpick]*(spin[up,jpick]+spin[down,jpick] + spin[ipick,left] + spin[ipick,right]) 
    
    #Aceptamos o rechazamos el cambio en el espin considerando el metodo de Metropolis Monte Carlo
    #Si el cambio en energia es negativo (la nueva configuracion tiene una energia mas baja) lo aceptamos.
    #Si el cambio en energia es positivo, elegimos un numero aleatorio entre 0 y 1 y comparamos con el factor de Boltzmann.
    #Si el numero es menor que ese factor aceptamos la configuracion.
    
        accept = 0
        if deltaenergy<=0:
            accept = 1
        elif random.random()<np.exp(-(deltaenergy/T)):
            accept = 1
    
        Etotold = totalenergy
    
        if accept==1:
    
    #Se se acepta la configuracion, se le da la vuelta al espin
    #Se calcula la nueva magnetizacion y se actualiza el numero de sucesos aceptados.
     
            vmag = vmag - spin[ipick,jpick]
            spin[ipick,jpick] = -spin[ipick,jpick]
            vmag = vmag + spin[ipick,jpick]
            isteplastsuccess=istep
            nsuccess+=1    
    
    #Se calcula la energia total del sistema
    
            totalenergy = 0.0
            for i in range(length):
                for j in range(length):
                    if i==1:
                        down = length-1
                    else:
                        down = i-1
                    if i==length-1:
                        up = 1
                    else:
                        up = i+1
                    if j==1:
                        left = length-1
                    else:
                        left = j-1
                    if j==length-1:
                        right = 1
                    else:
                        right = j+1
    
                    eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])
                    totalenergy = totalenergy + 0.5*eatom[i,j]
    
            Eave += totalenergy
        else:
    #Si no se acepta la configuracion nos quedamos con la anterior, actualizando los contadores.
            Eave += Etotold
            nfail += 1
            vmag = vmag
    
        magnetization.append(vmag) 
        energy.append(totalenergy) 
        success.append(nsuccess)
        
    Ecopia = energy.copy()
    Esli = Ecopia[-1000:]
    Emean = np.mean(Esli)
    Emeanlista.append(Emean)

    
    Mcop = magnetization.copy()
    Msli = Mcop[-1000:]
    Mmean = np.mean(Msli)
    Mmeanlista.append(Mmean)



plt.plot(Temp,Mmeanlista,'o',markersize = 4,label=r'$\overline{M}$')
plt.axvline(1.7, color = 'red', linestyle ='--',label = 'Temperatura crítica')
plt.axhline(-64,linestyle='--',color = 'blue', label= r'$\overline{M_{-}}$ a bajas temperaturas')
plt.axhline(64,linestyle='--',color = 'g', label= r'$\overline{M_{+}}$ a bajas temperaturas')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{M}}$')
plt.title(r'$\mathbf{\overline{M(T)}}$' + ' para los últimos 1000 pasos con N = 10000')
plt.axhline(0,color='black', linestyle='--')
plt.legend(loc = 'best',fontsize = 8)
plt.show()


plt.plot(Temp,Emeanlista,'o',markersize = 4,label=r'$\overline{E}$')
plt.axvline(1.7, color = 'red', linestyle ='--',label = 'Temperatura crítica')
plt.axhline(-128,linestyle='--',color = 'blue', label= r'$\overline{E}$ a bajas temperaturas')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{E}}$')
plt.title(r'$\mathbf{\overline{E(T)}}$' + ' para los últimos 1000 pasos con N = 10000')
plt.legend(loc = 'best')
plt.show()

print('Observando las gráficas podemos ver que cuando la temperatura se va acercando al valor 1.7 aproximadamente, se toman valores muy dispares, hasta estabilizarse de nuevo a partir de la temperatura con valor aproximado 4')
print('Para temperaturas muy bajas, la magnetización media varía entre -64 y 64, mientras que la energía media toma el valor de E = -128')

print('El valor medio de la energía para los 1000 últimos pasos es: ' + str(Emean) + ' ,Y el valor medio de la magnetizacion para los últimos 1000 pasos es:' + str(Mmean))

#################################Apartado 8####################################

Temp = np.linspace(0.1,5,250)
Emeanlista = []
Cvs = []
for T in Temp:
    print(round(T,4))
    istepmax = 10000
    length = 8
    spin = np.zeros((length,length))
    eatom = np.zeros((length,length))
    vmag0 = 0
    ntotal = 0
    energy = []
    magnetization = []
    success = []
    for i in range(length):
        for j in range(length):
            if random.random() > 0.5:
                spin[i,j] = 1
            else:
                spin[i,j] = -1
    
    #Calculamos la magnetizacion total, como la suma de todos los espines.
    
            vmag0 += spin[i,j]
    
    #Contamos el numero total de espines generados.
    
            ntotal+= 1

    #Calculamos la energia total del sistema.
    
    totalenergy = 0.0
    
    #Buscamos los vecinos de cada espin, considerando solo primeros vecinos y
    #con condiciones periodicas.
    #Este calculo se puede hacer de muchas maneras. La siguiente es solo una de ellas.
    #Vecinos son los espines que estan arriba, abajo, hacia la derecha o hacia la izquierda.
    #Utilizamos condiciones periodicas: un espin en la posicion length-1 tiene como vecino el del otro extremo (posicion 0)
    
    for i in range(length):
        for j in range(length):
            if i == 0:
                down = length -1
            else:
                down = i - 1
            if i == length-1:
                up = 0
            else:
                up = i + 1
            if j == 0:
                left = length - 1
            else:
                left = i - 1
            if i == length-1:
                right = 0
            else:
                right = i + 1
    
            eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])
    
            totalenergy += 0.5*eatom[i,j]
    
    
    #Comenzamos el bucle de Monte Carlo. 
    #Es lo mismo que hemos hecho hasta ahora pero moviendo uno de los espines cada vez de -1 a 1 o al reves.
    #Inicializamos variables
    #istep = numero total de pasos
    #nsuccess = numero de pasos de Monte Carlo aceptados.
    #nfail = numero de pasos de Monte Carlo no aceptados.
    #isteplastsucess = 
    #Eave = energia total promediada
    #vmag = magnetizacion promediada
    
    istep = 0
    nsuccess = 0
    nfail = 0
    isteplastsuccess = 0
    Eave = 0
    vmag = vmag0
    
    #Comenzamos el bucle de Monte Carlo para los pasos de simulacion que se han leido como dato de entrada.
    
    for istep in range(istepmax):
    
    #Seleccionamos uno de los espines de forma aleatoria, seleccionamos una posicion 
          
        ipick = random.randint(0,length-1)
        jpick = random.randint(0,length-1)
    
    #Calculamos el cambio en energia del sistema debido a este cambio en la energia del espin.
    #De esta forma no tenemos que calcular la energia de todo el sistema cada vez que cambiamos un espin.
    #Para hacer este cambio tenemos que utilizar de nuevo condiciones periodicas, en el caso de elegir un espin en el borde.
    
        if ipick==0:
            down = length-1
        else:
            down = ipick-1
        if ipick==length-1:
            up = 0
        else:
            up = ipick+1
        if jpick==0:
            left = length-1
        else:
            left = jpick-1
        if jpick==length-1:
            right = 0
        else:
            right = jpick+1
    
        deltaenergy = 2.*spin[ipick,jpick]*(spin[up,jpick]+spin[down,jpick] + spin[ipick,left] + spin[ipick,right]) 
    
    #Aceptamos o rechazamos el cambio en el espin considerando el metodo de Metropolis Monte Carlo
    #Si el cambio en energia es negativo (la nueva configuracion tiene una energia mas baja) lo aceptamos.
    #Si el cambio en energia es positivo, elegimos un numero aleatorio entre 0 y 1 y comparamos con el factor de Boltzmann.
    #Si el numero es menor que ese factor aceptamos la configuracion.
    
        accept = 0
        if deltaenergy<=0:
            accept = 1
        elif random.random()<np.exp(-(deltaenergy/T)):
            accept = 1
    
        Etotold = totalenergy
    
        if accept==1:
    
    #Se se acepta la configuracion, se le da la vuelta al espin
    #Se calcula la nueva magnetizacion y se actualiza el numero de sucesos aceptados.
     
            vmag = vmag - spin[ipick,jpick]
            spin[ipick,jpick] = -spin[ipick,jpick]
            vmag = vmag + spin[ipick,jpick]
            isteplastsuccess=istep
            nsuccess+=1    
    
    #Se calcula la energia total del sistema
    
            totalenergy = 0.0
            for i in range(length):
                for j in range(length):
                    if i==1:
                        down = length-1
                    else:
                        down = i-1
                    if i==length-1:
                        up = 1
                    else:
                        up = i+1
                    if j==1:
                        left = length-1
                    else:
                        left = j-1
                    if j==length-1:
                        right = 1
                    else:
                        right = j+1
    
                    eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])
                    totalenergy = totalenergy + 0.5*eatom[i,j]
    
            Eave += totalenergy
        else:
    #Si no se acepta la configuracion nos quedamos con la anterior, actualizando los contadores.
            Eave += Etotold
            nfail += 1
            vmag = vmag
    
        magnetization.append(vmag) 
        energy.append(totalenergy) 
        success.append(nsuccess)
        
    Ecopia = energy.copy()
    Esli = Ecopia[-1000:]
    Emean = np.mean(Esli)
    Emeanlista.append(Emean)


    Esli2 = Esli.copy()
    for i in range(len(Esli2)):
        Esli2[i] = Esli2[i]**2
    Esli2mean = np.mean(Esli2)
    cv = (Esli2mean - Emean**2)/(T**2)
    Cvs.append(cv)


plt.plot(Temp,Cvs,'o', color = 'b',markersize = 4,label=r'$C_{v}$')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{C_{v}}}$')
plt.title(r'$\mathbf{C_{v}}$' + ' en función de la temperatura')
plt.legend(loc = 'best')
plt.show()    

print('Aproximadamente, se empieza a ver un pico en el valor medio de la capacidad calorífica cuandot T = 2.5 aproximadamente')

####################################Apartado 9#################################

Temp = np.linspace(0.1,10,250)
Mmeanlista = []
Emeanlista = []
Cvs = []
for T in Temp:
    print(round(T,4))
    istepmax = 10000
    length = 8
    spin = np.zeros((length,length))
    eatom = np.zeros((length,length))
    vmag0 = 0
    ntotal = 0
    energy = []
    magnetization = []
    success = []
    for i in range(length):
        for j in range(length):
            if random.random() > 0.5:
                spin[i,j] = 1
            else:
                spin[i,j] = -1
    
    #Calculamos la magnetizacion total, como la suma de todos los espines.
    
            vmag0 += spin[i,j]
    
    #Contamos el numero total de espines generados.
    
            ntotal+= 1

    #Calculamos la energia total del sistema.
    
    totalenergy = 0.0
    
    #Buscamos los vecinos de cada espin, considerando solo primeros vecinos y
    #con condiciones periodicas.
    #Este calculo se puede hacer de muchas maneras. La siguiente es solo una de ellas.
    #Vecinos son los espines que estan arriba, abajo, hacia la derecha o hacia la izquierda.
    #Utilizamos condiciones periodicas: un espin en la posicion length-1 tiene como vecino el del otro extremo (posicion 0)
    
    for i in range(length):
        for j in range(length):
            if i == 0:
                down = length -1
            else:
                down = i - 1
            if i == length-1:
                up = 0
            else:
                up = i + 1
            if j == 0:
                left = length - 1
            else:
                left = i - 1
            if i == length-1:
                right = 0
            else:
                right = i + 1
    
            eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right]) -spin[i,j]
    
            totalenergy += 0.5*eatom[i,j]
    
    
    #Comenzamos el bucle de Monte Carlo. 
    #Es lo mismo que hemos hecho hasta ahora pero moviendo uno de los espines cada vez de -1 a 1 o al reves.
    #Inicializamos variables
    #istep = numero total de pasos
    #nsuccess = numero de pasos de Monte Carlo aceptados.
    #nfail = numero de pasos de Monte Carlo no aceptados.
    #isteplastsucess = 
    #Eave = energia total promediada
    #vmag = magnetizacion promediada
    
    istep = 0
    nsuccess = 0
    nfail = 0
    isteplastsuccess = 0
    Eave = 0
    vmag = vmag0
    
    #Comenzamos el bucle de Monte Carlo para los pasos de simulacion que se han leido como dato de entrada.
    
    for istep in range(istepmax):
    
    #Seleccionamos uno de los espines de forma aleatoria, seleccionamos una posicion 
          
        ipick = random.randint(0,length-1)
        jpick = random.randint(0,length-1)
    
    #Calculamos el cambio en energia del sistema debido a este cambio en la energia del espin.
    #De esta forma no tenemos que calcular la energia de todo el sistema cada vez que cambiamos un espin.
    #Para hacer este cambio tenemos que utilizar de nuevo condiciones periodicas, en el caso de elegir un espin en el borde.
    
        if ipick==0:
            down = length-1
        else:
            down = ipick-1
        if ipick==length-1:
            up = 0
        else:
            up = ipick+1
        if jpick==0:
            left = length-1
        else:
            left = jpick-1
        if jpick==length-1:
            right = 0
        else:
            right = jpick+1
    
        deltaenergy = 2.*spin[ipick,jpick]*(spin[up,jpick]+spin[down,jpick] + spin[ipick,left] + spin[ipick,right])-spin[ipick,jpick]
    
    #Aceptamos o rechazamos el cambio en el espin considerando el metodo de Metropolis Monte Carlo
    #Si el cambio en energia es negativo (la nueva configuracion tiene una energia mas baja) lo aceptamos.
    #Si el cambio en energia es positivo, elegimos un numero aleatorio entre 0 y 1 y comparamos con el factor de Boltzmann.
    #Si el numero es menor que ese factor aceptamos la configuracion.
    
        accept = 0
        if deltaenergy<=0:
            accept = 1
        elif random.random()<np.exp(-(deltaenergy/T)):
            accept = 1
    
        Etotold = totalenergy
    
        if accept==1:
    
    #Se se acepta la configuracion, se le da la vuelta al espin
    #Se calcula la nueva magnetizacion y se actualiza el numero de sucesos aceptados.
     
            vmag = vmag - spin[ipick,jpick]
            spin[ipick,jpick] = -spin[ipick,jpick]
            vmag = vmag + spin[ipick,jpick]
            isteplastsuccess=istep
            nsuccess+=1    
    
    #Se calcula la energia total del sistema
    
            totalenergy = 0.0
            for i in range(length):
                for j in range(length):
                    if i==1:
                        down = length-1
                    else:
                        down = i-1
                    if i==length-1:
                        up = 1
                    else:
                        up = i+1
                    if j==1:
                        left = length-1
                    else:
                        left = j-1
                    if j==length-1:
                        right = 1
                    else:
                        right = j+1
    
                    eatom[i,j] = -spin[i,j]*(spin[up,j]+spin[down,j]+spin[i,left]+spin[i,right])-spin[i,j]
                    totalenergy = totalenergy + 0.5*eatom[i,j]
    
            Eave += totalenergy
        else:
    #Si no se acepta la configuracion nos quedamos con la anterior, actualizando los contadores.
            Eave += Etotold
            nfail += 1
            vmag = vmag
    
        magnetization.append(vmag) 
        energy.append(totalenergy) 
        success.append(nsuccess)
            
    Ecop = energy.copy()
    Esli = Ecop[-1000:]
    Emean = np.mean(Esli)
    Emeanlista.append(Emean)
    
    
    Mcop = magnetization.copy()
    Msli = Mcop[-1000:]
    Mmean = np.mean(Msli)
    Mmeanlista.append(Mmean)

    Esli2 = Esli.copy()
    for i in range(len(Esli2)):
        Esli2[i] = Esli2[i]**2
    Esli2mean = np.mean(Esli2)
    cv = (Esli2mean - Emean**2)/(T**2)
    Cvs.append(cv)



plt.plot(Temp,Mmeanlista,'o',markersize = 4,label=r'$\overline{M}$')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{M}}$')
plt.title(r'$\mathbf{\overline{M(T)}}$' + ' para los últimos 1000 pasos con N = 10000')
plt.axhline(0,color='red', linestyle='--')
plt.legend(loc = 'best',fontsize = 8)
plt.show()

plt.plot(Temp,Emeanlista,'o',markersize = 4,label=r'$\overline{E}$')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{E}}$')
plt.title(r'$\mathbf{\overline{E(T)}}$' + ' para los últimos 1000 pasos con N = 10000')
plt.legend(loc = 'best')
plt.show()

plt.plot(Temp,Cvs,'o', color = 'blue',markersize = 4,label=r'$C_{v}$')
plt.xlabel('T',weight = 'bold')
plt.ylabel(r'$\mathbf{\overline{C_{v}}}$')
plt.title(r'$\mathbf{C_{v}}$' + ' en función de la temperatura')
plt.legend(loc = 'best')
plt.show()    

print('La magnetización no varía de signo a temperaturas bajas, cambia su valor dependiendo de la magnetizacion del campo')






