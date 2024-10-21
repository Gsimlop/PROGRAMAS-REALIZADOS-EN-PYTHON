import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)
    
def recolectador(nombre):
    with open('{}.txt'.format(nombre), 'r') as f:
        # Leer el contenido del archivo
        lines = f.readlines()
        
        # Encontramos en que punto esta alguna informacion interesante
        Encabezado = [i for i, line in enumerate(lines) if "Waveform Data" in line]
        irho = [i for i, line in enumerate(lines) if "Vertical Scale" in line]
        nrho = 0 
        if irho != []:
            partes = lines[irho[0]].split(',') 
            valor, potencia = partes[-2].split('E')
            nrho = float(valor) * pow(10, int(potencia))
        # sacamos los valores de la tabla
        inicio = Encabezado[0] + 1
        primera = []
    for line in lines[inicio:-1]:
        partes = line.split()
        for j in range(len(partes)):
            numero = float(partes[0].replace(',', '.')) 
            if j == 0: 
                primera.append(numero)
    return np.array(primera)*(nrho/25)

##############################
a1 = recolectador('A0001CH1')
a2 = recolectador('A0001CH2')
##############################
a3 = recolectador('A0002CH1')
a4 = recolectador('A0002CH2')
##############################
a5 = recolectador('A0003CH1')
a6 = recolectador('A0003CH2')
##############################
a7 = recolectador('A0004CH1')
a8 = recolectador('A0004CH2')
##############################
a9 = recolectador('A0005CH1')
a10 = recolectador('A0005CH2')
##############################
a11 = recolectador('A0007CH1')
a12 = recolectador('A0007CH2')
##############################
a13 = recolectador('A0008CH1')
a14 = recolectador('A0008CH2')
##############################
a15 = recolectador('A0009CH1')
a16 = recolectador('A0009CH2')
##############################


'Calculo B'

fact1 = (56000*9*10**-6)/(400*3.8*10**-4)
b1 = a2 * fact1
b2 = a4 * fact1
b3 = a6 * fact1
b4 = a8 * fact1
b5 = a10 * fact1
b6 = a12 * fact1
b7 = a14 * fact1
b8 = a16 * fact1

'Calculo H'

fact2 = 200/(29.4*10**-2 *20)
h1 = a1 * fact2
h2 = a3 * fact2
h3 = a5 * fact2
h4 = a7 * fact2
h5 = a9 * fact2
h6 = a11 * fact2
h7 = a13 * fact2
h8 = a15 * fact2

'Calculo maximos B'

maxb1 = np.max(b1)
maxb2 = np.max(b2)
maxb3 = np.max(b3)
maxb4 = np.max(b4)
maxb5 = np.max(b5)
maxb6 = np.max(b6)
maxb7 = np.max(b7)
maxb8 = np.max(b8)

'Calculo maximos H'

maxh1 = np.max(h1)
maxh2 = np.max(h2)
maxh3 = np.max(h3)
maxh4 = np.max(h4)
maxh5 = np.max(h5)
maxh6 = np.max(h6)
maxh7 = np.max(h7)
maxh8 = np.max(h8)

'Calculo B/H'

per1 = maxb1/maxh1
per2 = maxb2/maxh2
per3 = maxb3/maxh3
per4 = maxb4/maxh4
per5 = maxb5/maxh5
per6 = maxb6/maxh6
per7 = maxb7/maxh7
per8 = maxb8/maxh8


def area_encerrada(x_coords, y_coords):
    """
    Calcula el área encerrada por una gráfica dadas dos listas de coordenadas
    (x,y) utilizando el método del trapecio.
    """
    area = 0.0
    n = len(a15)

    for i in range(n-1):
        x1, y1 = a15[i], a16[i]
        x2, y2 = a15[i+1], a16[i+1]

        # Área del trapecio
        a = (x2 - x1) * (y1 + y2) / 2.0
        if a > 0:
            area += a
        else:
            area -= a

    return abs(area)


print("El área encerrada por la gráfica es:", "{0: .2f}". format(area_encerrada(a15, a16)),"m^2")

print("La potencia potencia consumida es:" , "{0: .2f}". format(area_encerrada(a15, a16)* 111.72*10**-6 * 50),"J/m^3")

'Graficamos V'   

fig1 = plt.figure(num=None, figsize=(9, 6), dpi=160, facecolor='w', edgecolor='k')
ax1 = fig1.add_subplot(111)
ax1.set_title('V_B frente V_H',fontsize=20)
##############################
ax1.plot(a1, a2,'r')
##############################
ax1.plot(a3, a4,'k')
##############################
ax1.plot(a5, a6,'tan')
##############################
ax1.plot(a7, a8,'y')
##############################
ax1.plot(a9, a10,'g')
##############################
ax1.plot(a11, a12,'c')
##############################
ax1.plot(a13, a14,'b')
##############################
ax1.plot(a15, a16,'pink')
##############################

ax1.set_xlabel('VH',fontsize=15)
ax1.set_ylabel('VB',fontsize=15)

ax1.grid()
fig1.gca().xaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
fig1.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
#ax1.legend(fontsize= 10, facecolor= 'w',edgecolor='k')
plt.show()


'Graficamos B-H'  

fig2 = plt.figure(num=None, figsize=(9, 6), dpi=160, facecolor='w', edgecolor='k')
ax2 = fig2.add_subplot(111)
ax2.set_title('H-B',fontsize=20)
##############################
ax2.plot(h1, b1,'r')
##############################
ax2.plot(h2, b2,'k')
##############################
ax2.plot(h3, b3,'tan')
##############################
ax2.plot(h4, b4,'y')
##############################
ax2.plot(h5, b5,'g')
##############################
ax2.plot(h6, b6,'c')
##############################
ax2.plot(h7, b7,'b')
##############################
ax2.plot(h8, b8,'pink')
##############################

ax2.set_xlabel('H',fontsize=15)
ax2.set_ylabel('B',fontsize=15)

ax2.grid()
fig2.gca().xaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
fig2.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
#ax2.legend(fontsize= 10, facecolor= 'w',edgecolor='k')
plt.show()

'Graficamos B/H'  

fig3 = plt.figure(num=None, figsize=(9, 6), dpi=160, facecolor='w', edgecolor='k')
ax3 = fig3.add_subplot(111)
ax3.set_title('Curva conmutación',fontsize=20)
##############################
ax3.plot(maxh1, maxb1,'o')
##############################
ax3.plot(maxh2, maxb2,'o')
##############################
ax3.plot(maxh3, maxb3,'o')
##############################
ax3.plot(maxh4, maxb4,'o')
##############################
ax3.plot(maxh5, maxb5,'o')
##############################
ax3.plot(maxh6, maxb6,'o')
##############################
ax3.plot(maxh7, maxb7,'o')
##############################
ax3.plot(maxh8, maxb8,'o')
##############################

ax3.set_xlabel('maxH',fontsize=15)
ax3.set_ylabel('maxB',fontsize=15)

ax3.grid()
fig3.gca().xaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
fig3.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
#ax3.legend(fontsize= 10, facecolor= 'w',edgecolor='k')
plt.show()

'Graficamos B/H frente H'  

fig4 = plt.figure(num=None, figsize=(9, 6), dpi=160, facecolor='w', edgecolor='k')
ax4 = fig4.add_subplot(111)
ax4.set_title('Curva permeabilidad',fontsize=20)
##############################
ax4.plot(maxh1, per1,'o')
##############################
ax4.plot(maxh2, per2,'o')
##############################
ax4.plot(maxh3, per3,'o')
##############################
ax4.plot(maxh4, per4,'o')
##############################
ax4.plot(maxh5, per5,'o')
##############################
ax4.plot(maxh6, per6,'o')
##############################
ax4.plot(maxh7, per7,'o')
##############################
ax4.plot(maxh8, per8,'o')
##############################

ax4.set_xlabel('Hb',fontsize=15)
ax4.set_ylabel('Bm/Hb',fontsize=15)

ax4.grid()
fig4.gca().xaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
fig4.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
#ax4.legend(fontsize= 10, facecolor= 'w',edgecolor='k')
plt.show()
