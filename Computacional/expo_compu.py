import numpy as np
import matplotlib.pyplot as plt


###############################################################################
#                               FUNCIONES
###############################################################################



def LennardJones(r):
    sigma = 2.3151; eps = 0.167; n = 12; m = 6
    return 4*eps*( (sigma/r)**n - 2*(sigma/r)**m )
#4*eps*( (sigma/r)**n - (sigma/r)**m )
def MorsePotential(r, De=0.167 , a=1, re=2.3151):
    return De*((np.exp(-2*a*(r-re)))- 2 * np.exp(-a*(re-r)))

#De*((np.exp(1-r/re))**12 - 2*(np.exp(1-r/re))**6)
#De*((np.exp(-2*a(r-re)))-2*np.exp(-a*(re-r)))
def CalcR(r):
    n = len(r)
    R = np.empty((n, n))
    for i in range(n):
        R[i] = np.sqrt( (r[:,0]-r[:,0][i])**2 + (r[:,1]-r[:,1][i])**2 + (r[:,2]-r[:,2][i])**2)
    return R

def CalcVM(r, De, a, re, factor_de_corte=3):
    R = CalcR(r)
    n = len(R)
    R_filtered = np.logical_and(R < factor_de_corte * re, R > 0)
    V_red = np.empty(n)
    for i in range(n):
        V_red = 0.5 * np.sum(MorsePotential(R[R_filtered], De, a, re))
    return np.sum(V_red)


def CalcV(r, sigma = 2.3151, factor_de_corte = 3, ConversionRate = 16):
    R = CalcR(r)
    n = len(R)
    R_filter = R<factor_de_corte*sigma
    R_filter2 = 0<R
    R_filtered = np.logical_and(R < factor_de_corte * sigma, R > 0) 
    V_red = np.empty(n)
    for i in range(n):
        #V_red[i] = 0.5*np.sum(LennardJones(R[i][R_filter[i]*R_filter2[i]]))
        V_red = 0.5 * np.sum(LennardJones(R[R_filtered]))
    return np.sum(V_red)

def random_pick_vect(N, dim):
    v = np.random.random_sample(dim)-0.5
    v = v/np.sqrt(np.sum(v**2))
    return  np.random.randint(0,N), v

def random_move(r, step, dim):
    global L
    while True:
        i, vect = random_pick_vect(len(r), dim)
        if (r[i]+vect*step < [L,L,L]).any() or (r[i]+vect*step > [0,0,0]).any():
             r[i] = r[i]+vect*step
             return r

def Metropolis(r, T, E_i):
    global kb, dim, step
    rj = random_move(r.copy(), step, dim)
    E_j = CalcV(rj)
    if E_j <= E_i:
        Pr = 1
    else:
        Pr = np.exp(-(E_j-E_i)/(kb*T))
    if np.random.choice((True, False), p = [Pr, 1-Pr]):
        return rj, E_j
    return r, E_i

def MetropolisMorse(r, T, E_i, De, a, re):
    global kb, dim, step
    rj = random_move(r.copy(), step, dim)
    E_j = CalcVM(rj, De, a, re)
    if E_j <= E_i:
        Pr = 1
    else:
        Pr = np.exp(-(E_j - E_i) / (kb * T))
    if np.random.choice((True, False), p=[Pr, 1 - Pr]):
        return rj, E_j
    return r, E_i




###############################################################################
#                               PARÁMETROS
###############################################################################
kb = 8.6181024e-5
dim = 3
a = 3.603
step = a/100
L = 2*a
sigma=2.3151
eps=0.167
n=12
m=6


###############################################################################
#                               INPUTS
###############################################################################
L = a*int(float(input('Lado del cubo de simulación:   ')))
n = int(float(input('Defina el número de partículas:   ')))
T = float(input('Defina la temperatura del sistema (Kelvin):    '))
n_p = int(float(input('Defina el número de pasos de montecarlo necesarios:   ')))
r = a*np.random.random_sample((n,dim))

###############################################################################
#                               GRÁFICAS
###############################################################################

Energy = np.empty(n_p+1)
Energy[0] = CalcV(r)
EnergyMor = np.empty(n_p+1)
EnergyMor[0] = CalcVM(r, De = 2.3151 , a= 1 , re = 2.3151)


fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Posición inicial.')
ax.scatter(r[:,0], r[:,1], r[:,2], c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

marker = int(n_p/20)
print('Progreso de la simulación 1: |', end = '')
for i in range(n_p):
    r, Energy[i+1] = Metropolis(r, T, Energy[i])
    if i%marker == 0:
        print('#', end = '')
print('|  Completada')
Energy /= n
print(f'Energía final:  {Energy[-1]}')

fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Posición final')
ax.scatter(r[:,0], r[:,1], r[:,2], c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

plt.figure()
plt.title(f'Metrópolis con {n_p} pasos, {n} partículas y T = {T}')
plt.plot(range(n_p+1), Energy, color = 'deepskyblue', linewidth = 0.7)
plt.xlabel('Nº de pasos')
plt.ylabel('Energía')
plt.show()



marker = int(n_p/20)
print('Progreso de la simulación 2: |', end = '')
for i in range(n_p):
    r, EnergyMor[i+1] = MetropolisMorse(r, T, EnergyMor[i], De = 2.3151 , a= 3.603 , re = 2.3151 )
    if i%marker == 0:
        print('#', end = '')
print('|  Completada')
Energy /= n
print(f'Energía final 2:  {EnergyMor[-1]}')

fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Posición final')
ax.scatter(r[:,0], r[:,1], r[:,2], c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

plt.figure()
plt.title(f'Metrópolis con {n_p} pasos, {n} partículas y T = {T}')
plt.plot(range(n_p+1), EnergyMor, color = 'deepskyblue', linewidth = 0.7)
plt.xlabel('Nº de pasos')
plt.ylabel('Energía')
plt.show()

