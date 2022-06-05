from cmath import *
from math import sqrt, pow, pi
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------- #
#                             Função Transformador                             #
# ---------------------------------------------------------------------------- #
def CalcularTransformador(Uf, Rc, k, f, R1, R2, L1, L2, C1, p):
    C2 = C1              # Capacitâncias
    w = 2*pi*f           # Frequência angular
    wr = 1/sqrt(L1*C1)   # Ressonância 
    M = k * sqrt(L1*L2)  # Indutância Mútua
    # -------------------------------- Reatâncias -------------------------------- #
    XL1 = 1j*w*L1        # reatância do Indutor 1
    XL2 = 1j*w*L2        # reatância do Indutor 2
    XC1 = 1/(1j*w*C1)    # reatância da Capacitor 1
    XC2 = 1/(1j*w*C2)    # reatância da Capacitor 2    
    XM = 1j*w*M          # reatância da Indutância Mútua
    # -------------------------------- Prints de Controle -------------------------------- #
    if p:    
      # Reatâncias da Indutância mútua
      print(f"A Indutäncia Mútua M:{M} Henry")
      print(f"A frequência de ressonancia Wr é: {wr/(2*pi)} Hz")
      print(f"A frequência angular W é: {w} Rad/s")
      print(f"Reatancia do indutor: {XL1}")
      print(f"Reatancia do capacitor: {XC1}")
      print(f"Reatancia dos transformadores: {XM}")
      print(f"-"*30)

    # ---------------------- Impedância Equivalente Sistema ---------------------- #
    Z = np.array([[XC1 + R1 + XL1, -XM], [-XM, XL1 + R2 + ( (XC1 * Rc)/(XC1 + Rc) ) ]])
    # ------------------- Impedância Equivalente na Saída (RC) ------------------- #
    Zeq = ( (XC2 *Rc)/(XC2 + Rc) )
    # ------------------------- Tensão em forma fasorial ------------------------- #
    V = np.array([Uf, 0])
    # ------------------- Correntes no primário e no secundário ------------------ #
    i = np.dot(np.linalg.inv(Z),V)
    # ----------------------- Tensão de Saída (Secundário) ----------------------- #
    v2 = Zeq * i[1]
    
    return i[0], i[1], v2

# ---------------------------------------------------------------------------- #
#                                  Parâmetros                                  #
# ---------------------------------------------------------------------------- #
# Entradas da função: Uf(volts), Rc(ohms), k, R1(ohms), L1(Henry), C1(Henry)
Uf = 10   # Tensão da fonte (pico)
Rc = 400     # Resistência do capacitor
k = 0.4714    # Coeficiente de acoplamento
f = 44646.464    # frequencia da fonte
Vt = (2/pi)*Uf*cos(2*pi*f)
R1 = 825e-3    # Resistência da bobina primaria
R2 = 670e-3    # Resistência da bobina secundaria
L1 = 135e-6 #Indutância da bobina primaria
L2 = 101.4e-6  #Indutância da bobina secundaria
C1 = 1e-7   # Capacitância dos capacitores
#rodando com os valores iniciais
i1, i2, v2 = CalcularTransformador(Uf, Rc, k, f, R1, L1, C1, True)
print(f"i1: {i1}")
print(f"i2: {i2}")

# ------------------------------ Tensão de Saída ----------------------------- #
print('Com resistor na saída de %.2f Ω:' %Rc)
print('i1 (pico) = %.6f A' %np.abs(i1))
print('i2 (pico) = %.6f A' %np.abs(i2))
print('V1 (pico) = %.6f V' %Uf)
print('V2 (pico) = %.6f V' %(np.abs(v2)))
print('Relação V1/V2 = %.6f' %(np.abs(Uf/(v2))))



#Construção de gráficos de I1 e I2 em função da frequência:
#Frequencia que vai de 10hz - 100khz
F = np.linspace(1000,50000,100)
I1 = []
I2 = []
print("modulos das correntes em funcão da frequência da fonte de um k = {}".format(k))

#Uf, Rc, k, f, R1, L1, C1, p
max_i1 = 0
max_i2 = 0

for f in F:
    i1, i2, v2 = CalcularTransformador(Uf, Rc, k, f, R1, L1, C1, False)
    mod_i1 = abs(i1)
    mod_i2 = abs(i2)

    I1.append(mod_i1)
    I2.append(mod_i2)

    if mod_i1 > max_i1:
        max_i1 = mod_i1
        f_max_i1 = f

    if mod_i2 > max_i2:
        max_i2 = mod_i2
        f_max_i2 = f


plt.plot(F, I1, label = "corrente I1", color = "b")
plt.plot(F, I2, label = "corrente I2", color = "r")
plt.title("Módulo das correntes em função das correntes")
plt.xlabel("Frequência $f$")
plt.ylabel("Corrente $I_1 I_2$ [A]")
plt.grid(True)
plt.show()

print("Corrente máxima I1 = {} para uma frequencia de {}".format(max_i1, f_max_i1))
print("-"*30)
print("Corrente máxima I2 = {} para uma frequencia de {}".format(max_i2, f_max_i2))



#Construção de gráficos das potências e da eficiência em função da frequência:

# Criando um vetor para variar a frequência de 10hz - 100khz
F = np.linspace(1000,50000,100)
L_s1 = [] 
L_s2 = [] 
print("modulos das correntes em funcão da frequência da fonte de um k = {}".format(k))

max_s1 = 0
max_s2 = 0

f_max_s1 = 0
f_max_s2 = 0


for f in F:
    i1, i2, v2 = CalcularTransformador(Uf, Rc, k, f, R1, L1, C1, False)

    
    i_rc =  v2/Rc

    # Potências:
    S1 = Uf * i1.conjugate()/2
    # print('Potência real = %.2f' %(S1.real))

    S2 = v2 * i_rc.conjugate()/2
    # print('Potência real = %.2f' %(S2.real))

    
    if S1 > max_s1:
        max_s1 = S1
        f_max_s1 = f

    if S2 > max_s2:
        max_s2 = S2
        f_max_s2 = f

    L_s1.append(abs(S1))
    L_s2.append(abs(S2))

plt.plot(F, L_s1, label = "S1")
plt.legend()

plt.plot(F, L_s2, label = "S2")
plt.legend()

plt.xlabel("Frequencia")
plt.ylabel("Potencia")

plt.grid(True)
plt.show()


print("Potencia máxima S1 = {} para uma frequencia de {}".format(max_s1, f_max_s1))
print("-"*30)
print("Potencia máxima S2 = {} para uma frequencia de {}".format(max_s2, f_max_s2))


K=np.linspace(0,1,1000)
V2=[]
f = 36636 # Frequencia fixa, usamos a de pico 
w = 2 * pi * f


for k in K:
    m = k*sqrt(L1*L1)
    xm = 1j*w*m
    i1, i2, v2 = CalcularTransformador(Uf, Rc, k, f, R1, L1, C1, False)
    V2.append(abs(Rc*i2))
    
plt.plot(K,V2)
plt.xlabel("Coeficiente de acoplamento $k$")
plt.ylabel("Tensão $V_2$ [V]")
plt.grid(True)
plt.show()