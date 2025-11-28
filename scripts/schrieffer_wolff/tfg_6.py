"""
Transformación de Schrieffer-Wolff con acoplamientos dependientes del tiempo.

Extensión de tfg_5.py que incluye términos de acoplamiento que dependen
de la frecuencia de excitación (V2 y V4), necesarios para el análisis
de excitación bicromática.

Incluye cálculos simbólicos de los elementos del Hamiltoniano de Floquet
en el espacio extendido con índices (α, n₂, n₄).
"""

from src.core import FloquetSW, generate_ls_indices
from sympy import Symbol, init_printing, latex, Matrix, simplify, collect, expand, cancel, factor, apart
import numpy as np
from numpy import conjugate as co
t = Symbol(r't')
Om = Symbol(r'\Omega')
w_z = Symbol(r'\omega_z')
dw_z = Symbol(r'\delta\omega_z')
U = Symbol(r'U')
e = Symbol(r'\epsilon')
G_L = Symbol(r'\gamma_L')
dG_L4 = Symbol(r'\delta\gamma_{L4}')
dG_L2 = Symbol(r'\delta\gamma_{L2}')
G_R = Symbol(r'\gamma_R')
dG_R4 = Symbol(r'\delta\gamma_{R4}')
dG_R2 = Symbol(r'\delta\gamma_{R2}')
hbar = Symbol(r'\hbar')
w2 = Symbol(r'\omega_2')
w4 = Symbol(r'\omega_4')
ep2 = Symbol(r'\epsilon_{P2}')
ep4 = Symbol(r'\epsilon_{P4}')
H0 = np.array([
    [-dw_z  , 0      , co(G_R)  , G_L      , -co(t) , -co(t)],
    [0      , dw_z   , co(G_L)  , G_R      , t      , t],
    [G_R    , G_L    , w_z     , 0        , -Om    , -Om],
    [co(G_L), co(G_R), 0        , -w_z    , -co(Om), -co(Om)],
    [-t   , co(t), -co(Om), -Om    , U-e      , 0],
    [-t   , co(t), -co(Om), -Om    , 0        , U+e]])


V4 = np.array([
    [0, 0, 0,0, 0, 0],
    [0, 0, 0,0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])
V2 = np.array([
    [0, 0,0,0, 0, 0],
    [0, 0,0,0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0, 0, 0, 0, -ep2, 0],
    [0, 0, 0, 0, 0, ep2]])
V2 = np.array([
    [0, 0, co(dG_R2), dG_L2, 0, 0],
    [0, 0, co(dG_L2), dG_R2, 0, 0],
    [dG_R2, dG_L2, 0, 0, 0, 0],
    [co(dG_L2), co(dG_R2), 0, 0, 0, 0],
    [0, 0, 0, 0, -ep2, 0],
    [0, 0, 0, 0, 0, ep2]])
V4 = np.array([
    [0, 0, co(dG_R4), dG_L4, 0, 0],
    [0, 0, co(dG_L4), dG_R4, 0, 0],
    [dG_R4, dG_L4, 0, 0, 0, 0],
    [co(dG_L4), co(dG_R4), 0, 0, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])

# Crear instancia de FloquetSW con los estados de carga como privilegiados
sw = FloquetSW(H0, V2, V4, w2, w4, privileged=[4, 5])
    
alpha = 3
beta = 0
mp = [alpha, 0, 0]
m = [beta, 0, -1]
mppsA = [
    m
]
lsAn2n4 = [
    [0, -1],
    [1, -1],
    [-1, -1],
    [0, 0],
    [0, -2]
]
lsA = generate_ls_indices(lsAn2n4, m, mp)
mppsB = [
    mp
]
lsBn2n4 = [
    [0, 0],
    [-1, 0],
    [1, 0],
    [0, -1],
    [0, 1]
]
lsB = generate_ls_indices(lsBn2n4, m, mp)



init_printing()
HSFW = sw.SWF_3(m, mp, lsA, lsB, mppsA, mppsB)
HSFW = simplify(HSFW.subs({w4:-dw_z+w_z, dG_R2:0, dG_L2:0}))
print(latex(HSFW))
