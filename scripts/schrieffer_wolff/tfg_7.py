"""
Análisis perturbativo de tercer orden con g-TMR.

Este script calcula las correcciones de tercer orden en la transformación
de Schrieffer-Wolff, incluyendo los efectos de g-TMR (g-tensor modulated
resonance) a través de los parámetros w_z4 y dw_z4.

Las expresiones simbólicas se simplifican y exportan a archivos LaTeX
para su inclusión en la memoria del TFG.
"""

from src.core import FloquetSW, generate_ls_indices
from sympy import Symbol, init_printing, latex, Matrix, simplify, collect, expand, cancel, factor, apart
import sympy as sp
import numpy as np
from numpy import conjugate as co
t = Symbol(r't')
Om = Symbol(r'\Omega')
w_z = Symbol(r'\omega_z')
w_z2 = Symbol(r'\omega_{z2}')
w_z4 = Symbol(r'\omega_{z4}')
dw_z = Symbol(r'\delta\omega_z')
dw_z2 = Symbol(r'\delta\omega_{z2}')
dw_z4 = Symbol(r'\delta\omega_{z4}')
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
""" V2 = np.array([
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
    [0, 0, 0, 0, 0, ep4]]) """
V4 = np.array([
    [-dw_z4, 0, co(dG_R4), dG_L4, 0, 0],
    [0, dw_z4, co(dG_L4), dG_R4, 0, 0],
    [dG_R4, dG_L4, w_z4, 0, 0, 0],
    [co(dG_L4), co(dG_R4), 0, -w_z4, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])
V2 = np.array([
    [0, 0,0,0, 0, 0],
    [0, 0,0,0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0,0, 0, 0, 0, 0],
    [0, 0, 0, 0, -ep2, 0],
    [0, 0, 0, 0, 0, ep2]])

# Crear instancia de FloquetSW
sw = FloquetSW(H0, V2, V4, w2, w4, privileged=[4, 5])
    
""" alpha = 1#3
beta = 2#0
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
lsA = []
for gamma in range(6):
    for n2, n4 in lsAn2n4:
        l = [gamma, n2, n4]
        if l != m and l != mp:
            lsA.append(l)
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
lsB = []
for gamma in range(6):
    for n2, n4 in lsBn2n4:
        l = [gamma, n2, n4]
        if l != m and l != mp:
            lsB.append(l)
Q1_P4 = sw.SWF_3(m, mp, lsA, lsB, mppsA, mppsB)
Q1_P4 = Q1_P4.subs({w4:-dw_z+w_z}) 
""" 
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
            
Q1P4 = sw.SWF_3(m, mp, lsA, lsB, mppsA, mppsB)
Q1P4 = Q1P4.subs({w4:-dw_z+w_z})



init_printing()

""" # write to text file
with open('tfg_6.txt', 'w') as f:
    f.write(latex(simplify(np.abs(Q1P4 - Q1_P4))))
    f.write('\n')
    f.write(latex(simplify(np.abs(Q1P4 + Q1_P4)))) """
# write to file Q1P4
Q1P4 = cancel(simplify(Q1P4))
numerator, denominator = Q1P4.as_numer_denom()

# remove terms in the numerator where order(e) + order(U) < 4
newnumerator = 0
for term in numerator.as_ordered_terms():
    if sp.degree(term, gen=e) + sp.degree(term, gen=U) >= 4:
        newnumerator += term
Q1P4 = newnumerator/denominator


with open('tfg_7.txt', 'w') as f:
    f.write(latex((simplify(Q1P4))))