"""
Cálculos simbólicos adicionales de Schrieffer-Wolff.

Este script calcula correcciones de segundo orden para las diagonales
del Hamiltoniano efectivo (desplazamientos de energía), complementando
los cálculos off-diagonal de los scripts anteriores.

Los términos D1, D2, D10, D20 corresponden a las correcciones de
energía para los estados relevantes en las transiciones Q1 y Q1_.

Estas correcciones son importantes para calcular el desplazamiento
de Bloch-Siegert de forma analítica.
"""

from src.core import FloquetSW, generate_ls_indices
from sympy import Symbol, init_printing, latex, Matrix, simplify, collect, expand, cancel, factor, apart
import sympy as sp
import numpy as np
from numpy import conjugate as co
import matplotlib.pyplot as plt
t = Symbol(r't')
Om = Symbol(r'\Omega')
w_z = Symbol(r'\omega_z', real=True)
w_z2 = Symbol(r'\omega_{z2}', real=True)
w_z4 = Symbol(r'\omega_{z4}', real=True)
dw_z = Symbol(r'\delta\omega_z', real=True)
dw_z2 = Symbol(r'\delta\omega_{z2}', real=True)
dw_z4 = Symbol(r'\delta\omega_{z4}', real=True)
U = Symbol(r'U', real=True)
e = Symbol(r'\epsilon', real=True)
G_L = Symbol(r'\gamma_L', real=True)
dG_L4 = Symbol(r'\delta\gamma_{L4}', real=True)
dG_L2 = Symbol(r'\delta\gamma_{L2}', real=True)
G_R = Symbol(r'\gamma_R', real=True)
dG_R4 = Symbol(r'\delta\gamma_{R4}', real=True)
dG_R2 = Symbol(r'\delta\gamma_{R2}', real=True)
hbar = Symbol(r'\hbar', real=True)
w2 = Symbol(r'\omega_2', real=True)
w4 = Symbol(r'\omega_4', real=True)
ep2 = Symbol(r'\epsilon_{P2}', real=True)
ep4 = Symbol(r'\epsilon_{P4}', real=True)

H0 = np.array([
    [-dw_z  , 0      , G_R  , G_L      , -co(t) , -co(t)],
    [0      , dw_z   , G_L  , G_R      , t      , t],
    [G_R    , G_L    , w_z     , 0        , -Om    , -Om],
    [G_L, G_R, 0        , -w_z    , -co(Om), -co(Om)],
    [-t   , co(t), -co(Om), -Om    , U-e      , 0],
    [-t   , co(t), -co(Om), -Om    , 0        , U+e]])

H0 = np.array([
    [-dw_z  , 0      , 0  , 0     , -co(t) , -co(t)],
    [0      , dw_z   , 0  ,0      , t      , t],
    [0    , 0    , w_z     , 0        , -Om    , -Om],
    [0, 0, 0        , -w_z    , -co(Om), -co(Om)],
    [-t   , co(t), -co(Om), -Om    , U-e      , 0],
    [-t   , co(t), -co(Om), -Om    , 0        , U+e]])

V4 = np.array([
    [-dw_z4, 0, dG_R4, dG_L4, 0, 0],
    [0, dw_z4, dG_L4, dG_R4, 0, 0],
    [dG_R4, dG_L4, w_z4, 0, 0, 0],
    [dG_L4, dG_R4, 0, -w_z4, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])
V4 = np.array([
    [-dw_z4, 0, 0, 0, 0, 0],
    [0, dw_z4, 0, 0, 0, 0],
    [0, 0, w_z4, 0, 0, 0],
    [0, 0, 0, -w_z4, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])
V4 = np.array([
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -ep4, 0],
    [0, 0, 0, 0, 0, ep4]])
V2 = np.array([
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0 , 0, 0, 0, 0],
    [0, 0, 0, 0, -ep2, 0],
    [0, 0, 0, 0, 0,  ep2]])

# Crear instancia de FloquetSW (sin estados privilegiados en este caso)
sw = FloquetSW(H0, V2, V4, w2, w4, privileged=[])

alpha = 3
beta = 0

m = [alpha, 0, 0]
mp = m
lsn2n4 = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1]
]
ls = generate_ls_indices(lsn2n4, m, mp)
D1 = sw.SWF_2(m, mp, ls)


m = [beta, 1, -1]
mp = m
lsn2n4 = [
    [1, -1],
    [2, -1],
    [1, 0],
    [0, -1],
    [1, -2]
]
ls = generate_ls_indices(lsn2n4, m, mp)
D2 = sw.SWF_2(m, mp, ls)

    

m = [alpha, 0, 0]
mp = m
lsn2n4 = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1]
]
ls = generate_ls_indices(lsn2n4, m, mp)
D10 = sw.SWF_2(m, mp, ls)
m = [beta, 0, -1]
mp = m
lsn2n4 = [
    [0, -1],
    [1, -1],
    [0, 0],
    [-1, -1],
    [0, -2]
]
ls = generate_ls_indices(lsn2n4, m, mp)
D20 = sw.SWF_2(m, mp, ls)
init_printing()
with open("tfg_18.txt", "w") as f:
    f.write(latex(simplify(D1)))
    f.write("\n")
    f.write(latex(simplify(D2)))
    f.write("\n")
    f.write(latex(simplify(D10)))
    f.write("\n")
    f.write(latex(simplify(D20)))