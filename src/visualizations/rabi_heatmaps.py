"""
Visualización de resultados: mapas de calor de frecuencias de Rabi.

Este script carga datos pre-calculados de frecuencias de Rabi y genera
mapas de calor 2D mostrando la dependencia con los parámetros g-TMR
(w_z4 y dw_z4).

Los datos se generan en tfg_15.py (monocromático) y tfg_16.py (bicromático).
"""

import numpy as np
import matplotlib.pyplot as plt

rabisQ1P2P4 = np.load('rabis2.npy')
rabisQ1_P2P4 = np.load('rabis_2.npy')
rabisQ1P4 = np.load('rabis.npy')
rabisQ1_P4 = np.load('rabis_.npy')
db = abs(rabisQ1P2P4 - rabisQ1_P2P4)
dm = abs(rabisQ1P4 - rabisQ1_P4)
w_z4s = np.load('w_z4s2.npy')
dw_z4s = np.load('dw_z4s2.npy')

# Calculate extent for imshow plots
extent = [dw_z4s[0]/(2*np.pi), dw_z4s[-1]/(2*np.pi), 
          w_z4s[0]/(2*np.pi), w_z4s[-1]/(2*np.pi)]

fig, ax = plt.subplots()
im = ax.imshow(db, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$db$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots()
im = ax.imshow(dm, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$dm$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots()
im = ax.imshow(rabisQ1P4, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$Q1^{P4}$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots()
im = ax.imshow(rabisQ1_P4, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$Q1\_^{P4}$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots()
im = ax.imshow(rabisQ1P2P4, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$Q1^{P2,P4}$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots()
im = ax.imshow(rabisQ1_P2P4, extent=extent, aspect='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.ax.tick_params(labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15, width=1.2, length=6)
ax.tick_params(axis='both', which='minor', labelsize=15, width=1.2, length=2)
ax.set_xlabel(r'$\delta \omega_{z4}$ (GHz)', fontsize=15)
ax.set_ylabel(r'$\omega_{z4}$ (GHz)', fontsize=15)
ax.set_title(r'$Q1^{-P2,P4}$ Rabi Frequency (MHz)', fontsize=15)
plt.tight_layout()

plt.show()
plt.show()