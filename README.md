"# Teoría de Floquet en Puntos Cuánticos

## Descripción

Este repositorio contiene el código del Trabajo de Fin de Grado (TFG) que estudia la dinámica de **puntos cuánticos** bajo excitación periódica utilizando **teoría de Floquet**. El proyecto analiza la evolución temporal de sistemas de pocos electrones en puntos cuánticos semiconductores cuando son excitados con campos electromagnéticos oscilantes a diferentes frecuencias.

### Objetivos del proyecto

1. Implementar el formalismo de Floquet para Hamiltonianos dependientes del tiempo
2. Estudiar resonancias de Rabi en puntos cuánticos dobles
3. Analizar excitación monocromática y bicromática
4. Calcular correcciones perturbativas mediante transformaciones de Schrieffer-Wolff

## Fundamentos Teóricos

### Teoría de Floquet

Para un Hamiltoniano periódico en el tiempo $H(t) = H(t+T)$, el teorema de Floquet garantiza que las soluciones de la ecuación de Schrödinger tienen la forma:

$$|\psi(t)\rangle = e^{-i\varepsilon_\alpha t/\hbar} |\phi_\alpha(t)\rangle$$

donde:
- $\varepsilon_\alpha$ son las **quasi-energías de Floquet**
- $|\phi_\alpha(t)\rangle$ son los **modos de Floquet** (periódicos con periodo $T$)

El **Hamiltoniano de Floquet** $H_F$ actúa en un espacio extendido que incluye los modos de Fourier, permitiendo tratar el problema dependiente del tiempo como un problema de autovalores estático.

### Sistema físico

El sistema estudiado consiste en un **punto cuántico doble** con:
- Estados de espín (↑, ↓) en cada punto
- Acoplamiento túnel entre puntos
- Interacción de Coulomb
- Excitación externa con campos AC

El Hamiltoniano base tiene la estructura:

$$H(t) = H_0 + V_1 \cos(\omega_1 t) + V_2 \cos(\omega_2 t)$$

donde las matrices de acoplamiento $V$ representan la modulación de los niveles de energía.

## Estructura del Proyecto

```
TFG/
├── README.md                 # Este archivo
├── src/                      # Código fuente principal
│   ├── __init__.py
│   ├── core/                 # Módulos fundamentales
│   │   ├── __init__.py
│   │   ├── floquet.py        # Teoría de Floquet
│   │   ├── hamiltonian.py    # Segunda cuantización
│   │   └── quantum_dots.py   # Utilidades para QDs
│   ├── simulations/          # Scripts de simulación
│   │   └── __init__.py
│   └── visualizations/       # Scripts de visualización
│       └── __init__.py
├── examples/                 # Ejemplos de uso
├── data/                     # Datos y resultados
├── floquet_theory.py         # [Legacy] Módulo original de Floquet
├── floquet_bi.py             # [Legacy] Floquet bicromático original
├── hamiltonian.py            # [Legacy] Hamiltoniano original
├── qutipDots.py              # [Legacy] Utilidades QuTiP original
├── plots.py                  # [Legacy] Visualización de resultados
└── tfg_*.py                  # [Legacy] Scripts de simulación numerados
```

## Módulos Principales

### `src/core/floquet.py`

Implementa el formalismo de Floquet:

- **`floquet_mono(H0, Nw, V, w, phi)`**: Construye el Hamiltoniano de Floquet para excitación monocromática
- **`floquet_bichrom(H0, Nw, V, Vbar, w, wbar, phi, phibar)`**: Hamiltoniano de Floquet bicromático
- **`get_U(vecs, engs, w1, w2, t, t0, N, Nw)`**: Operador de evolución temporal
- **`get_P(alpha, beta, vecs, N, Nw)`**: Probabilidad de transición promediada

### `src/core/hamiltonian.py`

Construcción de Hamiltonianos con segunda cuantización:

- **`c_internal`**: Clase para operadores fermiónicos
- **`calc_Hamiltonian(H, basis)`**: Construye matriz del Hamiltoniano
- **`Space`**: Define espacios de Hilbert y productos tensoriales

### `src/core/quantum_dots.py`

Utilidades específicas para puntos cuánticos:

- **`f_create/f_destroy`**: Operadores fermiónicos con Jordan-Wigner
- **`get_Liouville(Gamma, H)`**: Superoperador para dinámica disipativa
- **`red_H`**: Reducción a subespacios

## Scripts de Simulación (Legacy)

Los archivos `tfg_*.py` contienen diferentes aspectos del estudio:

| Script | Descripción |
|--------|-------------|
| `tfg_2.py` | Evolución temporal con QuTiP vs Floquet |
| `tfg_5.py`, `tfg_6.py` | Transformaciones de Schrieffer-Wolff simbólicas |
| `tfg_7.py`, `tfg_8.py` | Análisis perturbativo de tercer orden |
| `tfg_9.py` | Simulación de ground truth con QuTiP |
| `tfg_11.py`, `tfg_12.py` | Transformación unitaria y picture de interacción |
| `tfg_13.py` | Barrido de frecuencias (excitación con P2) |
| `tfg_14.py` | Barrido de frecuencias (excitación con P4) |
| `tfg_15.py` | Diagonalización de bloque con pymablock (monocromático) |
| `tfg_16.py`, `tfg_19.py` | Diagonalización de bloque (bicromático) |
| `tfg_17.py` | Optimización de parámetros |
| `tfg_18.py` | Cálculos simbólicos adicionales |

## Dependencias

```python
numpy>=1.20
scipy>=1.7
matplotlib>=3.4
qutip>=4.6
numba>=0.54
sympy>=1.9
pymablock>=0.1
tqdm>=4.62
scienceplots>=1.0  # Opcional, para estilo de gráficos
```

### Instalación

```bash
# Clonar el repositorio
git clone https://github.com/papipopapu/TFG.git
cd TFG

# Instalar dependencias
pip install numpy scipy matplotlib qutip numba sympy tqdm
pip install pymablock  # Para diagonalización de bloques

# Opcional: estilo de gráficos científicos
pip install SciencePlots
```

## Uso

### Ejemplo básico: Evolución con Floquet monocromático

```python
import numpy as np
from src.core.floquet import floquet_mono, get_U_mono

# Parámetros del sistema (en unidades de 2π GHz)
w_z = 2.00 * 2*np.pi    # Splitting Zeeman
dw_z = 0.439 * 2*np.pi  # Diferencia de g-factor
t = 4.38 * 2*np.pi      # Túnel
Om = 3.46 * 2*np.pi     # Otro acoplamiento
U = 619 * 2*np.pi       # Interacción de Coulomb
e = 442 * 2*np.pi       # Detuning de carga

# Hamiltoniano estático (6x6 para punto cuántico doble)
H0 = np.array([...])  # Ver tfg_9.py para la forma completa

# Matriz de acoplamiento
V = np.array([...])  # Excitación en los estados de carga

# Construir Hamiltoniano de Floquet
Nw = 10  # Modos de Fourier
w = 1.5 * 2*np.pi  # Frecuencia de excitación
HF = floquet_mono(H0, Nw, V, w)

# Diagonalizar
vals, vecs = np.linalg.eigh(HF)

# Evolucionar en el tiempo
t_final = 100  # ns
U_evol = get_U_mono(vecs, vals, w, t_final, 0, 6, Nw)

# Probabilidad de transición
psi0 = np.array([0, 0, 0, 1, 0, 0])  # Estado inicial
psi_t = U_evol @ psi0
prob = np.abs(psi_t)**2
```

### Ejemplo con QuTiP (comparación)

```python
import qutip as qt

H0_qt = qt.Qobj(H0)
V_qt = qt.Qobj(V)

def V_coeff(t, args):
    return np.cos(args['w'] * t)

H = [H0_qt, [V_qt, V_coeff]]
tlist = np.linspace(0, 100, 1000)

psi0_qt = qt.basis(6, 3)  # Estado |↓↓⟩
result = qt.mesolve(H, psi0_qt, tlist, [], [psi0_qt * psi0_qt.dag()], 
                    args={'w': w})
```

## Resultados Principales

El proyecto demuestra:

1. **Resonancias de Rabi**: Para frecuencias de excitación que cumplen condiciones de resonancia ($\omega \approx \varepsilon_f - \varepsilon_i$), se observan oscilaciones coherentes entre estados.

2. **Excitación bicromática**: Condiciones de resonancia extendidas:
   - $n_1\omega_1 + n_2\omega_2 = \varepsilon_f - \varepsilon_i$
   - $n_1\omega_1 - n_2\omega_2 = \varepsilon_f - \varepsilon_i$

3. **Efecto g-TMR**: La modulación del tensor g produce correcciones adicionales a las frecuencias de Rabi.

4. **Corrección de Bloch-Siegert**: Desplazamientos de frecuencia por efectos contra-rotantes.

## Referencias

- Shirley, J. H. (1965). *Solution of the Schrödinger Equation with a Hamiltonian Periodic in Time*. Physical Review.
- Sambe, H. (1973). *Steady States and Quasienergies of a Quantum-Mechanical System in an Oscillating Field*. Physical Review A.
- Artículos experimentales sobre puntos cuánticos dobles en GaAs/AlGaAs y Si/SiGe.

## Licencia

Este proyecto es parte de un Trabajo de Fin de Grado académico.

## Contacto

Para preguntas sobre el código o la física, contactar al autor del TFG." 
