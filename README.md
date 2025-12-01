# Teoría de Floquet en Puntos Cuánticos

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
Driven-Quantum-Dots/
├── README.md                 # Este archivo
├── src/                      # Código fuente principal
│   ├── __init__.py
│   ├── core/                 # Módulos fundamentales
│   │   ├── __init__.py
│   │   ├── schrieffer_wolff.py  # Transformación de Schrieffer-Wolff
│   │   ├── floquet.py        # Teoría de Floquet
│   │   ├── hamiltonian.py    # Segunda cuantización
│   │   └── quantum_dots.py   # Utilidades para QDs
│   ├── simulations/          # Scripts de simulación
│   │   ├── __init__.py
│   │   └── examples/         # Ejemplos de uso
│   │       ├── __init__.py
│   │       └── floquet_bichrom_example.py  # Ejemplo de Floquet bicromático
│   └── visualizations/       # Scripts de visualización
│       ├── __init__.py
│       └── rabi_heatmaps.py  # Mapas de calor de frecuencias de Rabi
└── scripts/                  # Scripts de simulación del TFG
    ├── comparison/           # Comparación de métodos (QuTiP vs Floquet)
    │   ├── tfg_2.py          # Evolución temporal QuTiP vs Floquet bicromático
    │   └── tfg_9.py          # Ground truth con QuTiP
    ├── schrieffer_wolff/     # Transformaciones de Schrieffer-Wolff
    │   ├── tfg_5.py          # SW simbólica básica
    │   ├── tfg_6.py          # SW con acoplamientos dependientes del tiempo
    │   ├── tfg_7.py          # Análisis perturbativo de 3er orden con g-TMR
    │   ├── tfg_8.py          # Cálculo numérico de frecuencias de Rabi
    │   └── tfg_18.py         # Correcciones diagonales de SW
    ├── interaction_picture/  # Transformaciones unitarias
    │   ├── tfg_11.py         # Transformación unitaria (versión inicial)
    │   └── tfg_12.py         # Picture de interacción (versión corregida)
    ├── frequency_sweeps/     # Barridos de frecuencias
    │   ├── tfg_13.py         # Barrido con excitación P2
    │   └── tfg_14.py         # Barrido con excitación P4
    ├── block_diagonalization/# Diagonalización con pymablock
    │   ├── tfg_15.py         # Excitación monocromática
    │   ├── tfg_16.py         # Excitación bicromática
    │   └── tfg_19.py         # Bicromática con Bloch-Siegert mejorado
    └── optimization/         # Optimización de parámetros
        └── tfg_17.py         # Optimización g-TMR
```

## Módulos Principales

### `src/core/schrieffer_wolff.py`

Implementa la transformación de Schrieffer-Wolff, una técnica perturbativa fundamental para:
- Proyectar Hamiltonianos a subespacios de baja energía
- Calcular correcciones perturbativas de segundo y tercer orden
- Derivar expresiones analíticas para frecuencias de Rabi y desplazamientos de Bloch-Siegert

Funciones principales:

- **`SW_transform_2(H, A_idcs, B_idcs)`**: Transformación de segundo orden estática
- **`FloquetSW`**: Clase para transformaciones SW en el espacio de Floquet
  - `SWF_2(m, mp, ls)`: Corrección de segundo orden en Floquet
  - `SWF_3(m, mp, lsA, lsB, mppsA, mppsB)`: Corrección de tercer orden en Floquet
- **`generate_ls_indices(base_n2n4_list, m, mp)`**: Genera índices de estados intermedios

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

## Scripts de Simulación

Los scripts de simulación están organizados en el directorio `scripts/` según su funcionalidad:

### Comparación de Métodos (`scripts/comparison/`)

| Script | Descripción |
|--------|-------------|
| `tfg_2.py` | Comparación de evolución temporal: QuTiP vs Floquet bicromático |
| `tfg_9.py` | Simulación de ground truth para validación con QuTiP (fsesolve) |

### Transformaciones de Schrieffer-Wolff (`scripts/schrieffer_wolff/`)

| Script | Descripción |
|--------|-------------|
| `tfg_5.py` | Implementación simbólica básica de la transformación de SW |
| `tfg_6.py` | SW con acoplamientos dependientes del tiempo (V2, V4) |
| `tfg_7.py` | Análisis perturbativo de tercer orden con efectos g-TMR |
| `tfg_8.py` | Cálculo numérico de frecuencias de Rabi con parámetros experimentales |
| `tfg_18.py` | Correcciones diagonales de segundo orden (Bloch-Siegert) |

### Picture de Interacción (`scripts/interaction_picture/`)

| Script | Descripción |
|--------|-------------|
| `tfg_11.py` | Transformación unitaria al picture de interacción (versión inicial) |
| `tfg_12.py` | Picture de interacción corregido (reproduce resultados de tfg_9) |

### Barridos de Frecuencias (`scripts/frequency_sweeps/`)

| Script | Descripción |
|--------|-------------|
| `tfg_13.py` | Barrido de frecuencias con modulación P2 (plunger 2) |
| `tfg_14.py` | Barrido de frecuencias con modulación P4 (plunger 4) |

### Diagonalización de Bloques con pymablock (`scripts/block_diagonalization/`)

| Script | Descripción |
|--------|-------------|
| `tfg_15.py` | Frecuencias de Rabi para excitación monocromática |
| `tfg_16.py` | Frecuencias de Rabi para excitación bicromática |
| `tfg_19.py` | Versión final con cálculo mejorado de Bloch-Siegert |

### Optimización (`scripts/optimization/`)

| Script | Descripción |
|--------|-------------|
| `tfg_17.py` | Optimización de parámetros g-TMR para ajustar datos experimentales |

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
git clone https://github.com/papipopapu/Driven-Quantum-Dots.git
cd Driven-Quantum-Dots

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
H0 = np.array([...])  # Ver scripts/comparison/tfg_9.py para la forma completa

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


## Licencia

Este proyecto es parte de un Trabajo de Fin de Grado académico.

## Contacto

Para preguntas sobre el código o la física, contactar al autor del TFG.
