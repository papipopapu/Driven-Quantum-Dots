"""
Módulo de utilidades para puntos cuánticos (quantum dots).

Este módulo proporciona funciones auxiliares para trabajar con sistemas
de puntos cuánticos, incluyendo:
- Operadores fermiónicos con la convención de Jordan-Wigner
- Construcción de estados de ocupación
- Operadores de Liouville para dinámica disipativa
- Reducción de Hamiltonianos a subespacios

Teoría:
-------
Los puntos cuánticos son sistemas de pocos electrones confinados en
potenciales semiconductores. El Hamiltoniano típico incluye:
- Energía de confinamiento
- Interacción de Coulomb
- Acoplamiento espín-órbita
- Interacción Zeeman con campo magnético externo

Para simulaciones con QuTiP, es necesario mapear los operadores
fermiónicos a operadores de espín usando la transformación de Jordan-Wigner.
"""

import numpy as np
from numpy.typing import ArrayLike
from typing import List
from qutip import tensor, basis, sigmaz, destroy, identity


def f_destroy(N: int, i: int):
    """
    Operador fermiónico de destrucción usando transformación de Jordan-Wigner.
    
    La transformación de Jordan-Wigner mapea operadores fermiónicos a
    operadores de espín, preservando las relaciones de anticonmutación:
    
        c_i = (∏_{j<i} σ_z^{(j)}) ⊗ σ_-^{(i)} ⊗ I ⊗ ... ⊗ I
    
    donde σ_- = (σ_x - iσ_y)/2 es el operador de bajada de espín.
    
    Parameters
    ----------
    N : int
        Número total de modos fermiónicos (orbitales × espines).
    i : int
        Índice del modo a destruir (0 a N-1).
    
    Returns
    -------
    Qobj
        Operador de destrucción fermiónico en representación de QuTiP.
    
    Notes
    -----
    El producto de σ_z antes del sitio i implementa la "cadena de Jordan-Wigner"
    que asegura las relaciones de anticonmutación fermiónicos.
    
    Example
    -------
    >>> c0 = f_destroy(2, 0)  # Destruye fermión en modo 0
    >>> c1 = f_destroy(2, 1)  # Destruye fermión en modo 1
    >>> # Verificar anticonmutación: {c0, c1†} = 0
    """
    return tensor([sigmaz()] * i + [destroy(2)] + [identity(2)] * (N - i - 1))


def f_create(N: int, i: int):
    """
    Operador fermiónico de creación usando transformación de Jordan-Wigner.
    
    Es simplemente el conjugado hermítico del operador de destrucción.
    
    Parameters
    ----------
    N : int
        Número total de modos fermiónicos.
    i : int
        Índice del modo a crear.
    
    Returns
    -------
    Qobj
        Operador de creación fermiónico (c†_i).
    """
    return f_destroy(N, i).dag()


def eqdot_state(occupations: List[bool]):
    """
    Construye un estado de Fock para un punto cuántico.
    
    Dado un patrón de ocupaciones, construye el estado producto
    correspondiente donde cada modo está ocupado (1) o vacío (0).
    
    Parameters
    ----------
    occupations : list of bool
        Lista indicando la ocupación de cada modo.
        True/1 = ocupado, False/0 = vacío.
    
    Returns
    -------
    Qobj
        Estado de Fock |n_1, n_2, ..., n_N⟩.
    
    Example
    -------
    >>> # Estado con electrón en modo 0, vacío en modo 1
    >>> psi = eqdot_state([True, False])
    >>> # Estado completamente vacío (2 modos)
    >>> vacuum = eqdot_state([False, False])
    """
    N = len(occupations)
    return tensor([basis(2, occupations[i]) for i in range(N)])


def get_Lambda(Gamma: ArrayLike) -> ArrayLike:
    """
    Calcula la matriz de decoherencia Lambda a partir de tasas de transición.
    
    La matriz Lambda aparece en la ecuación maestra de Lindblad y describe
    la decoherencia entre estados debida a transiciones incoherentes:
    
        Λ_{mn} = (1/2) Σ_{k≠n} (Γ_{km} + Γ_{kn})
    
    donde Γ_{ij} es la tasa de transición del estado j al estado i.
    
    Parameters
    ----------
    Gamma : ndarray
        Matriz de tasas de transición Γ_{ij} (dimensión DxD).
    
    Returns
    -------
    Lambda : ndarray
        Matriz de decoherencia (DxD).
    
    Notes
    -----
    Esta matriz aparece en los elementos off-diagonal de la ecuación
    maestra de Lindblad:
        ∂ρ_{mn}/∂t = ... - Λ_{mn}ρ_{mn}  (para m ≠ n)
    """
    D = Gamma.shape[0]
    Lambda = np.zeros((D, D))
    
    for m in range(D):
        for n in range(D):
            Lambda_mn = 0
            for k in range(D):
                if k != n:
                    Lambda_mn += 1/2 * (Gamma[k, m] + Gamma[k, n])
            Lambda[m, n] = Lambda_mn
    
    return Lambda


def get_Liouville(Gamma: ArrayLike, H: ArrayLike) -> ArrayLike:
    """
    Construye el superoperador de Liouville para dinámica disipativa.
    
    El superoperador de Liouville L gobierna la evolución de la matriz
    densidad en forma vectorizada:
        d|ρ⟩⟩/dt = L|ρ⟩⟩
    
    donde |ρ⟩⟩ es la vectorización de ρ (columnas apiladas).
    
    El superoperador incluye:
    - Evolución coherente: -i[H, ρ]
    - Transiciones incoherentes (diagonal): Σ_k Γ_{nk}ρ_{kk}
    - Decoherencia (off-diagonal): -Λ_{mn}ρ_{mn}
    
    Parameters
    ----------
    Gamma : ndarray
        Matriz de tasas de transición (DxD).
    H : ndarray
        Hamiltoniano del sistema (DxD).
    
    Returns
    -------
    Liouville : ndarray
        Superoperador de Liouville (D²×D²).
    
    Notes
    -----
    La convención de vectorización es: |ρ⟩⟩_{m+nD} = ρ_{mn}
    
    Esto corresponde a la ecuación maestra de Lindblad en la
    aproximación secular.
    """
    D = Gamma.shape[0]
    Liouville = np.zeros((D*D, D*D), dtype=np.complex128)
    I = np.eye(D)
    
    # Parte coherente: -i(H⊗I - I⊗H^T)
    Liouville += -1j * (np.kron(I, H) - np.kron(H.T, I))
    
    # Parte disipativa
    Lambda = get_Lambda(Gamma)
    
    for m in range(D):
        for n in range(D):
            mn_idx = m + n*D
            if m == n:
                # Poblaciones: dρ_{nn}/dt = Σ_k Γ_{nk}ρ_{kk} - Σ_k Γ_{kn}ρ_{nn}
                for k in range(D):
                    if k != n:
                        Liouville[mn_idx, k + k*D] += Gamma[n, k]  # Entrada
                        Liouville[mn_idx, n + n*D] -= Gamma[k, n]  # Salida
            else:
                # Coherencias: dρ_{mn}/dt = ... - Λ_{mn}ρ_{mn}
                Liouville[mn_idx, mn_idx] -= Lambda[m, n]
        
    return Liouville


def red_H_idx(H, allowed_idx):
    """
    Reduce un Hamiltoniano a un subespacio por índices.
    
    Proyecta el Hamiltoniano al subespacio definido por los índices dados,
    manteniendo solo los elementos de matriz correspondientes.
    
    Parameters
    ----------
    H : ndarray
        Hamiltoniano completo.
    allowed_idx : list of int
        Índices de los estados del subespacio.
    
    Returns
    -------
    H_red : ndarray
        Hamiltoniano reducido (len(allowed_idx) × len(allowed_idx)).
    
    Example
    -------
    >>> # Reducir a los primeros 3 estados
    >>> H_sub = red_H_idx(H_full, [0, 1, 2])
    """
    N = len(allowed_idx)
    H_red = np.zeros((N, N), dtype=np.complex128)
    
    for i in range(N):
        for j in range(N):
            H_red[i, j] = H[allowed_idx[i], allowed_idx[j]]
    
    return H_red


def red_H(H, states):
    """
    Reduce un Hamiltoniano a un subespacio por estados.
    
    Proyecta el Hamiltoniano usando un conjunto de estados (objetos Qobj).
    Calcula H_red[i,j] = ⟨states[i]|H|states[j]⟩.
    
    Parameters
    ----------
    H : Qobj
        Hamiltoniano como operador de QuTiP.
    states : list of Qobj
        Estados que definen el subespacio.
    
    Returns
    -------
    H_red : ndarray
        Matriz del Hamiltoniano en el subespacio.
    
    Notes
    -----
    Esta función es útil cuando los estados del subespacio no son
    simplemente estados base, sino superposiciones (ej. estados vestidos).
    """
    N = len(states)
    H_red = np.zeros((N, N), dtype=np.complex128)
    
    for i in range(N):
        for j in range(N):
            H_red[i, j] = (states[i].dag() * H * states[j])
    
    return H_red
