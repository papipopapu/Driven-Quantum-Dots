"""
Módulo de Teoría de Floquet para sistemas cuánticos periódicamente excitados.

Este módulo implementa la teoría de Floquet para sistemas cuánticos con
Hamiltonianos dependientes del tiempo de forma periódica. Incluye funciones
para construir el Hamiltoniano de Floquet tanto para excitación monocromática
como bicromática, así como el cálculo del operador de evolución temporal.

Teoría:
-------
Para un Hamiltoniano periódico H(t) = H(t+T), el teorema de Floquet establece
que las soluciones de la ecuación de Schrödinger se pueden escribir como:

    |ψ(t)⟩ = exp(-iε_α t/ℏ) |φ_α(t)⟩

donde ε_α son las quasi-energías de Floquet y |φ_α(t)⟩ son los modos de Floquet.

El Hamiltoniano de Floquet H_F se construye en un espacio extendido (espacio
de Floquet) que incluye los índices de los modos de Fourier del Hamiltoniano.

Referencias:
-----------
- Shirley, J. H. (1965). Solution of the Schrödinger Equation with a 
  Hamiltonian Periodic in Time. Physical Review.
- Sambe, H. (1973). Steady States and Quasienergies of a Quantum-Mechanical 
  System in an Oscillating Field. Physical Review A.
"""

import numpy as np
import numba as nb


def floquet_bichrom(H0, Nw, V, Vbar, w, wbar, phi=0, phibar=0):
    """
    Construye el Hamiltoniano de Floquet para excitación bicromática.
    
    Para un sistema con dos frecuencias de excitación ω y ω̄, el Hamiltoniano
    tiene la forma:
        H(t) = H₀ + V·cos(ωt + φ) + V̄·cos(ω̄t + φ̄)
    
    El Hamiltoniano de Floquet se construye en el espacio extendido con
    índices (α, n, n̄) donde α es el índice del sistema y n, n̄ son los
    índices de los modos de Fourier.
    
    Parameters
    ----------
    H0 : ndarray
        Hamiltoniano estático del sistema (matriz NxN).
    Nw : int
        Número de modos de Fourier a incluir (-Nw a +Nw).
    V : ndarray
        Matriz de acoplamiento para la primera frecuencia.
    Vbar : ndarray
        Matriz de acoplamiento para la segunda frecuencia.
    w : float
        Primera frecuencia de excitación ω.
    wbar : float
        Segunda frecuencia de excitación ω̄.
    phi : float, optional
        Fase de la primera excitación (default: 0).
    phibar : float, optional
        Fase de la segunda excitación (default: 0).
    
    Returns
    -------
    HF : ndarray
        Hamiltoniano de Floquet en el espacio extendido.
        Dimensión: (2*Nw+1)² * N x (2*Nw+1)² * N
    
    Notes
    -----
    La indexación del Hamiltoniano es [α, β, n₁, n₂] donde:
    - α, β son índices del sistema
    - n₁ corresponde a la frecuencia w
    - n₂ corresponde a la frecuencia wbar
    
    El Hamiltoniano de Floquet cumple la ecuación de eigenvalores:
        H_F |ε_α, n₁, n₂⟩ = ε_α |ε_α, n₁, n₂⟩
    
    donde las quasi-energías ε_α están definidas módulo ℏω y ℏω̄.
    """
    # Matrices de Kronecker delta
    delta_alpha_beta = np.eye(H0.shape[0], dtype=np.complex128)
    delta_n_k = np.eye(2*Nw+1, dtype=np.complex128)
    
    # Matrices de desplazamiento para acoplamientos n-k = ±1
    delta_nminusk_1 = np.diag(np.ones(2*Nw, dtype=np.complex128), k=-1)
    delta_nminusk_minus1 = np.diag(np.ones(2*Nw, dtype=np.complex128), k=1)
    
    # Contribuciones diagonales modificadas
    delta_alpha_beta_mod = H0
    delta_n_k_mod = np.diag(np.arange(-Nw, Nw+1, dtype=np.complex128) * w)
    delta_nbar_kbar_mod = np.diag(np.arange(-Nw, Nw+1, dtype=np.complex128) * wbar)

    # Construcción del Hamiltoniano de Floquet
    # Orden de productos de Kronecker: n_k x nbar_kbar x alpha_beta
    HF = np.zeros(((2*Nw+1)**2 * H0.shape[0], (2*Nw+1)**2 * H0.shape[0]), dtype=np.complex128)
    
    # Término del Hamiltoniano estático
    HF += np.kron(delta_n_k, np.kron(delta_n_k, delta_alpha_beta_mod))
    
    # Términos de quasi-energía
    HF += np.kron(delta_n_k, np.kron(delta_n_k_mod, delta_alpha_beta))
    HF += np.kron(delta_nbar_kbar_mod, np.kron(delta_n_k, delta_alpha_beta))
    
    # Términos de acoplamiento para frecuencia w (V)
    HF += np.kron(delta_n_k, np.kron(delta_nminusk_1, V/2 * np.exp(1j * phi)))
    HF += np.kron(delta_n_k, np.kron(delta_nminusk_minus1, V/2 * np.exp(-1j * phi)))
    
    # Términos de acoplamiento para frecuencia wbar (Vbar)
    HF += np.kron(delta_nminusk_1, np.kron(delta_n_k, Vbar/2 * np.exp(1j * phibar)))
    HF += np.kron(delta_nminusk_minus1, np.kron(delta_n_k, Vbar/2 * np.exp(-1j * phibar)))
    
    return HF


def floquet_mono(H0, Nw, V, w, phi=0):
    """
    Construye el Hamiltoniano de Floquet para excitación monocromática.
    
    Para un sistema con una frecuencia de excitación ω, el Hamiltoniano
    tiene la forma:
        H(t) = H₀ + V·cos(ωt + φ)
    
    Parameters
    ----------
    H0 : ndarray
        Hamiltoniano estático del sistema (matriz NxN).
    Nw : int
        Número de modos de Fourier a incluir (-Nw a +Nw).
    V : ndarray
        Matriz de acoplamiento para la excitación.
    w : float
        Frecuencia de excitación ω.
    phi : float, optional
        Fase de la excitación (default: 0).
    
    Returns
    -------
    HF : ndarray
        Hamiltoniano de Floquet en el espacio extendido.
        Dimensión: (2*Nw+1) * N x (2*Nw+1) * N
    
    Notes
    -----
    La indexación del Hamiltoniano es [α, β, n] donde:
    - α, β son índices del sistema
    - n corresponde al modo de Fourier
    
    Este es un caso especial del Hamiltoniano bicromático donde
    la segunda frecuencia es cero.
    """
    # Matrices de Kronecker delta
    delta_alpha_beta = np.eye(H0.shape[0], dtype=np.complex128)
    delta_n_k = np.eye(2*Nw+1, dtype=np.complex128)
    
    # Matrices de desplazamiento para acoplamientos n-k = ±1
    delta_nminusk_1 = np.diag(np.ones(2*Nw, dtype=np.complex128), k=-1)
    delta_nminusk_minus1 = np.diag(np.ones(2*Nw, dtype=np.complex128), k=1)
    
    # Contribuciones diagonales modificadas
    delta_alpha_beta_mod = H0
    delta_n_k_mod = np.diag(np.arange(-Nw, Nw+1) * w)
    
    # Construcción del Hamiltoniano de Floquet
    # Orden de productos de Kronecker: n_k x alpha_beta
    HF = np.zeros(((2*Nw+1) * H0.shape[0], (2*Nw+1) * H0.shape[0]), dtype=np.complex128)
    
    # Término del Hamiltoniano estático
    HF += np.kron(delta_n_k, delta_alpha_beta_mod)
    
    # Término de quasi-energía
    HF += np.kron(delta_n_k_mod, delta_alpha_beta)
    
    # Términos de acoplamiento
    HF += np.kron(delta_nminusk_1, V/2 * np.exp(1j * phi))
    HF += np.kron(delta_nminusk_minus1, V/2 * np.exp(-1j * phi))
    
    return HF


@nb.njit
def idx_mono(alpha, n, N, Nw):
    """
    Calcula el índice lineal para el espacio de Floquet monocromático.
    
    Parameters
    ----------
    alpha : int
        Índice del estado del sistema (0 a N-1).
    n : int
        Índice del modo de Fourier (-Nw a +Nw).
    N : int
        Dimensión del sistema.
    Nw : int
        Número máximo de modos de Fourier.
    
    Returns
    -------
    int
        Índice lineal en el espacio de Floquet.
    """
    return alpha + n * N + Nw * N


@nb.njit
def idx(alpha, n1, n2, N, Nw):
    """
    Calcula el índice lineal para el espacio de Floquet bicromático.
    
    La indexación sigue el orden: idx = α + n₁·N + n₂·N·(2Nw+1) + offset
    
    Parameters
    ----------
    alpha : int
        Índice del estado del sistema (0 a N-1).
    n1 : int
        Índice del primer modo de Fourier (-Nw a +Nw).
    n2 : int
        Índice del segundo modo de Fourier (-Nw a +Nw).
    N : int
        Dimensión del sistema.
    Nw : int
        Número máximo de modos de Fourier.
    
    Returns
    -------
    int
        Índice lineal en el espacio de Floquet bicromático.
    """
    return alpha + n1 * N + n2 * N * (2*Nw+1) + Nw * Nw * N


@nb.njit
def get_U_mono(vecs, engs, w, t, t0, N, Nw):
    """
    Calcula el operador de evolución temporal para excitación monocromática.
    
    El operador de evolución U(t, t₀) conecta estados en tiempos diferentes:
        |ψ(t)⟩ = U(t, t₀)|ψ(t₀)⟩
    
    Usando la descomposición de Floquet:
        U_{αβ}(t, t₀) = Σ_{n,γ} exp(-iε_γ(t-t₀)) u_{α,n}^{(γ)} [u_{β,0}^{(γ)}]* exp(inωt)
    
    Parameters
    ----------
    vecs : ndarray
        Autovectores del Hamiltoniano de Floquet.
    engs : ndarray
        Quasi-energías de Floquet.
    w : float
        Frecuencia de excitación.
    t : float
        Tiempo final.
    t0 : float
        Tiempo inicial.
    N : int
        Dimensión del sistema.
    Nw : int
        Número de modos de Fourier.
    
    Returns
    -------
    U : ndarray
        Matriz de evolución temporal (NxN).
    """
    U = np.zeros((N, N), dtype=np.complex128)
    Ntot = (2*Nw+1) * N
    
    for alpha in np.arange(N):
        for beta in np.arange(N):
            for n in np.arange(-Nw, Nw+1):
                for gamma_l in np.arange(Ntot):
                    # Factor de fase de quasi-energía
                    a = np.exp(-1j * engs[gamma_l] * (t-t0))
                    # Coeficiente del modo de Floquet
                    b = vecs[idx_mono(beta, n, N, Nw), gamma_l]
                    c = np.conj(vecs[idx_mono(alpha, 0, N, Nw), gamma_l])
                    # Factor de fase temporal
                    d = np.exp(1j * n * w * t)
                    U[beta, alpha] += a * b * c * d
                    
    return U


@nb.njit
def get_U(vecs, engs, w1, w2, t, t0, N, Nw):
    """
    Calcula el operador de evolución temporal para excitación bicromática.
    
    Generalización del caso monocromático a dos frecuencias:
        U_{αβ}(t, t₀) = Σ_{n₁,n₂,γ} exp(-iε_γ(t-t₀)) u_{α,n₁,n₂}^{(γ)} 
                        [u_{β,0,0}^{(γ)}]* exp(i(n₁ω₁ + n₂ω₂)t)
    
    Parameters
    ----------
    vecs : ndarray
        Autovectores del Hamiltoniano de Floquet bicromático.
    engs : ndarray
        Quasi-energías de Floquet.
    w1 : float
        Primera frecuencia de excitación.
    w2 : float
        Segunda frecuencia de excitación.
    t : float
        Tiempo final.
    t0 : float
        Tiempo inicial.
    N : int
        Dimensión del sistema.
    Nw : int
        Número de modos de Fourier.
    
    Returns
    -------
    U : ndarray
        Matriz de evolución temporal (NxN).
    """
    U = np.zeros((N, N), dtype=np.complex128)
    Ntot = (2*Nw+1) * (2*Nw+1) * N
    
    for alpha in np.arange(N):
        for beta in np.arange(N):
            for n1 in np.arange(-Nw, Nw+1):
                for n2 in np.arange(-Nw, Nw+1):
                    for gamma_k1k2 in np.arange(Ntot):
                        # Factor de fase de quasi-energía
                        a = np.exp(-1j * engs[gamma_k1k2] * (t-t0))
                        # Coeficiente del modo de Floquet
                        b = vecs[idx(alpha, n1, n2, N, Nw), gamma_k1k2]
                        c = np.conj(vecs[idx(beta, 0, 0, N, Nw), gamma_k1k2])
                        # Factor de fase temporal
                        d = np.exp(1j * (n1 * w1 + n2 * w2) * t)
                        U[beta, alpha] += a * b * c * d
                        
    return U


@nb.njit
def get_P(alpha, beta, vecs, N, Nw):
    """
    Calcula la probabilidad de transición promediada en tiempo (bicromático).
    
    Para transiciones periódicas, la probabilidad promediada sobre un ciclo
    completo se puede calcular directamente de los modos de Floquet:
        P_{α→β} = Σ_{n₁,n₂,γ} |u_{β,n₁,n₂}^{(γ)}|² |u_{α,0,0}^{(γ)}|²
    
    Parameters
    ----------
    alpha : int
        Estado inicial.
    beta : int
        Estado final.
    vecs : ndarray
        Autovectores del Hamiltoniano de Floquet.
    N : int
        Dimensión del sistema.
    Nw : int
        Número de modos de Fourier.
    
    Returns
    -------
    P : float
        Probabilidad de transición promediada.
    """
    P = 0
    Ntot = (2*Nw+1) * (2*Nw+1) * N
    
    for n1 in np.arange(-Nw, Nw+1):
        for n2 in np.arange(-Nw, Nw+1):
            for gamma_k1k2 in np.arange(Ntot):
                P += np.abs(vecs[idx(beta, n1, n2, N, Nw), gamma_k1k2] * 
                           np.conj(vecs[idx(alpha, 0, 0, N, Nw), gamma_k1k2]))**2
    return P


@nb.njit
def get_P_mono(alpha, beta, vecs, N, Nw):
    """
    Calcula la probabilidad de transición promediada en tiempo (monocromático).
    
    Versión monocromática del cálculo de probabilidad:
        P_{α→β} = Σ_{n,γ} |u_{β,n}^{(γ)}|² |u_{α,0}^{(γ)}|²
    
    Parameters
    ----------
    alpha : int
        Estado inicial.
    beta : int
        Estado final.
    vecs : ndarray
        Autovectores del Hamiltoniano de Floquet.
    N : int
        Dimensión del sistema.
    Nw : int
        Número de modos de Fourier.
    
    Returns
    -------
    P : float
        Probabilidad de transición promediada.
    """
    P = 0
    Ntot = (2*Nw+1) * N
    
    for n in np.arange(-Nw, Nw+1):
        for gamma_k1k2 in np.arange(Ntot):
            P += np.abs(vecs[idx_mono(beta, n, N, Nw), gamma_k1k2] * 
                       np.conj(vecs[idx(alpha, 0, 0, N, Nw), gamma_k1k2]))**2
    return P
