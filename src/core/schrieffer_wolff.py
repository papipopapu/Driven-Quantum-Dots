"""
Transformación de Schrieffer-Wolff para puntos cuánticos.

Este módulo implementa la transformación de Schrieffer-Wolff, una técnica
perturbativa que proyecta un Hamiltoniano a un subespacio de baja energía
eliminando los acoplamientos con estados de alta energía.

La transformación es fundamental para:
- Calcular Hamiltonianos efectivos de sistemas de pocos niveles
- Obtener correcciones perturbativas a las energías y acoplamientos
- Derivar expresiones analíticas para frecuencias de Rabi
- Calcular el desplazamiento de Bloch-Siegert

Teoría:
------
Para un Hamiltoniano H = H_A + H_B + V, donde H_A actúa en el subespacio A
(baja energía) y H_B en el subespacio B (alta energía), el Hamiltoniano
efectivo hasta segundo orden es:

H_eff = H_A + (1/2) Σ_l [V_{al} V_{la} (1/(E_a - E_l) + 1/(E_a' - E_l))]

Este módulo extiende esta idea al espacio de Floquet para excitación
periódica, donde los estados se indexan como (α, n₂, n₄) con α el estado
base y n₂, n₄ los índices de fotones.
"""

import numpy as np
from sympy import Matrix
from tqdm import tqdm


def SW_transform_2(H, A_idcs, B_idcs, show_progress=True):
    """
    Transformación de Schrieffer-Wolff de segundo orden (estática).
    
    Calcula la corrección de segundo orden al Hamiltoniano efectivo
    proyectado al subespacio A, integrando los estados del subespacio B.
    
    Parameters
    ----------
    H : np.ndarray
        Hamiltoniano completo (puede contener símbolos de SymPy)
    A_idcs : list
        Índices de los estados del subespacio A (baja energía)
    B_idcs : list
        Índices de los estados del subespacio B (alta energía)
    show_progress : bool
        Mostrar barra de progreso
        
    Returns
    -------
    Matrix
        Matriz de SymPy con la corrección de segundo orden en el subespacio A
        
    Example
    -------
    >>> H_sw = SW_transform_2(H, [0, 1, 2, 3], [4, 5])
    """
    H_sw = np.zeros(H.shape, dtype=object)
    iterator = tqdm(range(len(A_idcs))) if show_progress else range(len(A_idcs))
    
    for m in iterator:
        for mp in range(len(A_idcs)):
            for l in range(len(B_idcs)):
                Efact_m = 1/(H[A_idcs[m], A_idcs[m]] - H[B_idcs[l], B_idcs[l]])
                Efact_mp = 1/(H[A_idcs[mp], A_idcs[mp]] - H[B_idcs[l], B_idcs[l]])
                Efact = Efact_m + Efact_mp
                H_sw[A_idcs[m], A_idcs[mp]] += 1/2 * H[A_idcs[m], B_idcs[l]] * H[B_idcs[l], A_idcs[mp]] * Efact
 
    # Retornar solo la matriz en la base A
    H_sw_A = np.zeros((len(A_idcs), len(A_idcs)), dtype=object)
    for m in range(len(A_idcs)):
        for mp in range(len(A_idcs)):
            H_sw_A[m, mp] = H_sw[A_idcs[m], A_idcs[mp]]
    
    return Matrix(H_sw_A)


class FloquetSW:
    """
    Transformación de Schrieffer-Wolff en el espacio de Floquet.
    
    Esta clase implementa las transformaciones de SW para Hamiltonianos
    periódicos en el tiempo, donde los estados se indexan como (α, n₂, n₄)
    con α el estado base y n₂, n₄ los índices de modos de Fourier.
    
    Parameters
    ----------
    H0 : np.ndarray
        Hamiltoniano estático (6x6 típicamente para punto cuántico doble)
    V2 : np.ndarray
        Matriz de acoplamiento para frecuencia ω₂
    V4 : np.ndarray
        Matriz de acoplamiento para frecuencia ω₄
    w2 : Symbol
        Símbolo de SymPy para la frecuencia ω₂
    w4 : Symbol
        Símbolo de SymPy para la frecuencia ω₄
    privileged : list, optional
        Índices de estados "privilegiados" que no adquieren energía de fotones
        (típicamente estados de carga [4, 5])
    """
    
    def __init__(self, H0, V2, V4, w2, w4, privileged=None):
        self.H0 = H0
        self.V2 = V2
        self.V4 = V4
        self.w2 = w2
        self.w4 = w4
        self.privileged = privileged if privileged is not None else []
    
    def HFp_mmp(self, m, mp):
        """
        Elemento de matriz del Hamiltoniano de Floquet perturbado.
        
        Calcula H_F'[m, m'] donde m = (α, n₂, n₄) es un índice compuesto.
        Solo devuelve elementos off-diagonal y acoplamientos por fotones.
        
        Parameters
        ----------
        m : list
            Índice [α, n₂, n₄] del bra
        mp : list  
            Índice [α', n₂', n₄'] del ket
            
        Returns
        -------
        Elemento de matriz (puede ser simbólico)
        """
        alpha, n2, n4 = m
        beta, k2, k4 = mp
        
        # Acoplamientos por H0 (sin cambio de fotones)
        if n2 == k2 and n4 == k4 and alpha != beta:
            return self.H0[alpha, beta]
        
        # Acoplamientos por V2 (±1 fotón en modo 2)
        if n4 == k4 and (n2 == k2 + 1 or n2 == k2 - 1):
            return self.V2[alpha, beta] / 2
        
        # Acoplamientos por V4 (±1 fotón en modo 4)
        if n2 == k2 and (n4 == k4 + 1 or n4 == k4 - 1):
            return self.V4[alpha, beta] / 2
        
        return 0
    
    def E_m(self, m):
        """
        Energía diagonal del estado m = (α, n₂, n₄).
        
        Para estados privilegiados (típicamente carga), la energía no
        depende de los índices de fotones.
        """
        alpha, n2, n4 = m
        if alpha in self.privileged:
            return self.H0[alpha, alpha]
        return (n2 * self.w2 + n4 * self.w4) + self.H0[alpha, alpha]
    
    def diff_E_mmp(self, m, mp):
        """
        Diferencia de energías E_m - E_m'.
        
        Maneja correctamente los estados privilegiados.
        """
        alpha, n2, n4 = m
        beta, k2, k4 = mp
        alpha_priv = alpha in self.privileged
        beta_priv = beta in self.privileged
        
        if alpha_priv and not beta_priv:
            return self.H0[alpha, alpha]
        if beta_priv and not alpha_priv:
            return -self.H0[beta, beta]
        
        return (self.H0[alpha, alpha] - self.H0[beta, beta] + 
                (n2 - k2) * self.w2 + (n4 - k4) * self.w4)
    
    def SWF_2(self, m, mp, ls, show_progress=True):
        """
        Corrección de Schrieffer-Wolff de segundo orden en Floquet.
        
        Calcula el elemento de matriz efectivo entre los estados m y m'
        al segundo orden en teoría de perturbaciones.
        
        Parameters
        ----------
        m : list
            Estado [α, n₂, n₄] inicial
        mp : list
            Estado [α', n₂', n₄'] final
        ls : list
            Lista de estados intermedios [γ, n₂, n₄]
        show_progress : bool
            Mostrar barra de progreso
            
        Returns
        -------
        Expresión simbólica de la corrección
        """
        A = 0
        iterator = tqdm(ls) if show_progress else ls
        
        for l in iterator:
            num = self.HFp_mmp(m, l) * self.HFp_mmp(l, mp)
            den1 = self.diff_E_mmp(m, l)
            den2 = self.diff_E_mmp(mp, l)
            A += num * (1 / den1 + 1 / den2)
        
        return A / 2
    
    def SWF_3(self, m, mp, lsA, lsB, mppsA, mppsB, show_progress=True):
        """
        Corrección de Schrieffer-Wolff de tercer orden en Floquet.
        
        Calcula el elemento de matriz efectivo al tercer orden, incluyendo
        tres contribuciones: A (intermedios en A), B (intermedios en B),
        y C (mixto A-B).
        
        Parameters
        ----------
        m : list
            Estado [α, n₂, n₄] inicial
        mp : list
            Estado [α', n₂', n₄'] final
        lsA : list
            Estados intermedios del tipo A
        lsB : list
            Estados intermedios del tipo B
        mppsA : list
            Estados m'' del tipo A
        mppsB : list
            Estados m'' del tipo B
        show_progress : bool
            Mostrar barra de progreso
            
        Returns
        -------
        Expresión simbólica de la corrección de tercer orden
        """
        A = 0
        B = 0
        C = 0
        
        # Término A: intermedios tipo A
        iterator_A = tqdm(lsA) if show_progress else lsA
        for l in iterator_A:
            for mpp in mppsA:
                num = self.HFp_mmp(m, l) * self.HFp_mmp(l, mpp) * self.HFp_mmp(mpp, mp)
                den = self.diff_E_mmp(mp, l) * self.diff_E_mmp(mpp, l)
                A += num / den
        A = -A / 2
        
        # Término B: intermedios tipo B
        iterator_B = tqdm(lsB) if show_progress else lsB
        for l in iterator_B:
            for mpp in mppsB:
                num = self.HFp_mmp(m, l) * self.HFp_mmp(l, mpp) * self.HFp_mmp(mpp, mp)
                den = self.diff_E_mmp(mp, l) * self.diff_E_mmp(mpp, l)
                B += num / den
        B = -B / 2
        
        # Término C: mixto A-B
        iterator_C = tqdm(lsA) if show_progress else lsA
        for l in iterator_C:
            for lp in lsB:
                num = self.HFp_mmp(m, l) * self.HFp_mmp(l, lp) * self.HFp_mmp(lp, mp)
                den1 = self.diff_E_mmp(m, l) * self.diff_E_mmp(m, lp)
                den2 = self.diff_E_mmp(mp, l) * self.diff_E_mmp(mp, lp)
                C += num * (1 / den1 + 1 / den2)
        C = C / 2
        
        return A + B + C


def generate_ls_indices(base_n2n4_list, m, mp, n_states=6):
    """
    Genera lista de estados intermedios para sumas de SW.
    
    Parameters
    ----------
    base_n2n4_list : list
        Lista de pares [n₂, n₄] base
    m : list
        Estado m = [α, n₂, n₄] a excluir
    mp : list
        Estado m' = [α', n₂', n₄'] a excluir
    n_states : int
        Número de estados base (típicamente 6)
        
    Returns
    -------
    list
        Lista de estados [γ, n₂, n₄] intermedios válidos
    """
    ls = []
    for gamma in range(n_states):
        for n2, n4 in base_n2n4_list:
            l = [gamma, n2, n4]
            if l != m and l != mp:
                ls.append(l)
    return ls
