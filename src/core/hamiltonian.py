"""
Módulo de construcción de Hamiltonianos con segunda cuantización.

Este módulo implementa el formalismo de segunda cuantización para construir
Hamiltonianos de sistemas de muchos cuerpos, especialmente útil para
puntos cuánticos con múltiples orbitales y espines.

Teoría:
-------
En segunda cuantización, los operadores fermiónicos de creación (c†) y
aniquilación (c) satisfacen las relaciones de anticonmutación:
    {c_i, c†_j} = δ_{ij}
    {c_i, c_j} = {c†_i, c†_j} = 0

El valor esperado ⟨vac|c†_1...c†_n c_m...c_1|vac⟩ se calcula usando
el teorema de Wick para operadores fermiónicos.

Uso:
----
El módulo permite definir espacios de Hilbert de una partícula y construir
bases de muchos cuerpos. Los Hamiltonianos se definen como sumas de productos
de operadores de creación y aniquilación.

Ejemplo:
--------
>>> spin = Space('spin', ['up', 'down'])
>>> site = Space('site', [1, 2])
>>> full_space = spin * site
>>> cs = full_space.creations()
>>> # Ahora cs contiene los operadores c†_{spin,site}
"""

from typing import List, Tuple, Dict, Any
from numbers import Number
from numpy import conjugate as scalar_conjugate, zeros, ndarray


class c_internal:
    """
    Representación interna de un operador fermiónico de creación o aniquilación.
    
    Esta clase representa un operador c o c† actuando en un sitio específico
    del espacio de Hilbert, caracterizado por un conjunto de números cuánticos.
    
    Attributes
    ----------
    create : bool
        True si es operador de creación (c†), False si es de aniquilación (c).
    site : dict
        Diccionario con los números cuánticos que identifican el sitio.
        Ejemplo: {'spin': 'up', 'orbital': 1}
    
    Methods
    -------
    d()
        Retorna el operador conjugado hermítico (c↔c†).
    """
    
    def __init__(self, create: bool = False, site: Dict = None):
        """
        Inicializa un operador fermiónico.
        
        Parameters
        ----------
        create : bool, optional
            True para creación, False para aniquilación (default: False).
        site : dict, optional
            Números cuánticos del sitio (default: {'spin':'up', 'site':'down'}).
        """
        self.create = create
        self.site = site if site is not None else {'spin': 'up', 'site': 'down'}
    
    def d(self):
        """
        Retorna el operador hermítico conjugado.
        
        Returns
        -------
        c_internal
            Operador con create negado (c → c†, c† → c).
        """
        return c_internal(not self.create, self.site)
    
    def __str__(self):
        """Representación en string del operador."""
        return str(self.create) + str(self.site)


def delta_site(a: c_internal, b: c_internal) -> bool:
    """
    Verifica si dos operadores actúan en el mismo sitio (δ de Kronecker).
    
    Compara los números cuánticos de dos operadores. Retorna True si todos
    los números cuánticos compartidos son iguales.
    
    Parameters
    ----------
    a : c_internal
        Primer operador.
    b : c_internal
        Segundo operador.
    
    Returns
    -------
    bool
        True si los operadores actúan en el mismo sitio.
    """
    if len(a.site) > len(b.site):
        for key in b.site.keys():
            if key not in a.site.keys() or a.site[key] != b.site[key]:
                return False
    else:
        for key in a.site.keys():
            if key not in b.site.keys() or a.site[key] != b.site[key]:
                return False
            
    return True


def split_c_list(c_list: List[c_internal]) -> Tuple[List[Tuple[c_internal, int]], List[Tuple[c_internal, int]]]:
    """
    Separa una lista de operadores en creación y aniquilación.
    
    Parameters
    ----------
    c_list : list of c_internal
        Lista de operadores fermiónicos.
    
    Returns
    -------
    creation_list : list of (c_internal, int)
        Lista de tuplas (operador de creación, índice original).
    annihilation_list : list of (c_internal, int)
        Lista de tuplas (operador de aniquilación, índice original).
    """
    creation_list = []
    annihilation_list = []
    
    for i, c_i in enumerate(c_list):
        if c_i.create:
            creation_list.append((c_i, i))
        else:
            annihilation_list.append((c_i, i))
    
    return creation_list, annihilation_list


def calc_vac(c_list: List[c_internal]):
    """
    Calcula el valor esperado en el vacío de un producto de operadores.
    
    Utiliza el teorema de Wick para calcular:
        ⟨vac|c_{i_1} c_{i_2} ... c_{i_n}|vac⟩
    
    El resultado es no-nulo solo si hay igual número de operadores de
    creación y aniquilación, y se puede emparejar cada aniquilación con
    una creación a su derecha que actúe en el mismo sitio.
    
    Parameters
    ----------
    c_list : list of c_internal
        Producto de operadores fermiónicos (orden de izquierda a derecha).
    
    Returns
    -------
    result : int
        Valor esperado (0, +1, o -1 por signos fermiónicos).
    
    Notes
    -----
    El cálculo se realiza recursivamente contrayendo pares de operadores
    usando la regla de anticonmutación y acumulando los signos fermiónicos.
    """
    result = 0
    creation_list, anhilation_list = split_c_list(c_list)
    
    # El número de operadores de creación debe igualar al de aniquilación
    if len(creation_list) != len(anhilation_list) or creation_list == []:
        return 0
    
    # Caso base: un solo par
    elif len(creation_list) == 1:
        if delta_site(creation_list[0][0], anhilation_list[0][0]) and creation_list[0][1] > anhilation_list[0][1]:
            return 1
        else:
            return 0
    
    # Caso recursivo: contraer el primer operador de aniquilación
    c_i, i = anhilation_list[0]
    for c_j, j in creation_list:
        if delta_site(c_i, c_j) and j > i:
            # Contracción con signo fermiónico (-1)^(|i-j|-1)
            result += calc_vac(c_list[:i] + c_list[i+1:j] + c_list[j+1:]) * (-1)**(abs(i-j)-1)
            
    return result


def conjugate(c_list: List[c_internal]) -> List[c_internal]:
    """
    Calcula el conjugado hermítico de un producto de operadores.
    
    Para (A·B·C)† = C†·B†·A†, se conjugan todos los operadores y
    se invierte el orden.
    
    Parameters
    ----------
    c_list : list of c_internal
        Producto de operadores.
    
    Returns
    -------
    list of c_internal
        Producto conjugado hermítico.
    """
    result = []
    for c_i in c_list:
        result.append(c_i.d())
    return result[::-1]


def calc_Hamiltonian(H: List[Tuple[Number, List[c_internal]]], 
                     basis: List[List[Tuple[Number, List[c_internal]]]]) -> ndarray:
    """
    Calcula la representación matricial del Hamiltoniano en una base dada.
    
    El Hamiltoniano se especifica como suma de términos, cada uno con un
    coeficiente y un producto de operadores. La base también se especifica
    como superposiciones de estados de Fock.
    
    Parameters
    ----------
    H : list of (Number, list of c_internal)
        Hamiltoniano como lista de tuplas (coeficiente, producto de operadores).
    basis : list of list of (Number, list of c_internal)
        Base de estados, cada estado es una superposición de estados de Fock.
    
    Returns
    -------
    result : ndarray
        Matriz del Hamiltoniano H_{ij} = ⟨basis_i|H|basis_j⟩.
    
    Example
    -------
    >>> # Hamiltoniano de hopping: H = -t(c†_1 c_2 + c†_2 c_1)
    >>> c1_dag = c_internal(True, {'site': 1})
    >>> c1 = c_internal(False, {'site': 1})
    >>> c2_dag = c_internal(True, {'site': 2})
    >>> c2 = c_internal(False, {'site': 2})
    >>> H = [(-1, [c1_dag, c2]), (-1, [c2_dag, c1])]
    """
    result = zeros((len(basis), len(basis)), dtype=object)
    
    for i, bra_i in enumerate(basis):
        for factor_i, bra_i_el in bra_i:
            # Conjugar el bra: ⟨ψ| = (|ψ⟩)†
            bra_i_el = conjugate(bra_i_el)
            factor_i = scalar_conjugate(factor_i)
            
            for j, ket_j in enumerate(basis):
                for factor_j, ket_j_el in ket_j:
                    for factor_k, H_k in H:
                        # ⟨bra|H|ket⟩ = ⟨vac|bra_ops · H_ops · ket_ops|vac⟩
                        result[i][j] += factor_i * factor_j * calc_vac(bra_i_el + H_k + ket_j_el) * factor_k
    
    return result


class Space:
    """
    Representa un espacio de Hilbert de una partícula.
    
    El espacio se define por un nombre (o lista de nombres) para los
    números cuánticos y los valores que pueden tomar.
    
    Attributes
    ----------
    name : list of str
        Nombres de los números cuánticos.
    span : list of list
        Valores posibles para cada número cuántico.
    
    Methods
    -------
    __mul__(other)
        Producto tensorial de espacios.
    dim()
        Número de números cuánticos.
    creations()
        Lista de operadores de creación para todos los estados.
    annihilations()
        Lista de operadores de aniquilación para todos los estados.
    
    Example
    -------
    >>> spin = Space('spin', ['up', 'down'])
    >>> orbital = Space('orbital', [1, 2, 3])
    >>> full = spin * orbital  # 6 estados
    """
    
    def __init__(self, name: str | List[str], span: List | List[List]):
        """
        Inicializa un espacio de Hilbert.
        
        Parameters
        ----------
        name : str or list of str
            Nombre(s) del/los número(s) cuántico(s).
        span : list
            Valores posibles. Si es lista de listas, cada sublista corresponde
            a una combinación de valores para los números cuánticos.
        """
        if not isinstance(name, list):
            name = [name]
        self.name = name
        
        for i, el in enumerate(span):
            if not isinstance(el, list):
                span[i] = [span[i]]
        self.span = span

    def __mul__(self, other):
        """
        Producto tensorial de dos espacios.
        
        Parameters
        ----------
        other : Space
            Otro espacio de Hilbert.
        
        Returns
        -------
        Space
            Espacio producto con todos los pares de estados.
        """
        new_span = [] 
        for ei in self.span:
            for ej in other.span:
                new_span.append(ei + ej)      
        new_space = Space(self.name + other.name, new_span)
        return new_space
    
    def dim(self):
        """
        Retorna el número de números cuánticos.
        
        Returns
        -------
        int
            Número de grados de libertad.
        """
        return len(self.name)
    
    def creations(self) -> List[c_internal]:
        """
        Genera operadores de creación para todos los estados.
        
        Returns
        -------
        list of c_internal
            Operadores c† para cada estado del espacio.
        """
        cs = []
        for el in self.span:
            site = dict(zip(self.name, el))
            cs.append(c_internal(True, site))
        return cs
    
    def annihilations(self) -> List[c_internal]:
        """
        Genera operadores de aniquilación para todos los estados.
        
        Returns
        -------
        list of c_internal
            Operadores c para cada estado del espacio.
        """
        cs = []
        for el in self.span:
            site = dict(zip(self.name, el))
            cs.append(c_internal(False, site))
        return cs
