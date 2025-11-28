"""
Módulos core para teoría de Floquet en puntos cuánticos.

Este paquete contiene las funciones fundamentales para:
    - Construcción del Hamiltoniano de Floquet (monocromático y bicromático)
    - Cálculo del operador de evolución temporal
    - Transformaciones de Schrieffer-Wolff
    - Construcción de Hamiltonianos de segunda cuantización
"""

from .floquet import (
    floquet_bichrom,
    floquet_mono,
    idx,
    idx_mono,
    get_U,
    get_U_mono,
    get_P,
    get_P_mono,
)

from .hamiltonian import (
    c_internal,
    delta_site,
    split_c_list,
    calc_vac,
    conjugate,
    calc_Hamiltonian,
    Space,
)

from .schrieffer_wolff import (
    SW_transform_2,
    FloquetSW,
    generate_ls_indices,
)

# qutip es una dependencia opcional
try:
    from .quantum_dots import (
        f_destroy,
        f_create,
        eqdot_state,
        get_Lambda,
        get_Liouville,
        red_H_idx,
        red_H,
    )
    _HAS_QUTIP = True
except ImportError:
    _HAS_QUTIP = False

__all__ = [
    # floquet
    'floquet_bichrom',
    'floquet_mono',
    'idx',
    'idx_mono',
    'get_U',
    'get_U_mono',
    'get_P',
    'get_P_mono',
    # hamiltonian
    'c_internal',
    'delta_site',
    'split_c_list',
    'calc_vac',
    'conjugate',
    'calc_Hamiltonian',
    'Space',
    # schrieffer_wolff
    'SW_transform_2',
    'FloquetSW',
    'generate_ls_indices',
]

# Añadir exports de quantum_dots solo si qutip está disponible
if _HAS_QUTIP:
    __all__.extend([
        'f_destroy',
        'f_create',
        'eqdot_state',
        'get_Lambda',
        'get_Liouville',
        'red_H_idx',
        'red_H',
    ])
