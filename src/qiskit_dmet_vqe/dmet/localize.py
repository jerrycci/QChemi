
from __future__ import annotations
from typing import Literal, Tuple
import numpy as np
from pyscf.lo import pipek, boys

def localize_orbitals(mf, method: Literal['pipek','boys']='pipek') -> Tuple[np.ndarray, np.ndarray]:
    """Localize occupied/virtual separately and concatenate.
    Returns (C_loc, occ_mask).
    """
    mol = mf.mol
    C = np.asarray(mf.mo_coeff)
    occ = np.asarray(mf.mo_occ) > 1e-8
    occ_idx = np.where(occ)[0]
    vir_idx = np.where(~occ)[0]
    if method == 'pipek':
        C_occ = pipek.PM(mol, C[:, occ_idx]).kernel()
        C_vir = pipek.PM(mol, C[:, vir_idx]).kernel()
    elif method == 'boys':
        C_occ = boys.Boys(mol, C[:, occ_idx]).kernel()
        C_vir = boys.Boys(mol, C[:, vir_idx]).kernel()
    else:
        raise ValueError('unknown localization method')
    C_loc = np.hstack([C_occ, C_vir])
    occ_mask = np.r_[np.ones(C_occ.shape[1], dtype=bool), np.zeros(C_vir.shape[1], dtype=bool)]
    return C_loc, occ_mask
