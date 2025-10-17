
from __future__ import annotations
from typing import Sequence, Dict, Any
import numpy as np

def _infer_spin_scale(D_mo: np.ndarray) -> float:
    lam_max = np.linalg.eigvalsh(D_mo).max()
    return 2.0 if lam_max > 1.5 else 1.0

def build_bath_from_1rdm(dm1_ao: np.ndarray, C_loc: np.ndarray, frag_indices: Sequence[int], S_ao: np.ndarray, thresh: float=1e-6) -> Dict[str, Any]:
    """One-shot DMET bath from environment block of 1-RDM."""
    nmo = C_loc.shape[1]
    frag_idx = np.array(sorted(set(frag_indices)), dtype=int)
    env_idx = np.setdiff1d(np.arange(nmo), frag_idx)
    D_mo = C_loc.T @ (S_ao @ (dm1_ao @ (S_ao @ C_loc)))
    scale = _infer_spin_scale(D_mo)
    D_mo /= scale
    D_ee = D_mo[np.ix_(env_idx, env_idx)]
    evals, Ue = np.linalg.eigh(D_ee)
    frac_mask = (evals > thresh) & (evals < 1.0 - thresh)
    C_frag = C_loc[:, frag_idx]
    C_env  = C_loc[:, env_idx]
    C_bath = C_env @ Ue[:, frac_mask]
    # 🔧 保底：若沒有 fractional，就挑一個與 fragment 最耦合的環境軌域
    if C_bath.shape[1] == 0 and Ue.shape[1] > 0:
        # 用 D_mo 的 FE（fragment-env）塊作為耦合度量
        D_mo = C_loc.T @ (S_ao @ (dm1_ao @ (S_ao @ C_loc)))
        D_fe = D_mo[np.ix_(frag_idx, env_idx)]         # 形狀 (nfrag, nenv)
        j = np.argmax(np.linalg.norm(D_fe, axis=0))    # 取最相關的一個
        C_bath = C_env[:, [j]]

    return {
        'C_frag': C_frag,
        'C_bath': C_bath,
        'bath_eigs': evals[frac_mask],
        'env_frac_mask': frac_mask
    }
