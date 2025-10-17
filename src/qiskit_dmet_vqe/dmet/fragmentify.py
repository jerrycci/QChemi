
from __future__ import annotations
import numpy as np

def transform_one_body(h1_ao: np.ndarray, C: np.ndarray) -> np.ndarray:
    return C.T @ (h1_ao @ C)

def transform_two_body(h2_ao: np.ndarray, C: np.ndarray) -> np.ndarray:
    g = h2_ao
    g = np.einsum('pqrs,pi->iqrs', g, C, optimize=True)
    g = np.einsum('iqrs,qj->ijrs', g, C, optimize=True)
    g = np.einsum('ijrs,rk->ijks', g, C, optimize=True)
    g = np.einsum('ijks,sl->ijkl', g, C, optimize=True)
    return g

def project_to_fragment(h1_ao: np.ndarray, h2_ao: np.ndarray, C_cluster: np.ndarray):
    return transform_one_body(h1_ao, C_cluster), transform_two_body(h2_ao, C_cluster)
