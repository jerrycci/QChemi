"""qiskit_dmet_vqe: minimal skeleton for DMET+VQE binding-energy workflow.

This package provides a thin, testable surface that mirrors the Kirsopp et al. IJQC (2022)
pipeline: MM field -> DMET fragments -> (HF + VQE on [NH2–C–NH+]) -> energies -> ranking.

Each module is deliberately lightweight so you can swap in real chemistry back-ends later.
"""

__all__ = [
    "energies", "dmet", "qpu"
]
