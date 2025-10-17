
from __future__ import annotations
import numpy as np
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.operators import ElectronicIntegrals
from qiskit_nature.second_q.problems import ElectronicStructureProblem

def problem_from_integrals(h1_path: str, h2_path: str, enuc_path: str) -> ElectronicStructureProblem:
    h1 = np.load(h1_path)                     # (norb, norb)
    h2 = np.load(h2_path)                     # (norb, norb, norb, norb)
    enuc = float(np.load(enuc_path))

    # 防呆（型別/維度）
    assert h1.ndim == 2 and h1.shape[0] == h1.shape[1], "h1 must be square (norb,norb)"
    assert h2.ndim == 4 and all(s == h1.shape[0] for s in h2.shape), "h2 must be (norb,norb,norb,norb)"

    integrals = ElectronicIntegrals.from_raw_integrals(h1, h2)
    ham = ElectronicEnergy(integrals)
    ham.nuclear_repulsion_energy = enuc

    problem = ElectronicStructureProblem(ham)

    # 🔧 明確設定空間軌域數（避免被誤判成 1）
    problem.num_spatial_orbitals = int(h1.shape[0])

    return problem
