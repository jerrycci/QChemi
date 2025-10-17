from __future__ import annotations
from typing import Dict, Any, Tuple
import math
import random

def solve_fragment_hf_stub(fragment: Dict[str, Any], basis: str = "sto-3g") -> Tuple[float, Dict[str, Any]]:
    """Very cheap closed-form 'HF' on a fragment: returns a deterministic toy energy.

    Replace with: build 1e/2e integrals for the fragment (PySCF), run RHF/UHF.
    """
    n_orb = int(fragment.get("n_orb", 2))
    n_elec = int(fragment.get("n_elec", 2))
    # Simple convex function as a stand-in
    e = -0.3 * n_orb - 0.02 * (n_elec - n_orb/2.0)**2
    return e, {"method": "HF-stub", "basis": basis, "n_orb": n_orb, "n_elec": n_elec}

def solve_fragment_vqe_or_stub(fragment: Dict[str, Any],
                               basis: str = "sto-3g",
                               seed: int = 7,
                               use_vqe: bool = True) -> Tuple[float, Dict[str, Any]]:
    """Try a real Qiskit VQE on a tiny 2e/4so Hamiltonian if Qiskit is available; otherwise stub.

    The *API* is stable no matter what backend you have, so the example script works for everyone.
    """
    if not use_vqe:
        return _vqe_stub(fragment, seed=seed)

    try:
        # Lazy imports so that the file can be imported without qiskit installed
        from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
        from qiskit_nature.second_q.operators import ElectronicIntegrals
        from qiskit_nature.second_q.mappers import JordanWignerMapper
        from qiskit_algorithms.utils import algorithm_globals
        from qiskit_algorithms.optimizers import COBYLA
        from qiskit.primitives import Estimator
        from qiskit.circuit.library import TwoLocal
        import numpy as np

        algorithm_globals.random_seed = seed

        # Build a tiny 2e/4so toy Hamiltonian that vaguely resembles an active-space fragment
        h1 = np.array([[ -1.0, -0.1, 0.0, 0.0],
                       [ -0.1, -0.9, 0.0, 0.0],
                       [  0.0,  0.0, -0.5, -0.05],
                       [  0.0,  0.0, -0.05, -0.45]])
        h2 = np.zeros((4,4,4,4))
        # Minimal repulsion on first pair
        h2[0,0,0,0] = 0.6
        h2[1,1,1,1] = 0.58

        elec_ints = ElectronicIntegrals.from_raw_integrals(h1, h2)
        elec_h = ElectronicEnergy.from_raw_integrals(elec_ints)
        second_q_op = elec_h.second_q_op()

        mapper = JordanWignerMapper()
        qubit_op = mapper.map(second_q_op)

        ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cx", reps=1, insert_barriers=False)
        est = Estimator()
        optimizer = COBYLA(maxiter=50)

        # Parameter-shift-free manual loop (since VQE class signatures changed across Qiskit 2.x)
        import numpy as np
        params = np.zeros(ansatz.num_parameters)

        def energy(theta):
            job = est.run([(ansatz.bind_parameters(theta), qubit_op, ())])
            val = job.result().values[0] if hasattr(job.result(), "values") else job.result().values
            return float(val)

        from math import isfinite
        best_e = None
        best_x = params.copy()

        # A tiny COBYLA-like loop (optimizer.minimize requires SciPy in some stacks)
        step = 0.1
        for _ in range(60):
            trial = best_x + (np.random.rand(ansatz.num_parameters) - 0.5) * step
            e = energy(trial)
            if (best_e is None) or (e < best_e):
                best_e, best_x = e, trial
            step *= 0.98

        meta = {"method": "VQE", "basis": basis, "n_qubits": qubit_op.num_qubits, "ansatz": "TwoLocal(ry,cx)"}
        return float(best_e), meta

    except Exception as exc:
        # Fall back to stub if qiskit not present or errors occur
        return _vqe_stub(fragment, seed=seed, note=f"fallback: {type(exc).__name__}: {exc}")

def _vqe_stub(fragment: Dict[str, Any], seed: int = 7, note: str = "") -> Tuple[float, Dict[str, Any]]:
    random.seed(seed)
    n_orb = int(fragment.get("n_orb", 4))
    # Produces a slightly better energy than HF stub to mimic correlation
    e = -0.35 * n_orb - 0.03
    return e, {"method": "VQE-stub", "note": note, "n_orb": n_orb}
