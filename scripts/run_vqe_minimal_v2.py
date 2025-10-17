
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit.primitives import StatevectorEstimator as Estimator
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize
import numpy as np

def main():
    driver = PySCFDriver(
        atom="H 0 0 0; H 0 0 0.735",
        basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM
    )
    problem = driver.run()
    mapper = ParityMapper(num_particles=problem.num_particles)
    tmapper = problem.get_tapered_mapper(mapper)
    ferm_op = problem.hamiltonian.second_q_op()
    qubit_op: SparsePauliOp = tmapper.map(ferm_op)
    e_nuc = float(problem.hamiltonian.nuclear_repulsion_energy)

    n_orb = problem.num_spatial_orbitals
    n_part = problem.num_particles
    hf = HartreeFock(n_orb, n_part, tmapper)
    ansatz = UCCSD(num_spatial_orbitals=n_orb, num_particles=n_part,
                   qubit_mapper=tmapper, initial_state=hf)

    estimator = Estimator()
    theta0 = np.zeros(ansatz.num_parameters)

    def energy(theta):
        bound = ansatz.assign_parameters(theta, inplace=False)
        job = estimator.run([(bound, qubit_op, [])])
        pub = job.result()[0]
        evs = np.asarray(pub.data.evs)
        ev = float(evs if evs.ndim == 0 else evs.flat[0])
        return ev + e_nuc

    res = minimize(energy, x0=theta0, method="COBYLA", options={"maxiter": 200})
    print("E_total (Ha):", res.fun)

if __name__ == "__main__":
    main()
