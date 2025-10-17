
import argparse, os, yaml
import numpy as np
from scipy.optimize import minimize

from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator as Estimator

from src.qiskit_dmet_vqe.integrals_to_problem import problem_from_integrals


def problem_from_config(sys_cfg_path: str):
    with open(sys_cfg_path, 'r') as f:
        s = yaml.safe_load(f)['system']
    atom = open(s['xyz_file']).read() if os.path.exists(s['xyz_file']) else s['xyz_file']
    driver = PySCFDriver(
        atom=atom,
        basis=s.get('basis', 'sto3g'),
        charge=int(s.get('charge', 0)),
        spin=int(s.get('spin', 0)),
        unit=DistanceUnit.ANGSTROM if str(s.get('unit', 'angstrom')).lower().startswith('ang') else DistanceUnit.BOHR
    )
    return driver.run()


def vqe_energy_v2(problem, maxiter=200, no_taper=False):
    # 若使用者沒有提供 num_particles（如 integrals 模式），先嘗試建不帶粒子數的 mapper
    mapper = ParityMapper(num_particles=getattr(problem, 'num_particles', None))
    if no_taper:
        tmapper = mapper
    else:
        try:
            tmapper = problem.get_tapered_mapper(mapper)
        except Exception as e:
            print(f"[warn] tapering disabled: {e}")
            tmapper = mapper

    ferm_op = problem.hamiltonian.second_q_op()
    qubit_op: SparsePauliOp = tmapper.map(ferm_op)
    e_nuc = float(problem.hamiltonian.nuclear_repulsion_energy)

    n_orb = problem.num_spatial_orbitals
    n_part = getattr(problem, 'num_particles', None)

    
    # ✅ Debug 輸出
    n_orb = problem.num_spatial_orbitals
    n_part = getattr(problem, 'num_particles', None)
    print(f"[dbg] num_spatial_orbitals={n_orb}, num_particles={n_part}")

    # 若未設置 num_particles，UCCSD 仍可在某些映射下工作，但建議指定 --nelec
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

    res = minimize(energy, x0=theta0, method="COBYLA", options={"maxiter": maxiter})
    return float(res.fun)


def main():
    ap = argparse.ArgumentParser(description='VQE (V2 primitives) on a fragment: driver or integrals mode.')
    ap.add_argument('--mode', choices=['driver','integrals'], required=True)
    ap.add_argument('--config', help='configs/system.yaml (for driver mode)')
    ap.add_argument('--h1', help='h1.npy (for integrals mode)')
    ap.add_argument('--h2', help='h2.npy (for integrals mode)')
    ap.add_argument('--enuc', help='enuc.npy (for integrals mode)')
    ap.add_argument('--nelec', help="總電子數 (例如 2) 或 'nalpha,nbeta' (例如 1,1)")
    ap.add_argument('--no-taper', action='store_true', help='關閉 Z2 tapering')
    ap.add_argument('--maxiter', type=int, default=200)
    ap.add_argument('--out', default='results/energy.txt')
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.mode == 'driver':
        assert args.config, '--config is required for driver mode'
        problem = problem_from_config(args.config)
    else:
        assert args.h1 and args.h2 and args.enuc, 'integrals mode requires --h1 --h2 --enuc'
        problem = problem_from_integrals(args.h1, args.h2, args.enuc)
        # 🔧 明確指定空間軌域數（依 h1 的維度）
        import numpy as np
        n_orb_from_h1 = int(np.load(args.h1).shape[0])
        problem.num_spatial_orbitals = n_orb_from_h1


    # 如果使用者提供 --nelec，就設定 problem.num_particles
    if args.nelec:
        if ',' in args.nelec:
            na, nb = map(int, args.nelec.split(','))
        else:
            ne = int(args.nelec)
            na = ne // 2
            nb = ne - na
        problem.num_particles = (na, nb)

    E = vqe_energy_v2(problem, maxiter=args.maxiter, no_taper=args.no_taper)
    print('E_total (Ha):', E)
    with open(args.out, 'w') as f:
        f.write(str(E))

if __name__ == '__main__':
    main()
