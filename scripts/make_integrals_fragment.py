
import os, argparse, yaml, numpy as np
from pyscf import gto, scf, ao2mo
from qiskit_nature.units import DistanceUnit

from src.qiskit_dmet_vqe.dmet.localize import localize_orbitals
from src.qiskit_dmet_vqe.dmet.build_bath import build_bath_from_1rdm
from src.qiskit_dmet_vqe.dmet.fragmentify import project_to_fragment


def _infer_spin_scale(D: np.ndarray) -> float:
    lam_max = np.linalg.eigvalsh(D).max()
    return 2.0 if lam_max > 1.5 else 1.0


def main():
    ap = argparse.ArgumentParser(description='Generate fragment+bath MO integrals via PySCF+DMET.')
    ap.add_argument('--config', default='configs/system.yaml')
    ap.add_argument('--frag-indices', default='0', help='逗號分隔的 LMO 索引，例如 0 或 0,1')
    ap.add_argument('--h1', default='frag_h1_mo.npy')
    ap.add_argument('--h2', default='frag_h2_mo.npy')
    ap.add_argument('--enuc', default='frag_enuc.npy')
    ap.add_argument('--emit-nelec', default='frag_nelec.txt', help='輸出電子數的檔案 (nelec 與 nalpha,nbeta)')
    args = ap.parse_args()

    frag_idx = [int(x) for x in args.frag_indices.split(',') if x.strip()]

    # 讀幾何
    with open(args.config, 'r') as f:
        s = yaml.safe_load(f)['system']
    atom = open(s['xyz_file']).read() if os.path.exists(s['xyz_file']) else s['xyz_file']
    unit = s.get('unit', 'angstrom').lower()
    mol = gto.M(atom=atom, basis=s.get('basis', 'sto3g'),
                charge=int(s.get('charge', 0)), spin=int(s.get('spin', 0)),
                unit='Angstrom' if unit.startswith('ang') else 'Bohr')
    mf = scf.RHF(mol).run()

    S = mol.intor_symmetric('int1e_ovlp')
    dm1_ao = mf.make_rdm1()

    # 局域化
    #C_loc, occ_mask = localize_orbitals(mf, method='pipek')
    # DMET bath
    #bath = build_bath_from_1rdm(dm1_ao, C_loc, frag_idx, S, thresh=1e-6)
    # 1) 局域化（原本用 pipek，改用 boys 試試）
    C_loc, occ_mask = localize_orbitals(mf, method='boys')

    # 2) 建 bath 時，放寬 fractional 選擇（例如 1e-3）
    bath = build_bath_from_1rdm(dm1_ao, C_loc, frag_idx, S, thresh=1e-3)

    C_cluster = np.hstack([bath['C_frag'], bath['C_bath']])

    # AO 積分
    hcore_ao = mf.get_hcore()
    nao = mol.nao_nr()
    eri_ao = ao2mo.restore(1, ao2mo.full(mol, np.eye(nao)), nao)

    # 投影到 cluster-MO
    h1_mo, h2_mo = project_to_fragment(hcore_ao, eri_ao, C_cluster)

    # 估算 cluster 電子數（spin-summed → 轉成 spinless 再取跡）
    D_cluster = C_cluster.T @ (S @ (dm1_ao @ (S @ C_cluster)))
    scale = _infer_spin_scale(D_cluster)
    D_spinless = D_cluster / scale
    nelec_est = float(np.trace(D_spinless))
    # 預設封殼拆分
    nalpha = int(round(nelec_est/2))
    nbeta  = int(round(nelec_est - nalpha))

    # 儲存積分與核斥能
    np.save(args.h1, h1_mo)
    np.save(args.h2, h2_mo)
    np.save(args.enuc, np.array(float(mol.energy_nuc())))

    # 輸出電子數資訊
    with open(args.emit_nelec, 'w') as f:
        f.write(f"nelec={nelec_est}")
        f.write(f"nalpha={nalpha},nbeta={nbeta}")
    print(f"[OK] Saved: {args.h1}, {args.h2}, {args.enuc}; n_cluster={C_cluster.shape[1]}; nelec~{nelec_est:.6f} -> ({nalpha},{nbeta})")

if __name__ == '__main__':
    main()
