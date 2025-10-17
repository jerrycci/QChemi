
import os, argparse, yaml, numpy as np
from pyscf import gto, scf, ao2mo

def main():
    ap = argparse.ArgumentParser(description='Generate full-system MO integrals via PySCF.')
    ap.add_argument('--config', default='configs/system.yaml')
    ap.add_argument('--h1', default='frag_h1_mo.npy')
    ap.add_argument('--h2', default='frag_h2_mo.npy')
    ap.add_argument('--enuc', default='frag_enuc.npy')
    ap.add_argument('--df', action='store_true', help='use density-fitting RHF')
    args = ap.parse_args()

    with open(args.config, 'r') as f:
        s = yaml.safe_load(f)['system']
    atom = open(s['xyz_file']).read() if os.path.exists(s['xyz_file']) else s['xyz_file']
    unit = s.get('unit', 'angstrom').lower()
    mol = gto.M(atom=atom, basis=s.get('basis', 'sto3g'),
                charge=int(s.get('charge', 0)), spin=int(s.get('spin', 0)),
                unit='Angstrom' if unit.startswith('ang') else 'Bohr')
    mf = (scf.RHF(mol).density_fit() if args.df else scf.RHF(mol)).run()

    C = mf.mo_coeff
    hcore_ao = mf.get_hcore()
    enuc = mol.energy_nuc()
    h1_mo = C.T @ hcore_ao @ C
    norb = C.shape[1]
    eri_mo_2idx = ao2mo.full(mol, C)
    h2_mo = ao2mo.restore(1, eri_mo_2idx, norb)

    np.save(args.h1, h1_mo)
    np.save(args.h2, h2_mo)
    np.save(args.enuc, np.array(float(enuc)))
    print(f"[OK] Saved: {args.h1}, {args.h2}, {args.enuc} | norb={norb} | DF={args.df}")

if __name__ == '__main__':
    main()
