#!/usr/bin/env python3
import os, sys, subprocess, random
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = Path(".")
DATA_DIR = BASE_DIR / "data"
OUT_BASE = BASE_DIR / "output"
PH_VALUE = 7.0

OUT_BASE.mkdir(exist_ok=True, parents=True)

class LigandSelect(Select):
    def accept_residue(self, residue):
        het = residue.id[0].strip()
        return bool(het and het not in ("W", "HOH"))

class ProteinSelect(Select):
    def accept_residue(self, residue):
        het = residue.id[0].strip()
        return bool(not het or het in ("W", "HOH"))

def run_cmd(cmd):
    print(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def detect_guanidinium_and_protonate(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            Ns = [n for n in atom.GetNeighbors() if n.GetSymbol() == "N"]
            if len(Ns) >= 2:
                Ns[0].SetFormalCharge(+1)
                print(f"  → Protonated guanidinium-like center at C {atom.GetIdx()}")
                return True
    return False

def generate_for_pdb(pdb_path: Path):
    name = pdb_path.stem
    outdir = OUT_BASE / name
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Processing {name} ===")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, str(pdb_path))
    ligand_pdb = outdir / "ligand.pdb"
    protein_pdb = outdir / "protein.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(ligand_pdb), LigandSelect())
    io.save(str(protein_pdb), ProteinSelect())
    print("✅ Split into ligand.pdb & protein.pdb")

    total_charge = 0
    try:
        mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=False)
        if mol is None:
            raise ValueError("RDKit failed to parse ligand.pdb")
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        mol = Chem.AddHs(mol)
        has_guanidinium = detect_guanidinium_and_protonate(mol)

        try:
            conf = mol.GetConformer()
            _ = conf.GetAtomPosition(0)
        except Exception:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.UFFOptimizeMolecule(mol)

        total_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        Chem.MolToMolFile(mol, str(outdir / "ligand.sdf"))
        print(f"✅ ligand.sdf generated (formal charge: {total_charge:+d}, guanidinium={has_guanidinium})")
    except Exception as e:
        print(f"⚠️ RDKit failed for {name}: {e}")

    protein_pqr = outdir / "protein.pqr"
    try:
        run_cmd(["pdb2pqr", f"--ff=amber", f"--with-ph={PH_VALUE}", str(protein_pdb), str(protein_pqr)])
    except subprocess.CalledProcessError as e:
        print(f"❌ pdb2pqr failed for {name}: {e}")
        return

    charges_txt = outdir / "charges.txt"
    n_atoms, q_total = 0, 0.0
    with open(protein_pqr) as f_in, open(charges_txt, "w") as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 9:
                    x, y, z, q = parts[5:9]
                    qf = float(q)
                    f_out.write(f"{float(x):10.3f} {float(y):10.3f} {float(z):10.3f} {qf:8.4f}\n")
                    q_total += qf
                    n_atoms += 1
    print(f"✅ charges.txt written ({n_atoms} atoms, total charge = {q_total:+.3f} e)")

    if abs(q_total) > 0.5:
        n_na = int(round(q_total))
        with open(charges_txt, "a") as f:
            for _ in range(abs(n_na)):
                x, y, z = random.uniform(-25,25), random.uniform(-25,25), random.uniform(-25,25)
                q = -1.0 if q_total > 0 else +1.0
                f.write(f"{x:10.3f}{y:10.3f}{z:10.3f}{q:8.4f}\n")
        print(f"⚠️ Added {abs(n_na)} Na+ ions for neutralization")

    q_now = 0.0
    with open(charges_txt) as f:
        for ln in f:
            if ln.strip():
                q_now += float(ln.split()[-1])
    if abs(q_now) > 1e-4:
        corr = -q_now
        with open(charges_txt, "a") as f:
            f.write(f"{50.000:10.3f}{50.000:10.3f}{50.000:10.3f}{corr:8.4f}\n")
        print(f"✅ Added correction point charge {corr:+.4f} e → total charge ≈ 0")
    else:
        print("✅ System neutral after ion placement")

    print(f"📊 Summary: ligand charge={total_charge:+d}, atoms={n_atoms}, final total charge={q_now:+.3f} e")

def main():
    pdb_files = sorted(DATA_DIR.glob("*.pdb"))
    if not pdb_files:
        print(f"❌ No PDB files found in {DATA_DIR}")
        return
    print(f"Found {len(pdb_files)} PDB complexes in {DATA_DIR}")
    for pdb in pdb_files:
        try:
            generate_for_pdb(pdb)
        except Exception as e:
            print(f"⚠️ Skipped {pdb.name}: {e}")
    print("\n🎉 Batch generation complete!")

if __name__ == "__main__":
    main()
