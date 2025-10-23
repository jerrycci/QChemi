#!/usr/bin/env python3
from rdkit import Chem
import sys
from pathlib import Path

def has_headgroup(mol):
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        Ns = [n for n in a.GetNeighbors() if n.GetSymbol() == "N"]
        if len(Ns) >= 2 and any(n.GetFormalCharge() == 1 for n in Ns):
            return True
    return False

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 validate_sdf_topology.py <ligand_dir>")
        sys.exit(1)

    path = Path(sys.argv[1]) / "ligand.sdf"
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)

    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    print(f"Formal charge (SDF) = {charge}")
    print(f"Has NH2-C(=NH+)-NH2 headgroup? {'Yes' if has_headgroup(mol) else 'No'}")
    print("✅ ligand.sdf checks complete.")

if __name__ == "__main__":
    main()
