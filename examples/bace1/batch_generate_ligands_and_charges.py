#!/usr/bin/env python3
"""
Batch generator for ligand.sdf and charges.txt from multiple protein–ligand PDB complexes.
Implements the open-source equivalent of the MOE + AMBER10:EHT workflow (Kirsopp et al., 2022).
"""

import os, subprocess, sys
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem

# ============================================================
# Config
# ============================================================
BASE_DIR = Path(".")
DATA_DIR = BASE_DIR / "data"
OUT_BASE = BASE_DIR
PH_VALUE = 7.0

if not DATA_DIR.exists():
    sys.exit(f"❌ Data folder not found: {DATA_DIR}")

# ============================================================
# Selectors
# ============================================================
class LigandSelect(Select):
    def accept_residue(self, residue):
        het = residue.id[0].strip()
        return bool(het and het not in ("W", "HOH"))

class ProteinSelect(Select):
    def accept_residue(self, residue):
        het = residue.id[0].strip()
        return bool(not het or het in ("W", "HOH"))

# ============================================================
# Helper
# ============================================================
def run_cmd(cmd):
    print(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def generate_for_pdb(pdb_path: Path):
    """Generate ligand.sdf + charges.txt for one complex"""
    name = pdb_path.stem
    outdir = OUT_BASE / name
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Processing {name} ===")

    # --- Step 1: Split ligand & protein
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, str(pdb_path))

    ligand_pdb = outdir / "ligand.pdb"
    protein_pdb = outdir / "protein.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(ligand_pdb), LigandSelect())
    io.save(str(protein_pdb), ProteinSelect())
    print(f"✅ Split into ligand.pdb & protein.pdb")

    # --- Step 2: ligand.sdf via RDKit
    try:
        print("🔹 Generating ligand.sdf with RDKit...")
        mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
        if mol is None:
            raise ValueError("RDKit failed to parse ligand.pdb (no valid atoms)")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        Chem.MolToMolFile(mol, str(outdir / "ligand.sdf"))
        print(f"✅ ligand.sdf generated")
    except Exception as e:
        print(f"⚠️ RDKit failed for {name}: {e}")

    # --- Step 3: AMBER charges via pdb2pqr
    protein_pqr = outdir / "protein.pqr"
    try:
        run_cmd([
            "pdb2pqr",
            f"--ff=amber",
            f"--with-ph={PH_VALUE}",
            str(protein_pdb),
            str(protein_pqr)
        ])
    except subprocess.CalledProcessError as e:
        print(f"❌ pdb2pqr failed for {name}: {e}")
        return

    # --- Step 4: extract charges.txt
    charges_txt = outdir / "charges.txt"
    n_atoms = 0
    try:
        with open(protein_pqr) as f_in, open(charges_txt, "w") as f_out:
            for line in f_in:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) >= 9:
                        x, y, z, q = parts[5:9]
                        f_out.write(f"{float(x):10.3f} {float(y):10.3f} {float(z):10.3f} {float(q):8.4f}\n")
                        n_atoms += 1
        print(f"✅ charges.txt written ({n_atoms} atoms)")
    except Exception as e:
        print(f"⚠️ Failed writing charges.txt for {name}: {e}")

# ============================================================
# Main
# ============================================================
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
