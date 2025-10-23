#!/usr/bin/env python3
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 validate_charges_relaxed.py <ligand_dir>")
        sys.exit(1)

    case = Path(sys.argv[1])
    charges = case / "charges.txt"
    pqr = case / "protein.pqr"

    qs = []
    with open(charges) as f:
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.split()
            if len(parts) == 4:
                *xyz, q = map(float, parts)
                qs.append(q)

    with open(pqr) as f:
        pqr_atoms = [ln for ln in f if ln.startswith(("ATOM", "HETATM"))]

    diff = len(qs) - len(pqr_atoms)
    if diff > 0:
        print(f"Note: {diff} extra point charges (ions or correction) present.")
    elif diff < 0:
        print(f"⚠️ Missing charges for {abs(diff)} atoms.")
    else:
        print("✅ Atom counts consistent.")

    total_q = sum(qs)
    print(f"Total charge from charges.txt = {total_q:.4f} e")
    print("✅ charges.txt basic checks complete.")

if __name__ == "__main__":
    main()
