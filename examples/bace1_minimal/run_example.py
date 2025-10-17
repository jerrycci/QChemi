#!/usr/bin/env python3
"""Minimal end-to-end run for the skeleton binding-energy pipeline.

This script *runs without Qiskit installed* by using internal stubs.
If Qiskit 2.x + Nature 0.7 are available, the [NH2–C–NH+] fragment will try a tiny VQE.
"""
import json, os, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from qiskit_dmet_vqe.energies.binding_energy import compute_Ebind

ROOT = pathlib.Path(__file__).resolve().parent
lig_xyz = ROOT / "data" / "ligand.xyz"
prot_csv = ROOT / "data" / "protein_charges_ref.csv"
wat_xyz = ROOT / "data" / "waters_ref.xyz"

res = compute_Ebind(str(lig_xyz), str(prot_csv), str(wat_xyz), use_vqe=True, seed=7)
print("E(lig in protein) [Ha]:", res.e_ligand_in_protein)
print("E(lig in solvent) [Ha]:", res.e_ligand_in_solvent)
print("E_bind [Ha]:", res.e_bind)
print("\n-- Debug JSON (trimmed) --")
print(json.dumps({k: v for k, v in res.extras.items()}, indent=2)[:1200])
