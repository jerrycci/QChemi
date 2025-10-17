from __future__ import annotations
from typing import Dict, Iterable, Tuple

def parity(bitstring: str) -> int:
    """Return 0 if even number of 1s, else 1 (mod-2 parity)."""
    return sum(1 for b in bitstring if b == '1') & 1

def pmsv_filter(counts: Dict[str, int],
                symmetries: Iterable[Tuple[str, int]]) -> Dict[str, int]:
    """Partitioned measurement symmetry verification (PMSV) filter.

    Args:
        counts: mapping from bitstring -> shots
        symmetries: iterable of (mask, expected_parity), where mask is a bitmask-like string
                    of '0'/'1'/'X' (X=don't care). Parity is computed over bits where mask=='1'.

    Returns:
        Filtered counts dict with shots violating *any* symmetry removed.
    """
    def masked_parity(bits: str, mask: str) -> int:
        return sum(1 for b, m in zip(bits, mask) if m == '1' and b == '1') & 1

    out = {}
    for bits, shots in counts.items():
        ok = True
        for mask, expected in symmetries:
            if masked_parity(bits, mask) != (expected & 1):
                ok = False
                break
        if ok:
            out[bits] = out.get(bits, 0) + shots
    return out
