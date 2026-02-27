"""
validate_charges.py — Charge Method Sensitivity Analysis
=========================================================
Runs the same NH3 GCMC simulation on a subset of MOFs using both
PACMOF2 and MEPO-Qeq charge methods, then reports the difference.

This addresses the known validation gap: PACMOF2 predicts near-DFT DDEC6
charges but has not been benchmarked specifically on NH3 adsorption in MOFs.
This script produces a side-by-side comparison and flags MOFs where the
two methods disagree by more than a configurable threshold.

Usage:
    python validate_charges.py --input ./raw_nh3_core/cifs --max_mofs 10
    python validate_charges.py --input ./raw_nh3_core/cifs --threshold 0.5
"""

import json
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from nh3_simulation import nh3_uptake_pure


def run_both_methods(cif_file, rundir, rewrite_raspa_input=False):
    """
    Run pure NH3 GCMC with PACMOF2 and MEPO-Qeq on the same CIF.

    Returns
    -------
    dict with keys: uid, pacmof2, mepo_qeq, delta_mmol_g, delta_pct, flagged
    """
    uid = Path(cif_file).stem
    entry = {"uid": uid, "pacmof2": None, "mepo_qeq": None}

    for method in ("pacmof2", "mepo_qeq"):
        print(f"  [{uid}] Running {method}...")
        try:
            result = nh3_uptake_pure(
                cif_file,
                calc_charges=True,
                charge_method=method,
                rundir=str(Path(rundir) / method),
                rewrite_raspa_input=rewrite_raspa_input,
            )
            entry[method] = result
        except Exception as exc:
            print(f"  [{uid}] {method} failed: {exc}")
            entry[method] = {"error": str(exc)}

    return entry


def compute_comparison(entry, threshold=0.5):
    """
    Compare PACMOF2 vs MEPO-Qeq NH3 uptake for one MOF.

    threshold : float
        Flag the MOF if |delta| > threshold mmol/g.
    """
    p = entry.get("pacmof2") or {}
    m = entry.get("mepo_qeq") or {}

    p_val = p.get("NH3_uptake_mmol_g")
    m_val = m.get("NH3_uptake_mmol_g")

    if p_val is None or m_val is None:
        return {**entry, "delta_mmol_g": None, "delta_pct": None, "flagged": True}

    delta = p_val - m_val
    pct = 100.0 * delta / m_val if m_val != 0 else float("nan")
    flagged = abs(delta) > threshold

    return {**entry, "delta_mmol_g": round(delta, 5),
            "delta_pct": round(pct, 2), "flagged": flagged}


def print_summary(comparisons, threshold):
    """Print a human-readable comparison table."""
    header = f"{'MOF':<30} {'PACMOF2':>10} {'MEPO-Qeq':>10} {'Δ mmol/g':>10} {'Δ %':>8} {'Flag':>6}"
    print("\n" + "=" * len(header))
    print("CHARGE METHOD COMPARISON — NH3 UPTAKE @ 298 K, 1 bar")
    print(f"Flag threshold: |Δ| > {threshold} mmol/g")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    for c in comparisons:
        p_val = (c.get("pacmof2") or {}).get("NH3_uptake_mmol_g", "ERR")
        m_val = (c.get("mepo_qeq") or {}).get("NH3_uptake_mmol_g", "ERR")
        delta  = c.get("delta_mmol_g")
        pct    = c.get("delta_pct")
        flag   = "⚠" if c.get("flagged") else ""

        p_str = f"{p_val:.4f}" if isinstance(p_val, float) else str(p_val)
        m_str = f"{m_val:.4f}" if isinstance(m_val, float) else str(m_val)
        d_str = f"{delta:.4f}" if delta is not None else "N/A"
        pct_str = f"{pct:.1f}" if pct is not None else "N/A"

        print(f"{c['uid']:<30} {p_str:>10} {m_str:>10} {d_str:>10} {pct_str:>8} {flag:>6}")

    print("-" * len(header))
    flagged_n = sum(1 for c in comparisons if c.get("flagged"))
    deltas = [c["delta_mmol_g"] for c in comparisons if c.get("delta_mmol_g") is not None]
    if deltas:
        print(f"\nMean |Δ|: {np.mean(np.abs(deltas)):.4f} mmol/g  |  "
              f"Max |Δ|: {np.max(np.abs(deltas)):.4f} mmol/g  |  "
              f"MOFs flagged: {flagged_n}/{len(comparisons)}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Compare PACMOF2 vs MEPO-Qeq partial charges via NH3 GCMC uptake"
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Directory of CIF files to test")
    parser.add_argument("--max_mofs", type=int, default=10,
                        help="Number of MOFs to include (default: 10)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Flag threshold in mmol/g (default: 0.5)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: auto in input parent)")
    parser.add_argument("--rewrite_raspa_input", action="store_true",
                        help="Fix RASPA input line-endings (Windows)")
    args = parser.parse_args()

    input_dir = Path(args.input)
    cif_files = sorted(input_dir.glob("*.cif"))[: args.max_mofs]

    if not cif_files:
        print(f"No CIF files found in {input_dir}")
        sys.exit(1)

    print(f"Charge validation: {len(cif_files)} MOFs, threshold = {args.threshold} mmol/g\n")

    rundir = input_dir.parent / "gcmc_runs" / "charge_validation"
    rundir.mkdir(parents=True, exist_ok=True)

    raw = [run_both_methods(str(f), rundir, args.rewrite_raspa_input) for f in cif_files]
    comparisons = [compute_comparison(e, args.threshold) for e in raw]

    print_summary(comparisons, args.threshold)

    # Save JSON
    out_path = Path(args.output) if args.output else (
        input_dir.parent / "gcmc_runs" / "charge_validation_results.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(comparisons, f, indent=2)
    print(f"Results saved to: {out_path}")


if __name__ == "__main__":
    main()
