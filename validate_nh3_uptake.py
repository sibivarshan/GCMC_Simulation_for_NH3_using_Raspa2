"""
Validation Script: Compare GCMC NH3 uptake against raw_nh3_core dataset
========================================================================
Runs pure NH3 GCMC simulations on a subset of MOFs from the validation
dataset and compares simulated uptake against expected values.

Usage:
    python validate_nh3_uptake.py --n_mofs 3 --charge_method pacmof2
"""

import json
import argparse
import sys
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from nh3_simulation import nh3_uptake_pure


def load_validation_data(csv_path):
    """Load the raw_nh3_core validation dataset."""
    data = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            refcode = row["refcode"]
            uptake = float(row["NH3_uptake_298K_1bar [mmol/g]"])
            cif_path = row["cif_path"]
            data[refcode] = {"expected_uptake": uptake, "cif_path": cif_path}
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Validate NH3 GCMC against raw_nh3_core dataset"
    )
    parser.add_argument(
        "--dataset_dir", type=str,
        default=str(Path(__file__).parent / "raw_nh3_core"),
        help="Path to raw_nh3_core directory"
    )
    parser.add_argument(
        "--n_mofs", type=int, default=3,
        help="Number of MOFs to validate (default: 3)"
    )
    parser.add_argument(
        "--charge_method", type=str, default="pacmof2",
        choices=["pacmof2", "mepo_qeq"],
        help="Charge calculation method"
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output JSON file (default: auto)"
    )
    args = parser.parse_args()

    dataset_dir = Path(args.dataset_dir)
    csv_path = dataset_dir / "nh3_core_structured.csv"

    if not csv_path.exists():
        print(f"Error: Validation dataset not found at {csv_path}")
        sys.exit(1)

    # Load validation data
    val_data = load_validation_data(csv_path)
    print(f"Loaded {len(val_data)} MOFs from validation dataset")

    # Select a diverse subset: pick from low, medium, and high uptake ranges
    sorted_mofs = sorted(val_data.items(), key=lambda x: x[1]["expected_uptake"])

    if args.n_mofs >= 3:
        # Pick from low, middle, and high uptake ranges
        n = len(sorted_mofs)
        indices = [0, n // 2, n - 1]  # low, mid, high
        # Fill remaining with evenly spaced
        remaining = args.n_mofs - 3
        if remaining > 0:
            step = n // (remaining + 1)
            for i in range(1, remaining + 1):
                idx = i * step
                if idx not in indices:
                    indices.append(idx)
        indices = sorted(set(indices))[:args.n_mofs]
    else:
        indices = list(range(args.n_mofs))

    selected = [sorted_mofs[i] for i in indices]

    print(f"\nSelected {len(selected)} MOFs for validation:")
    for refcode, info in selected:
        print(f"  {refcode}: expected = {info['expected_uptake']:.4f} mmol/g")

    # Run GCMC simulations
    rundir = dataset_dir.parent / "gcmc_runs" / "validation"
    rundir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, (refcode, info) in enumerate(selected):
        cif_path = dataset_dir / info["cif_path"]
        if not cif_path.exists():
            print(f"\n[{i+1}/{len(selected)}] SKIP {refcode}: CIF not found at {cif_path}")
            continue

        print(f"\n[{i+1}/{len(selected)}] Simulating {refcode}...")
        print(f"  Expected uptake: {info['expected_uptake']:.4f} mmol/g")

        try:
            sim_result = nh3_uptake_pure(
                str(cif_path),
                calc_charges=True,
                charge_method=args.charge_method,
                rundir=str(rundir),
            )

            simulated = sim_result["NH3_uptake_mmol_g"]
            expected = info["expected_uptake"]
            error_pct = abs(simulated - expected) / expected * 100 if expected > 0 else float("inf")

            result = {
                "refcode": refcode,
                "expected_mmol_g": expected,
                "simulated_mmol_g": simulated,
                "error_pct": error_pct,
                "Qst_kJ_mol": sim_result.get("NH3_heat_of_adsorption_kJ_mol"),
                "status": "success",
            }
            print(f"  Simulated:  {simulated:.4f} mmol/g")
            print(f"  Error:      {error_pct:.1f}%")
            print(f"  Qst:        {result['Qst_kJ_mol']:.2f} kJ/mol" if result['Qst_kJ_mol'] else "")

        except Exception as e:
            result = {
                "refcode": refcode,
                "expected_mmol_g": info["expected_uptake"],
                "simulated_mmol_g": None,
                "error_pct": None,
                "status": f"error: {str(e)}",
            }
            print(f"  Error: {e}")

        results.append(result)

    # Save results
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = rundir / "validation_results.json"

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"\n{'=' * 60}")
    print("VALIDATION SUMMARY")
    print(f"{'=' * 60}")

    successful = [r for r in results if r["status"] == "success"]
    if successful:
        errors = [r["error_pct"] for r in successful if r["error_pct"] is not None]
        avg_error = sum(errors) / len(errors) if errors else float("inf")

        print(f"  MOFs tested:  {len(results)}")
        print(f"  Successful:   {len(successful)}")
        print(f"  Avg error:    {avg_error:.1f}%")
        print(f"\n  {'Refcode':<16} {'Expected':>10} {'Simulated':>10} {'Error':>8}")
        print(f"  {'-'*14}  {'-'*10} {'-'*10} {'-'*8}")
        for r in successful:
            print(
                f"  {r['refcode']:<16} {r['expected_mmol_g']:>10.4f} "
                f"{r['simulated_mmol_g']:>10.4f} {r['error_pct']:>7.1f}%"
            )
    else:
        print("  No successful simulations.")

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
