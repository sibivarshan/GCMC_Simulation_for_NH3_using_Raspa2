"""
NH3 GCMC Screening Script
==========================
Entry point for screening MOFs for NH3 capture from humid air.
Adapted from MOFDiff's mofdiff/scripts/gcmc_screen.py.

Usage:
    # Pure NH3 adsorption (validation against raw_nh3_core):
    python gcmc_screen_nh3.py --input ./raw_nh3_core/cifs --mode pure --max_mofs 5

    # 3-component air removal with vacuum swing:
    python gcmc_screen_nh3.py --input ./path/to/cifs --mode air_removal --max_mofs 5

    # Using MEPO-Qeq charges instead of PACMOF2:
    python gcmc_screen_nh3.py --input ./path/to/cifs --charge_method mepo_qeq
"""

import json
import argparse
import sys
from pathlib import Path

# Add the GCMC directory to path
sys.path.insert(0, str(Path(__file__).parent))

from nh3_simulation import nh3_uptake_pure, working_capacity_nh3_air_removal


def screen_single_mof(args_tuple):
    """Screen a single MOF. Used for both sequential and parallel execution."""
    ciffile, mode, charge_method, rundir, rewrite_raspa_input = args_tuple
    uid = ciffile.stem

    try:
        if mode == "pure":
            result = nh3_uptake_pure(
                str(ciffile),
                calc_charges=True,
                charge_method=charge_method,
                rundir=str(rundir),
                rewrite_raspa_input=rewrite_raspa_input,
            )
        elif mode == "air_removal":
            result = working_capacity_nh3_air_removal(
                str(ciffile),
                calc_charges=True,
                charge_method=charge_method,
                rundir=str(rundir),
                rewrite_raspa_input=rewrite_raspa_input,
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")

        info = "success"
    except Exception as e:
        print(f"Error processing {ciffile}: {e}")
        result = None
        info = str(e)

    return {"uid": uid, "info": info, "results": result}


def main():
    parser = argparse.ArgumentParser(
        description="Screen MOFs for NH3 capture using GCMC simulations"
    )
    parser.add_argument(
        "--input", type=str, required=True,
        help="Directory containing CIF files to screen"
    )
    parser.add_argument(
        "--mode", type=str, default="air_removal",
        choices=["pure", "air_removal"],
        help="Simulation mode: 'pure' (validation) or 'air_removal' (3-component VSA)"
    )
    parser.add_argument(
        "--charge_method", type=str, default="pacmof2",
        choices=["pacmof2", "mepo_qeq"],
        help="Charge calculation method"
    )
    parser.add_argument(
        "--max_mofs", type=int, default=5,
        help="Maximum number of MOFs to screen (default: 5)"
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output JSON file path (default: auto-generated in input dir)"
    )
    parser.add_argument(
        "--rewrite_raspa_input", action="store_true",
        help="Rewrite RASPA input files to fix line-ending issues"
    )
    parser.add_argument(
        "--parallel", action="store_true",
        help="Use parallel processing (requires p_tqdm)"
    )
    parser.add_argument(
        "--ncpu", type=int, default=4,
        help="Number of CPUs for parallel processing"
    )
    args = parser.parse_args()

    input_dir = Path(args.input)
    if not input_dir.exists():
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    # Find CIF files
    all_files = sorted(input_dir.glob("*.cif"))
    if not all_files:
        print(f"Error: No CIF files found in {input_dir}")
        sys.exit(1)

    # Limit number of MOFs
    if args.max_mofs > 0:
        all_files = all_files[: args.max_mofs]

    print(f"=" * 60)
    print(f"NH3 GCMC Screening")
    print(f"=" * 60)
    print(f"  Mode:           {args.mode}")
    print(f"  Charge method:  {args.charge_method}")
    print(f"  MOFs to screen: {len(all_files)}")
    print(f"  Input:          {input_dir}")
    print(f"=" * 60)

    # Set up run directory
    rundir = input_dir.parent / "gcmc_runs" / args.mode
    rundir.mkdir(parents=True, exist_ok=True)

    # Prepare arguments for each MOF
    task_args = [
        (f, args.mode, args.charge_method, rundir, args.rewrite_raspa_input)
        for f in all_files
    ]

    # Run screening
    if args.parallel:
        try:
            from p_tqdm import p_umap
            results = p_umap(screen_single_mof, task_args, num_cpus=args.ncpu)
        except ImportError:
            print("Warning: p_tqdm not installed, falling back to sequential processing")
            results = [screen_single_mof(ta) for ta in task_args]
    else:
        results = []
        for i, ta in enumerate(task_args):
            print(f"\n[{i+1}/{len(task_args)}] Processing {ta[0].stem}...")
            results.append(screen_single_mof(ta))

    # Save results
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_dir.parent / "gcmc_runs" / f"nh3_{args.mode}_results.json"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"RESULTS SUMMARY")
    print(f"{'=' * 60}")
    successful = [r for r in results if r["info"] == "success"]
    failed = [r for r in results if r["info"] != "success"]
    print(f"  Successful: {len(successful)}/{len(results)}")
    print(f"  Failed:     {len(failed)}/{len(results)}")

    for r in successful:
        res = r["results"]
        if args.mode == "pure":
            print(
                f"  {r['uid']}: NH3 uptake = {res['NH3_uptake_mmol_g']:.4f} mmol/g"
            )
        elif args.mode == "air_removal":
            print(
                f"  {r['uid']}: WC = {res['working_capacity_nh3_mmol_g']:.4f} mmol/g, "
                f"NH3/N2 sel = {res.get('NH3_N2_selectivity', 'N/A')}"
            )

    if failed:
        print(f"\nFailed MOFs:")
        for r in failed:
            print(f"  {r['uid']}: {r['info']}")

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
