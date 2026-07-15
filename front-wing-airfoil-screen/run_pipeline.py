#!/usr/bin/env python3
"""
Run the whole front-wing pipeline, in order, in one command.

WHY THIS EXISTS
---------------
The stages have a hard dependency chain: nothing downstream can run until
airfoil_screen.py has written results.xlsx, run_config.json and
kulfan_params.json. A VS Code "compound" launch configuration cannot express
that - compounds start every configuration in PARALLEL, so compare.py would
fire before the screen had produced anything for it to read.

This driver runs the stages sequentially, fails fast on the first error, and
prints timings. It is also the thing to point a debugger at when you want to
reproduce a whole run end to end.

THE CHAIN
---------
  0. aero_targets.py       vehicle load budget -> what the front wing must make.
                           Independent of the screen; run first because it is
                           what justifies the CL targets the screen then uses.
  1. airfoil_screen.py     screen + gate + rank        -> results.xlsx
  2. compare.py            side-by-side sheets         -> plots/compare_*.png
  3. robustness.py         transition sweep            -> plots/robustness_*.png
  4. weight_sensitivity.py is the ranking weight-driven? -> plots/weight_*.png
  5. mesh_section.py       2-D RANS mesh in ground effect -> mesh_out/

All five run by default. gmsh (stage 5) is in requirements.txt, so there is
nothing extra to install - and a CFD mesh nobody builds is a CFD mesh nobody
knows is broken. `--skip-mesh` opts out; `--quick` implies it, because a
150-airfoil smoke run's rank-1 sections are meaningless and meshing them just
wastes a minute.

Usage:
  python run_pipeline.py                       # everything -> screen_results/
  python run_pipeline.py --quick               # 150-airfoil smoke test (no mesh)
  python run_pipeline.py --skip-mesh           # stop after stage 4
  python run_pipeline.py --results out_dir     # somewhere else
  python run_pipeline.py --only 2 3            # re-run just some stages
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# (id, label, script, arg-builder). The arg-builder takes the parsed args and
# returns the CLI arguments for that stage - note the screen uses --out while
# everything downstream uses --results.
STAGES = [
    (0, "Vehicle load budget", "aero_targets.py",
     lambda a: []),
    (1, "Airfoil screen", "airfoil_screen.py",
     lambda a: ["--out", a.results] + (["--quick"] if a.quick else [])
               + (["--no-variants"] if a.no_variants else [])
               + (["--seed", str(a.seed)])),
    (2, "Comparison sheets", "compare.py",
     lambda a: ["--results", a.results]),
    (3, "Transition robustness", "robustness.py",
     lambda a: ["--results", a.results]),
    (4, "Weight sensitivity", "weight_sensitivity.py",
     lambda a: ["--results", a.results]),
    (5, "XFoil validation", "validate.py",
     lambda a: ["--results", a.results, "--element", "all",
                "--skip-if-no-xfoil"]),
    (6, "CFD mesh (ground effect)", "mesh_section.py",
     lambda a: ["--results", a.results]),
]


def run(stage_id, label, script, argv):
    cmd = [sys.executable, os.path.join(HERE, script)] + argv
    bar = "=" * 72
    print(f"\n{bar}\n[{stage_id}] {label}\n    {script} {' '.join(argv)}\n{bar}",
          flush=True)
    t0 = time.time()
    # Inherit stdout/stderr so tqdm bars and prints stream live, and so a
    # debugger attached to this process still shows the child's output.
    proc = subprocess.run(cmd, cwd=HERE)
    dt = time.time() - t0
    if proc.returncode != 0:
        print(f"\n!! [{stage_id}] {label} FAILED (exit {proc.returncode}) "
              f"after {dt:.0f}s", flush=True)
        if script == "mesh_section.py":
            print("   The CFD mesh stage needs gmsh. It is in requirements.txt:")
            print("     pip install -r requirements.txt")
            print("   If gmsh installs but fails to IMPORT with a missing")
            print("   libGLU.so.1 / libXcursor, you are on Linux and need the")
            print("   system GL libs too:")
            print("     sudo apt-get install libglu1-mesa libxcursor1 libxinerama1")
            print("   (Windows and macOS wheels bundle these - no extra step.)")
            print("   Or skip it:  python run_pipeline.py --skip-mesh")
        sys.exit(proc.returncode)
    print(f"[{stage_id}] {label} ok  ({dt:.0f}s)", flush=True)
    return dt


def main():
    ap = argparse.ArgumentParser(
        description="Run the front-wing pipeline end to end.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="screen_results",
                    help="output dir (screen writes it, the rest read it)")
    ap.add_argument("--quick", action="store_true",
                    help="150-airfoil smoke test - use this to debug the chain")
    ap.add_argument("--no-variants", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-mesh", action="store_true",
                    help="do not build the CFD mesh (stage 6)")
    ap.add_argument("--only", type=int, nargs="+", default=None,
                    help="run only these stage ids, e.g. --only 2 3 4")
    args = ap.parse_args()

    MESH_STAGE = 6
    # The mesh stage runs by DEFAULT. gmsh is in requirements.txt, and a CFD
    # mesh that nobody builds is a CFD mesh nobody knows is broken. --skip-mesh
    # opts out; --quick implies it, because a 150-airfoil smoke run's rank-1
    # sections are meaningless and meshing them wastes a minute.
    skip_mesh = args.skip_mesh or args.quick
    wanted = [s for s in STAGES
              if (args.only is None or s[0] in args.only)
              and (s[0] != MESH_STAGE or not skip_mesh
                   or (args.only and MESH_STAGE in args.only))]

    print("FRONT-WING PIPELINE")
    print(f"  results dir : {args.results}")
    print(f"  mode        : {'QUICK smoke test' if args.quick else 'full run'}")
    print(f"  stages      : {', '.join(f'{s[0]}:{s[1]}' for s in wanted)}")
    if skip_mesh and args.only is None:
        why = "--quick implies it" if args.quick else "--skip-mesh"
        print(f"  (CFD mesh stage skipped: {why})")

    total = 0.0
    times = []
    for sid, label, script, mkargs in wanted:
        dt = run(sid, label, script, mkargs(args))
        times.append((sid, label, dt))
        total += dt

    print("\n" + "=" * 72)
    print("PIPELINE COMPLETE")
    print("=" * 72)
    for sid, label, dt in times:
        print(f"  [{sid}] {label:<28s} {dt:6.0f}s")
    print(f"  {'TOTAL':<33s}{total:6.0f}s")
    print(f"\n  results  -> {args.results}/results.xlsx")
    print(f"  figures  -> {args.results}/plots/")
    if args.quick:
        print("\n  NOTE: --quick screens only 150 airfoils. The rankings from a")
        print("  quick run are NOT valid - it exists to test that the chain runs.")


if __name__ == "__main__":
    main()
