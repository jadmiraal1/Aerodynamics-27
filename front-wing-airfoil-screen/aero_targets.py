#!/usr/bin/env python3
"""
Front-wing load target from the vehicle downforce budget.

Answers the question the airfoil screen cannot: "what CL should the front wing
actually run at?" - by taking the car's total downforce and aero balance
targets and working backwards through the moment balance to the front wing.

WHY THIS FILE EXISTS
--------------------
`main_center`'s cl_target used to be 1.10, justified by a code comment reading
"lightly loaded: protect the undertray inlet". It had no derivation. Meanwhile
it silently determined which airfoil won: sweeping the target from 0.5 to 1.3
produced a completely different winning section at every step, with the
winner's camber tracking the target almost linearly. An undocumented constant
was choosing the airfoil.

THE TRAP THIS FILE EXISTS TO PREVENT
------------------------------------
"Front downforce distribution = 40%" is the fraction of downforce carried by
the FRONT AXLE. It is NOT the fraction produced by the FRONT WING. Those differ
by a large factor, because devices act through moment arms:

  - The front wing sits AHEAD of the front axle. Every newton it makes puts
    MORE than a newton on the front axle (~1.5x) and LIFTS the rear axle.
  - The rear wing sits BEHIND the rear axle, so it actively UNLOADS the front.
  - The undertray sits between the axles and splits roughly by position.

Size the front wing to make 40% of total downforce and you build a wing about
twice as large as you need, and push the car into deep understeer at speed.

WHAT IS DERIVED VS ASSUMED
--------------------------
DERIVED (given the inputs): device load factors, front wing share, F_fw, ClA_fw,
CL_fw, and the center/outboard CL split.

ASSUMED (edit these - they are the whole answer): everything in VEHICLE below.
The undertray share in particular swings the front wing target about as much as
the entire front-balance range does. It is on the critical path.

NOT MODELED: ground effect, wake interaction between devices, CoP migration with
ride height/pitch, and the difference between a 2-D section CL and the 3-D wing
CL. See the caveats printed at the end, and DESIGN_JUSTIFICATION.md.

Usage:
  python aero_targets.py                       # nominal + sensitivity
  python aero_targets.py --cla 3.5 --beta 0.40 --undertray 0.50
"""

import argparse

import numpy as np

RHO = 1.225                     # air density [kg/m^3], sea level ISA

# ---------------------------------------------------------------------------
# VEHICLE ASSUMPTIONS - EDIT THESE. Every number below is an input, not a
# result. Each carries a provenance tag so it can be defended (or challenged).
# ---------------------------------------------------------------------------
VEHICLE = dict(
    # --- targets (from the team's aero goals) ---
    CLA_total=4.0,          # [m^2] TARGET total ClA.            <- TEAM TARGET
    beta_front=0.40,        # front AXLE share of downforce.     <- TEAM TARGET
                            #    (0.35-0.45 band; see sensitivity below)

    # --- downforce split between devices ---
    undertray_share=0.45,   # undertray's share of TOTAL downforce.  <- ASSUMED
                            #    THE DOMINANT UNKNOWN. Pin this down.

    # --- geometry [m] ---
    wheelbase=1.55,         # <- ASSUMED. Check against this year's chassis.
    a_front_wing=0.73,      # front wing CoP AHEAD of the front axle.
                            #    Derived from the packaging study:
                            #    bulkhead->front wheel = 21.5 in = 0.546 m, and
                            #    the wing's CoP sits ~0.19 m ahead of the
                            #    bulkhead (roughly 40% back along the 0.317 m
                            #    aero package). <- SEMI-DERIVED
    x_undertray=0.90,       # undertray CoP BEHIND the front axle.   <- ASSUMED
    x_rear_wing=1.85,       # rear wing CoP BEHIND the front axle.   <- ASSUMED
                            #    (i.e. 0.30 m behind the rear axle)

    # --- front wing planform ---
    S_front_wing=0.44,      # [m^2] front wing reference area.       <- ASSUMED
                            #    ~1.40 m span x 0.317 m packaging chord.
    center_span_frac=0.30,  # fraction of span that is main-plane-only (center,
                            #    feeding the undertray).             <- CONCEPT

    # --- operating point ---
    v_design=20.0,          # [m/s] design speed, from the TR26 lap sim
                            #    aero-weighted speed distribution.   <- DERIVED
)


def load_factors(v):
    """Front-axle load produced per newton of each device's downforce.

    For a device at longitudinal position x (measured rearward from the front
    axle, so x < 0 means AHEAD of it), a static moment balance about the rear
    axle gives its share of the front-axle load as (L - x)/L.

      front wing:  x = -a  ->  (L + a)/L  > 1   (loads front, LIFTS rear)
      undertray:   0 < x < L  ->  0..1           (splits by position)
      rear wing:   x > L      ->  negative       (UNLOADS the front axle)
    """
    L = v["wheelbase"]
    return {
        "front_wing": (L + v["a_front_wing"]) / L,
        "undertray": (L - v["x_undertray"]) / L,
        "rear_wing": (L - v["x_rear_wing"]) / L,
    }


def solve_front_wing(v):
    """Front wing's share of total downforce, from the moment balance.

    Two constraints, three device forces (as fractions f, u, r of total):

        f + u + r = 1                                    (they sum to the total)
        k_f*f + k_u*u + k_r*r = beta                     (they produce the balance)

    Underdetermined by one, so the undertray share `u` is taken as an input -
    it is the device whose size is most constrained by the chassis, and the
    one you can measure first in CFD. Eliminating r = 1 - f - u and solving:

        f = (beta - k_r - u*(k_u - k_r)) / (k_f - k_r)
    """
    k = load_factors(v)
    u = v["undertray_share"]
    f = ((v["beta_front"] - k["rear_wing"] - u * (k["undertray"] - k["rear_wing"]))
         / (k["front_wing"] - k["rear_wing"]))
    r = 1.0 - f - u
    return f, u, r


def report(v):
    q = 0.5 * RHO * v["v_design"] ** 2          # dynamic pressure [Pa]
    D = q * v["CLA_total"]                      # total downforce [N]
    k = load_factors(v)
    f, u, r = solve_front_wing(v)

    F_fw = f * D
    cla_fw = F_fw / q
    cl_fw = cla_fw / v["S_front_wing"]

    print("=" * 72)
    print("FRONT WING LOAD TARGET")
    print("=" * 72)
    print(f"  ClA_total {v['CLA_total']:.2f}   beta_front {v['beta_front']:.0%}   "
          f"undertray {u:.0%}   v {v['v_design']:.0f} m/s")
    print(f"  Total downforce: {D:.0f} N ({D/9.81:.0f} kgf)\n")

    print("  Front-axle load per newton of device force (the moment arms):")
    print(f"    front wing {k['front_wing']:+5.2f}   "
          f"(1 N of wing -> {k['front_wing']:.2f} N front, "
          f"{1-k['front_wing']:+.2f} N rear: it LIFTS the rear)")
    print(f"    undertray  {k['undertray']:+5.2f}")
    print(f"    rear wing  {k['rear_wing']:+5.2f}   (NEGATIVE: UNLOADS the front)\n")

    print("  Downforce split that satisfies the balance:")
    print(f"    front wing {f:6.1%}   ({F_fw:5.0f} N)")
    print(f"    undertray  {u:6.1%}   ({u*D:5.0f} N)   [input]")
    print(f"    rear wing  {r:6.1%}   ({r*D:5.0f} N)")
    if r < 0:
        print("    !! NEGATIVE rear wing: this balance is unreachable with these")
        print("       CoP positions. The front wing + undertray already over-front")
        print("       the car. Move the undertray CoP rearward or cut beta_front.")
    print()
    print(f"  => FRONT WING must produce {F_fw:.0f} N at {v['v_design']:.0f} m/s")
    print(f"     ClA_fw = {cla_fw:.2f} m^2")
    print(f"     CL_fw  = {cl_fw:.2f}   (over S = {v['S_front_wing']:.2f} m^2, 3-D wing CL)")
    print()

    # --- spanwise split: center (main plane only) vs outboard (3 elements) ---
    c = v["center_span_frac"]
    print("-" * 72)
    print("SPANWISE SPLIT - what the center loading costs the outboard")
    print("-" * 72)
    print(f"  center = {c:.0%} of span (main plane only, feeds the undertray)")
    print(f"  outboard = {1-c:.0%} of span (3 elements)")
    print(f"  The area-averaged wing CL must equal {cl_fw:.2f} either way, so:\n")
    print(f"    {'center CL':>10s}  {'-> outboard CL required':>24s}")
    for cl_c in (0.7, 0.9, 1.1, 1.3, 1.5):
        cl_o = (cl_fw - c * cl_c) / (1 - c)
        flag = "  <-- outboard likely INFEASIBLE (>2.4 in GE)" if cl_o > 2.4 else ""
        print(f"    {cl_c:10.2f}  {cl_o:24.2f}{flag}")
    print()
    print("  READ THIS BACKWARDS. The outboard's MAXIMUM achievable CL sets the")
    print("  center's MINIMUM loading. Unloading the center to protect the")
    print("  undertray inlet is only affordable if the outboard can absorb the")
    print("  deficit. That - not the inlet wake - is what bounds the center from")
    print("  below.")
    print()


def sensitivity(v):
    q = 0.5 * RHO * v["v_design"] ** 2
    D = q * v["CLA_total"]
    print("-" * 72)
    print("SENSITIVITY - which unknown actually moves the answer")
    print("-" * 72)
    print(f"  {'beta':>6s} {'undertray':>10s} {'fw share':>9s} {'F_fw':>8s} "
          f"{'ClA_fw':>8s} {'CL_fw':>7s}")
    for beta in (0.35, 0.40, 0.45):
        for u in (0.35, 0.45, 0.55):
            vv = dict(v, beta_front=beta, undertray_share=u)
            f, _, r = solve_front_wing(vv)
            F = f * D
            print(f"  {beta:6.2f} {u:10.2f} {f:8.1%} {F:7.0f} N "
                  f"{F/q:8.2f} {F/q/v['S_front_wing']:7.2f}"
                  + ("   (rear wing < 0!)" if r < 0 else ""))
    print()
    print("  The undertray share swings the front wing target about as hard as the")
    print("  whole 0.35-0.45 balance range does. It is an ASSUMPTION. Pin it down")
    print("  before trusting any single number above.")
    print()
    print("  Everything scales linearly with ClA_total: if you land at 3.0 instead")
    print("  of 4.0, every front wing force and CL above drops by 25%.")
    print()


def caveats():
    print("=" * 72)
    print("WHAT THIS DOES NOT KNOW (say these before a judge says them to you)")
    print("=" * 72)
    print("""
  1. CL_fw above is a 3-D WING CL. The airfoil screen's cl_target is a 2-D
     SECTION CL. In ground effect and inside a multi-element stack these are
     NOT the same number. Mapping one to the other needs CFD. So this script
     narrows the box; it does not hand the screen its target.

  2. Ground effect is not modeled. A front wing at h/c ~ 0.1-0.3 makes far more
     downforce than free-air numbers suggest, and the CoP moves with ride
     height. Both the load factors and CL_fw are free-air statements.

  3. Device interaction is not modeled. The front wing's wake lands on the
     undertray inlet; the undertray's share is therefore partly a FUNCTION of
     the front wing, not an independent input. This is circular, and the only
     way out is CFD. Bootstrap, then iterate.

  4. The moment balance is static. It ignores CoP migration under pitch/heave,
     which for an undertray car is a real driveability effect.

  5. beta_front is a target, not a physical constraint. Whether it is the RIGHT
     target needs a 2-axle (ideally two-track) lap sim with load-sensitive
     tires. A point-mass sim structurally cannot produce an aero balance.
""")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cla", type=float, default=None, help="total ClA target [m^2]")
    ap.add_argument("--beta", type=float, default=None, help="front AXLE downforce share")
    ap.add_argument("--undertray", type=float, default=None, help="undertray share of total")
    ap.add_argument("--v", type=float, default=None, help="design speed [m/s]")
    ap.add_argument("--no-sensitivity", action="store_true")
    args = ap.parse_args()

    v = dict(VEHICLE)
    if args.cla is not None:
        v["CLA_total"] = args.cla
    if args.beta is not None:
        v["beta_front"] = args.beta
    if args.undertray is not None:
        v["undertray_share"] = args.undertray
    if args.v is not None:
        v["v_design"] = args.v

    report(v)
    if not args.no_sensitivity:
        sensitivity(v)
    caveats()


if __name__ == "__main__":
    main()
