#!/usr/bin/env python3
"""
2D multi-element front-wing section mesher (gmsh) - CFD stage, step 1.

Takes shortlist sections from airfoil_screen.py, rigs them into the outboard
three-element stack (main + flap1 + flap2), and produces a solver-ready 2D
RANS mesh in ground effect:

  - quad boundary layers on every element, first-cell height from a y+ target
    (default y+=1 at V_hi: transition-SST-ready, since the whole screen's
    logic rests on laminar/transition physics that fully-turbulent wall
    functions would throw away)
  - moving-ground-plane domain at parametric ride height
  - slot / wake / ground refinement zones
  - named boundary zones (inlet/outlet/top/ground/main/flap1/flap2) so BC
    assignment is scriptable in any solver

WHY GMSH AND MULTI-FORMAT EXPORT
--------------------------------
Compute resources for this year are not locked in, so the mesh stage must not
be. gmsh is free, fully scriptable, and exports SU2 (.su2), Fluent-importable
CGNS/UNV, legacy VTK, and gmsh .msh v2.2 (OpenFOAM's gmshToFoam). One mesh
script therefore serves Fluent, SU2, and OpenFOAM; the solver decision can
wait for the license situation without blocking mesh development. Every export
is attempted independently and the report records which succeeded.

GEOMETRY SOURCE
---------------
Prefers kulfan_params.json (the EXACT geometry the screen ranked - see
common.py rule 2). Falls back to shortlist_dat/*.dat with a warning, because
the .dat files are a lossy XFLR5 export: a coordinate round-trip yields a
slightly different section than the one that was ranked.

RIGGING CONVENTIONS (all parametric - these are the optimizer's variables)
--------------------------------------------------------------------------
Rigging is specified in the standard lift-up airfoil frame (the frame the
screen ranked in), then the whole assembly is mirrored into car orientation
(suction side toward the ground). Signs: POSITIVE = MORE LOAD.

  alpha_main   main-element incidence [deg]
  delta1       flap1 deflection relative to the MAIN chord line [deg]
  delta2       flap2 deflection relative to the FLAP1 chord line [deg]
  gap1/gap2    slot offset, fraction of PARENT chord, normal to parent chord
  overlap1/2   flap LE upstream of parent TE, fraction of PARENT chord
  ride_height  lowest point of the assembly above ground [m]

gap is the TRUE minimum slot gap (min surface-to-surface distance between
the flap and its parent), not a placement offset. The flap LE starts at
parent_TE - overlap*c_hat + gap*n_hat (c_hat = parent chord direction, n_hat
= parent pressure-side normal), is rotated about its own LE, then translated
along n_hat until the measured min gap equals the request (bisection; an
analytic placement alone closed the slot to <1 mm once the flap rotated).
Overlap is preserved by that translation (it is chordwise by construction).

Defaults are conventional FSAE starting values (main +3 deg, flaps +20/+25,
~2.5%c slots, ~3%c overlap, 50 mm ride height), meant as the CENTER of the
later optimization, not as an answer.

SCOPE: outboard stack only for now. The center main-only section is deferred
until the 2D nose+undertray interaction study defines its environment.

Usage (after a screen run):
  python mesh_section.py                          # rank-1 sections, defaults
  python mesh_section.py --airfoils MAIN F1 F2    # explicit sections
  python mesh_section.py --ride-height 0.035 --delta1 24 --gap1 0.02
  python mesh_section.py --no-ground              # freestream domain (debug)
Outputs (in --out, default mesh_out/):
  <tag>.su2 / .cgns / .unv / .msh / .vtk   whichever formats succeed
  <tag>_preview.png                        geometry + mesh close-ups
  <tag>_report.json                        counts, quality, y+ basis, gaps
"""

import argparse
import json
import math
import os
import sys

import numpy as np

from common import CONFIG, NU, PROFILES, V_DESIGN, load_kulfan, load_run_config

RHO = 1.225                       # sea-level density [kg/m^3]

# element -> (profile key, points per side for repaneling)
ELEMENTS = [("main", "main_outboard", 160),
            ("flap1", "flap1", 130),
            ("flap2", "flap2", 110)]

TE_COLLAPSE = 3.0e-4              # collapse blunter TEs than this [m]? no:
                                  # collapse SHARPER (base < this) to a point

RIG_DEFAULTS = dict(alpha_main=3.0, delta1=20.0, delta2=25.0,
                    gap1=0.025, overlap1=0.030,
                    gap2=0.025, overlap2=0.030,
                    ride_height=0.050)


# ---------------------------------------------------------------------------
# Geometry loading
# ---------------------------------------------------------------------------
def load_section(name, results_dir, profile_key=None):
    """(N,2) unit-chord coordinates, TE->upper->LE->lower->TE (Selig order).

    kulfan_params.json first (exact ranked geometry); .dat fallback + warning.

    THE TRAILING EDGE IS BLUNTED HERE, and it must be.
    ---------------------------------------------------
    The screen does not rank the airfoil as drawn. These wings are hot-wire cut
    from foam, and a knife-edge TE cannot be cut - so airfoil_screen.py sets
    TE_thickness = te_min_mm / chord for EACH ELEMENT before running any polar,
    and ranks the blunted section. (A 1.5 mm edge is 0.55% of the 275 mm main
    chord but 1.36% of the 110 mm flap2 chord, so the blunt is per-element and
    cannot be baked into a single per-name geometry file - which is exactly why
    kulfan_params.json still stores the AS-DRAWN section.)

    If we meshed the as-drawn geometry, the CFD would solve a knife-edge section
    that was never screened, never buildable, and whose aft loading differs from
    the one that won. So we re-apply the identical blunt here, using te_min_mm
    from the run's own run_config.json.
    """
    import aerosandbox as asb

    te_applied = None
    if profile_key is not None:
        chord_mm = PROFILES[profile_key]["chord"] * 1000.0
        te_needed = CONFIG["te_min_mm"] / chord_mm          # fraction of chord

    try:
        kulfan = load_kulfan(results_dir)
        if name in kulfan:
            kp = dict(kulfan[name])
            if profile_key is not None:
                te_as_drawn = float(kp["TE_thickness"])
                te_applied = max(te_as_drawn, te_needed)
                kp["TE_thickness"] = te_applied
                print(f"  {name} [{profile_key}]: TE {te_as_drawn*chord_mm:.2f} mm "
                      f"-> {te_applied*chord_mm:.2f} mm (blunted, as screened)")
            return asb.KulfanAirfoil(name=name, **kp)
    except FileNotFoundError:
        pass

    # Fallback: the .dat files are now written per element and ALREADY blunted
    # (shortlist_dat/<name>__<element>.dat), so prefer that one.
    datd = os.path.join(results_dir, "shortlist_dat")
    dat = os.path.join(datd, f"{name}__{profile_key}.dat") if profile_key else None
    if dat is None or not os.path.exists(dat):
        dat = os.path.join(datd, f"{name}.dat")          # legacy, as-drawn
    if not os.path.exists(dat):
        sys.exit(f"'{name}' not in kulfan_params.json and no .dat in {datd}.")
    print(f"  ! {name}: using .dat fallback (kulfan_params.json missing) - "
          f"coordinate round-trip, section differs slightly from the ranked one.")
    return asb.Airfoil(name=name, coordinates=np.loadtxt(dat, skiprows=1))


def rank1_sections(results_dir):
    """Default test sections: rank-1 of each element sheet in results.xlsx."""
    import pandas as pd
    xlsx = os.path.join(results_dir, "results.xlsx")
    picks = []
    for _, profile, _ in ELEMENTS:
        top = pd.read_excel(xlsx, sheet_name=f"{profile}_top")
        picks.append(str(top["name"].iloc[0]))
    return picks


def unit_coords(af, n_per_side):
    """Repaneled closed coordinates, tiny TE base collapsed to a point."""
    c = np.array(af.repanel(n_points_per_side=n_per_side).coordinates, float)
    # Selig order: first & last point are the TE (upper/lower). If the base
    # is tiny it would force micro-elements at the TE; collapse to midpoint.
    if np.linalg.norm(c[0] - c[-1]) < TE_COLLAPSE:  # unit chord here; scaled
        mid = 0.5 * (c[0] + c[-1])                  # later, but chords ~O(0.1
        c[0] = c[-1] = mid                          # -0.3 m) so same order
        c = c[:-1]                                  # drop duplicate: loop
        closed = True                               # closes last->first
    else:
        closed = False                              # keep blunt TE base edge
    return c, closed


# ---------------------------------------------------------------------------
# Rigging
# ---------------------------------------------------------------------------
def rot(pts, deg):
    a = math.radians(deg)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return pts @ R.T


def rig(sections, rig_p):
    """Place/scale/rotate the three elements; mirror to car frame.

    sections: {"main": (coords, closed), ...} unit chord, lift-up frame.
    Returns {"main": (N,2) car-frame coords, ...}, rigging_info dict.

    The gap parameter is enforced as the MEASURED min slot gap: the flap is
    first placed analytically, then translated along the parent's pressure-
    side normal until min surface-to-surface distance hits the request.
    (The analytic placement alone under-delivers badly - the rotated flap's
    suction surface swings up toward the parent TE and closes the slot.)
    """
    chords = {e: PROFILES[p]["chord"] for e, p, _ in ELEMENTS}
    theta = {"main": rig_p["alpha_main"],
             "flap1": rig_p["alpha_main"] + rig_p["delta1"]}
    theta["flap2"] = theta["flap1"] + rig_p["delta2"]

    placed, le = {}, {}
    le["main"] = np.zeros(2)
    placed["main"] = rot(sections["main"][0] * chords["main"], -theta["main"])
    for elem, parent, gap_k, ov_k in [("flap1", "main", "gap1", "overlap1"),
                                      ("flap2", "flap1", "gap2", "overlap2")]:
        tp = math.radians(theta[parent])
        c_hat = np.array([math.cos(tp), -math.sin(tp)])   # LE->TE, incidence
        n_hat = np.array([-math.sin(tp), -math.cos(tp)])  # pressure side
        te_parent = le[parent] + chords[parent] * c_hat
        cp = chords[parent]
        target = rig_p[gap_k] * cp
        le0 = te_parent - rig_p[ov_k] * cp * c_hat + target * n_hat
        flap0 = rot(sections[elem][0] * chords[elem], -theta[elem]) + le0

        def f(t):
            return min_gap(placed[parent], flap0 + t * n_hat) - target

        # Bracket the NEAR-SIDE root only. f is not monotonic: pushing the
        # flap up through the parent also "opens a gap" on the far side, so
        # naive bracket doubling can converge on a rig with the flap parked
        # 200 mm above the wing. +t moves away from the parent along n_hat;
        # gap grows monotonically for t > contact.
        if f(0.0) < 0:                              # too tight (usual case)
            t_lo, t_hi = 0.0, 0.005
            while f(t_hi) < 0 and t_hi < 0.10:
                t_hi += 0.005
        else:                                       # too open: walk down in
            t_hi = 0.0                              # small steps, stop at the
            t_lo = -0.002                           # first f<0 (near contact)
            while f(t_lo) > 0 and t_lo > -0.04:
                t_hi = t_lo
                t_lo -= 0.002
        for _ in range(30):                         # bisection to ~ um
            t_mid = 0.5 * (t_lo + t_hi)
            if f(t_mid) < 0:
                t_lo = t_mid
            else:
                t_hi = t_mid
        t = 0.5 * (t_lo + t_hi)
        le[elem] = le0 + t * n_hat
        placed[elem] = flap0 + t * n_hat

    # Mirror to car frame (downforce), set ride height, main LE at x=0.
    y_min = min(-p[:, 1].max() for p in placed.values())
    for elem, p in placed.items():
        q = p.copy()
        q[:, 1] = -q[:, 1] + (rig_p["ride_height"] - y_min)
        q[:, 0] -= le["main"][0]
        placed[elem] = q

    info = dict(theta_deg=theta, chords=chords,
                le_car={e: [float(placed[e][:, 0].min()),
                            float(placed[e][np.argmin(placed[e][:, 0]), 1])]
                        for e in placed})
    return placed, info


def min_gap(pa, pb):
    """Min distance between two closed polylines (vertex-to-segment, both ways)."""
    def v2s(pts, poly):
        a, b = poly, np.roll(poly, -1, axis=0)
        ab = b - a                                    # (M,2)
        d2 = np.einsum("ij,ij->i", ab, ab)
        best = np.inf
        for p in pts:
            t = np.clip(np.einsum("ij,j->i", ab, p) - np.einsum("ij,ij->i", ab, a), 0, d2)
            t = np.divide(t, d2, out=np.zeros_like(d2), where=d2 > 0)
            proj = a + ab * t[:, None]
            best = min(best, float(np.min(np.linalg.norm(proj - p, axis=1))))
        return best
    return min(v2s(pa, pb), v2s(pb, pa))


# ---------------------------------------------------------------------------
# Boundary-layer sizing
# ---------------------------------------------------------------------------
def first_layer(chord, v, yplus):
    """First-cell height for a y+ target (turbulent flat-plate estimate)."""
    re = v * chord / NU
    cf = 0.058 * re ** -0.2
    u_tau = math.sqrt(0.5 * cf) * v
    return yplus * NU / u_tau, re


def bl_spec(v_ref, yplus, slot_gap, ratio=1.2):
    """One global BL spec (first height h1, layers n, total thickness T).

    h1 from the LARGEST-chord element (largest Re -> thinnest first cell:
    conservative for all three). Total thickness targets the main element's
    turbulent delta99 but is capped at 40% of the measured minimum slot gap
    so opposing BL extrusions cannot collide inside the slot; the region
    between BL cap and delta99 is covered by the isotropic slot refinement.
    """
    c_main = PROFILES["main_outboard"]["chord"]
    h1, re = first_layer(c_main, v_ref, yplus)
    delta99 = 0.37 * c_main / re ** 0.2
    t_total = min(delta99, 0.40 * slot_gap)

    # FLOOR, not ceil. n = ceil() then recomputing the geometric sum OVERSHOOTS
    # t_total - which defeats the entire purpose of the slot-gap cap. Measured:
    # cap 1.70 mm, actual 1.83 mm, so two opposing BL extrusions ate 85% of the
    # 4.3 mm flap1->flap2 slot instead of the intended 80%, leaving 0.64 mm of
    # core. The cap exists to stop the boundary layers colliding inside the
    # slot; a cap that is exceeded by construction is not a cap.
    n = max(1, math.floor(math.log(1 + t_total * (ratio - 1) / h1)
                          / math.log(ratio)))
    t_actual = h1 * (ratio ** n - 1) / (ratio - 1)
    assert t_actual <= t_total * 1.001, "BL still exceeds its slot-gap cap"

    if n < 8:
        print(f"  ! BL has only {n} layers: the {slot_gap*1000:.1f} mm slot gap "
              f"caps total BL thickness at {t_total*1000:.2f} mm, which y+={yplus} "
              f"at ratio {ratio} cannot fill in 8+ layers.")
        print(f"    The near-wall BL is UNDER-RESOLVED. Open the slot "
              f"(--gap1/--gap2), raise --yplus, or accept wall functions - but "
              f"do not run transition-SST on this mesh and trust it.")

    return dict(h1=h1, ratio=ratio, n_layers=n, thickness=t_actual,
                delta99_main=delta99, Re_main=re, cap=t_total,
                capped_by="slot_gap" if 0.40 * slot_gap < delta99 else "delta99")


# ---------------------------------------------------------------------------
# gmsh
# ---------------------------------------------------------------------------
def build_mesh(placed, sections, spec, dom, out_base, formats):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    geo = gmsh.model.geo
    loops, all_curves, fan_pts = {}, [], []

    for elem, pts in placed.items():
        closed = sections[elem][1]
        n = len(pts)
        seg = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
        size = 0.5 * (seg + np.roll(seg, 1))          # local repanel spacing
        tags = [geo.addPoint(x, y, 0, s) for (x, y), s in zip(pts, size)]
        lines = [geo.addLine(tags[i], tags[i + 1]) for i in range(n - 1)]
        lines.append(geo.addLine(tags[-1], tags[0]))  # close (TE pt or base)
        loops[elem] = geo.addCurveLoop(lines)
        all_curves += lines
        fan_pts += [tags[0]] if closed else [tags[0], tags[-1]]  # TE fan(s)

    x0, x1, y0, y1 = dom["x0"], dom["x1"], dom["y0"], dom["y1"]
    s_far = dom["far_size"]
    c = [geo.addPoint(x0, y0, 0, s_far), geo.addPoint(x1, y0, 0, s_far),
         geo.addPoint(x1, y1, 0, s_far), geo.addPoint(x0, y1, 0, s_far)]
    ground = geo.addLine(c[0], c[1])
    outlet = geo.addLine(c[1], c[2])
    top = geo.addLine(c[2], c[3])
    inlet = geo.addLine(c[3], c[0])
    outer = geo.addCurveLoop([ground, outlet, top, inlet])
    surf = geo.addPlaneSurface([outer] + list(loops.values()))
    geo.synchronize()

    pg = gmsh.model.addPhysicalGroup
    pg(1, [inlet], name="inlet")
    pg(1, [outlet], name="outlet")
    pg(1, [top], name="top")
    pg(1, [ground], name="ground")
    # per-element curve slices (all_curves was filled in element order)
    i = 0
    for elem, pts in placed.items():
        n_lines = len(pts)
        pg(1, all_curves[i:i + n_lines], name=elem)
        i += n_lines
    pg(2, [surf], name="fluid")

    fld = gmsh.model.mesh.field
    # near-surface isotropic refinement (also resolves the slots)
    f_dist = fld.add("Distance"); fld.setNumbers(f_dist, "CurvesList", all_curves)
    fld.setNumber(f_dist, "Sampling", 200)
    f_near = fld.add("Threshold")
    fld.setNumber(f_near, "InField", f_dist)
    fld.setNumber(f_near, "SizeMin", dom["near_size"])
    fld.setNumber(f_near, "SizeMax", s_far)
    fld.setNumber(f_near, "DistMin", spec["thickness"])
    fld.setNumber(f_near, "DistMax", dom["near_dist"])
    # wake box
    f_wake = fld.add("Box")
    fld.setNumber(f_wake, "VIn", dom["wake_size"]); fld.setNumber(f_wake, "VOut", s_far)
    fld.setNumber(f_wake, "XMin", dom["wake_x0"]); fld.setNumber(f_wake, "XMax", dom["wake_x1"])
    fld.setNumber(f_wake, "YMin", y0); fld.setNumber(f_wake, "YMax", dom["wake_y1"])
    fld.setNumber(f_wake, "Thickness", 0.3)
    # ground strip under/behind the wing
    f_gnd = fld.add("Box")
    fld.setNumber(f_gnd, "VIn", dom["gnd_size"]); fld.setNumber(f_gnd, "VOut", s_far)
    fld.setNumber(f_gnd, "XMin", dom["gnd_x0"]); fld.setNumber(f_gnd, "XMax", dom["gnd_x1"])
    fld.setNumber(f_gnd, "YMin", y0); fld.setNumber(f_gnd, "YMax", y0 + dom["gnd_h"])
    fld.setNumber(f_gnd, "Thickness", 0.05)
    f_min = fld.add("Min")
    fld.setNumbers(f_min, "FieldsList", [f_near, f_wake, f_gnd])
    fld.setAsBackgroundMesh(f_min)

    # quad boundary layers on the airfoils
    f_bl = fld.add("BoundaryLayer")
    fld.setNumbers(f_bl, "CurvesList", all_curves)
    fld.setNumber(f_bl, "Size", spec["h1"])
    fld.setNumber(f_bl, "Ratio", spec["ratio"])
    fld.setNumber(f_bl, "Thickness", spec["thickness"])
    fld.setNumber(f_bl, "Quads", 1)
    if fan_pts:
        fld.setNumbers(f_bl, "FanPointsList", fan_pts)
    fld.setAsBoundaryLayer(f_bl)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)        # frontal-Delaunay
    gmsh.model.mesh.generate(2)

    # ---- stats ----
    types, tags_, _ = gmsh.model.mesh.getElements(2)
    counts = {}
    for t, tg in zip(types, tags_):
        name = gmsh.model.mesh.getElementProperties(t)[0]
        counts[name] = len(tg)
    n_cells = sum(counts.values())
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    # Quality PER TYPE: BL quads are intentionally anisotropic (AR ~ 100+),
    # which SICN scores near zero by construction. Judge mesh health on the
    # triangles; judge the quads only against their own population.
    qual = {}
    for t, tg in zip(types, tags_):
        name = gmsh.model.mesh.getElementProperties(t)[0]
        q = gmsh.model.mesh.getElementQualities(tg, "minSICN")
        qual[name] = dict(min=float(np.min(q)), mean=float(np.mean(q)),
                          frac_below_0p1=float(np.mean(q < 0.1)))
    stats = dict(cells_by_type=counts, n_cells=int(n_cells),
                 n_nodes=int(n_nodes), quality_minSICN=qual)

    written = {}
    for ext in formats:
        try:
            if ext == "msh":
                gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(f"{out_base}.{ext}")
            written[ext] = True
        except Exception as e:
            written[ext] = f"FAILED: {e}"
    stats["exports"] = written

    # node/element arrays for the preview, before gmsh goes away
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    xy = np.array(ncoords, float).reshape(-1, 3)[:, :2]
    remap = np.zeros(int(ntags.max()) + 1, dtype=np.int64)
    remap[np.asarray(ntags, dtype=np.int64)] = np.arange(len(ntags))
    elem_nodes = []
    types, tags_, nodess = gmsh.model.mesh.getElements(2)
    for t, nod in zip(types, nodess):
        nn = gmsh.model.mesh.getElementProperties(t)[3]
        elem_nodes.append(remap[np.asarray(nod, dtype=np.int64)].reshape(-1, nn))
    gmsh.finalize()
    return stats, xy, elem_nodes


# ---------------------------------------------------------------------------
# Preview
# ---------------------------------------------------------------------------
def preview(placed, xy, elem_nodes, dom, gaps, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    bx0 = min(p[:, 0].min() for p in placed.values())
    bx1 = max(p[:, 0].max() for p in placed.values())
    by1 = max(p[:, 1].max() for p in placed.values())
    pad = 0.06

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))

    ax[0].add_patch(plt.Rectangle((dom["x0"], dom["y0"]),
                                  dom["x1"] - dom["x0"], dom["y1"] - dom["y0"],
                                  fill=False, lw=1))
    for e, p in placed.items():
        ax[0].fill(p[:, 0], p[:, 1], lw=0.5, label=e)
    ax[0].axhline(dom["y0"], color="k", lw=2)
    ax[0].set_title("domain (ground plane at y=0)")
    ax[0].legend(fontsize=8); ax[0].set_aspect("equal")

    def mesh_panel(a, x0, x1, y0, y1, title):
        for nodes in elem_nodes:
            cx = xy[nodes[:, 0], 0]
            cy = xy[nodes[:, 0], 1]
            m = (cx > x0 - 0.02) & (cx < x1 + 0.02) & (cy > y0 - 0.02) & (cy < y1 + 0.02)
            if m.any():
                a.add_collection(PolyCollection(
                    xy[nodes[m]], facecolors="none", edgecolors="0.4", lw=0.15))
        for p in placed.values():
            a.plot(p[:, 0], p[:, 1], "k-", lw=0.7)
        a.axhline(dom["y0"], color="k", lw=1.5)
        a.set_xlim(x0, x1); a.set_ylim(y0, y1)
        a.set_aspect("equal"); a.set_title(title, fontsize=9)

    mesh_panel(ax[1], bx0 - pad, bx1 + pad, dom["y0"], by1 + pad,
               "mesh: full stack + ground")
    f1 = placed["flap1"]
    lex, ley = f1[np.argmin(f1[:, 0])]
    z = 0.030
    mesh_panel(ax[2], lex - z, lex + z, ley - z, ley + z,
               f"slot 1 zoom (measured gaps: "
               f"{1000*gaps['main_flap1']:.1f} / {1000*gaps['flap1_flap2']:.1f} mm)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="3-element front-wing 2D mesher")
    ap.add_argument("--results", default="screen_results")
    # Default None -> <results>/mesh_out. The mesh belongs WITH the run that
    # produced it: same directory as results.xlsx, run_config.json and
    # kulfan_params.json. Writing it to a cwd-relative 'mesh_out' orphaned it -
    # no link back to the config it was built from, and a re-screen silently
    # left a stale mesh sitting next to fresh results.
    ap.add_argument("--out", default=None,
                    help="mesh output dir (default: <results>/mesh_out)")
    ap.add_argument("--airfoils", nargs=3, metavar=("MAIN", "FLAP1", "FLAP2"),
                    default=None, help="section names; default: rank-1 per element")
    for k, v in RIG_DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_','-')}", type=float, default=v)
    ap.add_argument("--yplus", type=float, default=1.0)
    # Imported from common.V_DESIGN, NOT hardcoded. This was 27.0 - which
    # happened to equal V_DESIGN[-1], so it looked right and was a trap: change
    # the speed band in common.py and the mesh would keep sizing its boundary
    # layer for the old speed, silently. Same bug class as the MODEL_SIZE
    # bind-by-value: a constant duplicated instead of imported.
    ap.add_argument("--v-ref", type=float, default=V_DESIGN[-1],
                    help=f"sizing speed [m/s]; default V_hi = {V_DESIGN[-1]} "
                         f"(highest Re -> thinnest first cell -> conservative)")
    ap.add_argument("--near-size", type=float, default=2.0e-3,
                    help="isotropic size at the surface/slots [m]")
    ap.add_argument("--no-ground", action="store_true",
                    help="freestream domain (mirror-image debugging only)")
    ap.add_argument("--formats", default="su2,msh,unv,cgns,vtk")
    args = ap.parse_args()

    # Load the config the SCREEN actually ran with, before touching geometry.
    # te_min_mm lives here, and the TE blunt must match the screened section
    # exactly - otherwise we solve a shape that was never ranked.
    load_run_config(args.results)

    names = args.airfoils or rank1_sections(args.results)
    rig_p = {k: getattr(args, k) for k in RIG_DEFAULTS}
    print("Sections:", dict(zip([e for e, _, _ in ELEMENTS], names)))
    print("Rigging:", rig_p)
    print(f"TE blunt: {CONFIG['te_min_mm']} mm (from run_config.json) - "
          f"applied per element, as screened")

    sections = {}
    for (elem, profile, npts), name in zip(ELEMENTS, names):
        af = load_section(name, args.results, profile_key=profile)
        sections[elem] = unit_coords(af, npts)

    placed, info = rig(sections, rig_p)
    gaps = dict(main_flap1=min_gap(placed["main"], placed["flap1"]),
                flap1_flap2=min_gap(placed["flap1"], placed["flap2"]))
    ground_clear = min(p[:, 1].min() for p in placed.values())
    print(f"Measured slot gaps: main-flap1 {1000*gaps['main_flap1']:.1f} mm, "
          f"flap1-flap2 {1000*gaps['flap1_flap2']:.1f} mm; "
          f"ground clearance {1000*ground_clear:.1f} mm")

    spec = bl_spec(args.v_ref, args.yplus, min(gaps.values()))
    print(f"BL: first cell {spec['h1']*1e6:.1f} um, {spec['n_layers']} layers "
          f"@ ratio {spec['ratio']}, total {spec['thickness']*1000:.2f} mm "
          f"(delta99_main ~ {spec['delta99_main']*1000:.1f} mm)")

    c = PROFILES["main_outboard"]["chord"]
    bx1 = max(p[:, 0].max() for p in placed.values())
    by1 = max(p[:, 1].max() for p in placed.values())
    y_bottom = 0.0 if not args.no_ground else -10 * c
    dom = dict(x0=-8 * c, x1=16 * c, y0=y_bottom, y1=10 * c,
               far_size=0.15, near_size=args.near_size,
               near_dist=0.10,
               wake_size=6.0e-3, wake_x0=-0.05, wake_x1=bx1 + 4 * c,
               wake_y1=by1 + 0.05,
               gnd_size=2.5e-3, gnd_x0=-0.3, gnd_x1=bx1 + 3 * c, gnd_h=0.02)

    outdir = args.out or os.path.join(args.results, "mesh_out")
    os.makedirs(outdir, exist_ok=True)
    tag = "__".join(names) + "__stack"
    out_base = os.path.join(outdir, tag)
    print(f"Mesh output -> {outdir}/")
    stats, xy, elem_nodes = build_mesh(placed, sections, spec, dom, out_base,
                                       [f.strip() for f in args.formats.split(",")])

    report = dict(sections=dict(zip([e for e, _, _ in ELEMENTS], names)),
                  rigging=rig_p, rigging_resolved=info,
                  measured_slot_gaps_m=gaps,
                  ground_clearance_m=float(ground_clear),
                  bl=spec, yplus_target=args.yplus, v_ref=args.v_ref,
                  domain=dom, mesh=stats)
    with open(f"{out_base}_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    preview(placed, xy, elem_nodes, dom, gaps, f"{out_base}_preview.png")

    print(f"\nMesh: {stats['n_cells']} cells / {stats['n_nodes']} nodes "
          f"({stats['cells_by_type']})")
    for tname, qq in stats["quality_minSICN"].items():
        note = "  (BL quads: low SICN = intended anisotropy)" if "Quad" in tname else ""
        print(f"Quality minSICN [{tname}]: min {qq['min']:.3f} "
              f"mean {qq['mean']:.3f} (<0.1: {100*qq['frac_below_0p1']:.1f}%){note}")
    print(f"Exports: {stats['exports']}")
    print(f"Wrote {out_base}_report.json / _preview.png")
    if stats["n_cells"] + stats["n_nodes"] > 512_000:
        print("  ! cells+nodes exceeds the 512k Ansys Student limit - "
              "coarsen --near-size / wake sizes if that is your license")


if __name__ == "__main__":
    main()  # noqa
