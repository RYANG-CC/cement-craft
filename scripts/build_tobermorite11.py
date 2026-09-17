#!/usr/bin/env python3
"""
Build bulk tobermorite 11 A - normal and anomalous forms - from
Merlino, Bonaccorsi & Armbruster (2001), Eur. J. Mineral. 13, 577-590, "The
real structure of tobermorite 11 A: normal and anomalous forms, OD character
and polytypic modifications".

Sources (both CC0 on the Crystallography Open Database, both MDO2 / B11m):
  anomalous  COD 9005498 (Wessels mine, South Africa). Ca2 Si3 H5 O11,
             Ca/Si 0.667. No partial occupancy: the anomalous form is missing
             the Ca2 interlayer site entirely, not disordered on it.
  normal     COD 9005499 (Basenov, Urals, Russia). Ca2.25 Si3 H7 O11,
             Ca/Si 0.75. Ca2 (occ 0.25) and the Wat1/Wat3 waters (occ 0.5)
             are genuinely disordered - resolved explicitly below (see
             `_resolve_disorder`), one ordering, not interpolated or ignored.

Both refinements share the same 4 listed symmetry operators (a mirror at
z=0 combined with the B-centring translation), so one parser and one
disorder resolver serve both files; the anomalous file simply has nothing
partial to resolve.
"""
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from cement_charge import ChargeModel, formal_charge          # noqa: E402
from cement_core import (                                     # noqa: E402
    Structure, add_water_hydrogens, cell_from_parameters, close_contacts,
    detect_groups, transform_cell, write_xyz, _align_lone_pair,
)
from cement_surface import FRAMEWORK_METALS, _add_h           # noqa: E402
from cement_phases import O_H                                 # noqa: E402

CIF_DIR = ROOT / "assets" / "cif_sources"
OUT_DIR = ROOT / "assets" / "structures"

# Shared by both CIFs: identity, B-centring translation, the mirror at z=0,
# and the mirror combined with centring.
SYMOPS = [
    lambda x, y, z: (x, y, z),
    lambda x, y, z: (0.5 + x, y, 0.5 + z),
    lambda x, y, z: (x, y, -z),
    lambda x, y, z: (0.5 + x, y, 0.5 - z),
]

TETRAHEDRAL_ANGLE = 109.47   # sp3 initialisation angle; CELL_OPT relaxes it


def _parse_atom_site_loop(text):
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip() != "loop_":
            continue
        j = i + 1
        headers = []
        while j < len(lines) and lines[j].strip().startswith("_"):
            headers.append(lines[j].strip())
            j += 1
        if "_atom_site_label" not in headers:
            continue
        rows = []
        while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(("_", "loop_", "#")):
            rows.append(lines[j].split())
            j += 1
        return headers, rows
    raise ValueError("no _atom_site_label loop found")


def read_tobermorite_cif(path):
    """Parse one of these two specific CIFs: cell, symmetry-expanded sites.

    Deliberately not routed through the general-purpose `read_cif` - that
    collapses every site to its bare element symbol, which is exactly the
    distinction (Ca1/Ca3 vs Ca2, Wat2/Wat6 vs Wat1/Wat3) this disorder
    resolution needs to keep.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")

    def num(key):
        import re
        m = re.search(rf"^{key}\s+([-\d.]+)", text, re.MULTILINE)
        return float(m.group(1))

    a, b, c = num("_cell_length_a"), num("_cell_length_b"), num("_cell_length_c")
    al, be, ga = num("_cell_angle_alpha"), num("_cell_angle_beta"), num("_cell_angle_gamma")
    cell = cell_from_parameters(a, b, c, al, be, ga)

    headers, rows = _parse_atom_site_loop(text)
    idx = {h: k for k, h in enumerate(headers)}
    i_lab = idx["_atom_site_label"]
    i_x, i_y, i_z = idx["_atom_site_fract_x"], idx["_atom_site_fract_y"], idx["_atom_site_fract_z"]
    i_occ = idx.get("_atom_site_occupancy")
    i_elem = idx["_atom_site_type_symbol"]
    i_nh = idx.get("_atom_site_attached_hydrogens")

    sites = []
    for row in rows:
        lab = row[i_lab]
        frac = tuple(float(row[k]) for k in (i_x, i_y, i_z))
        occ = float(row[i_occ]) if i_occ is not None else 1.0
        elem = row[i_elem]
        nh = int(float(row[i_nh])) if i_nh is not None else 0
        sites.append(dict(label=lab, frac=frac, occ=occ, elem=elem, nh=nh))
    return cell, sites


def _mic_frac_dist(f1, f2, cell):
    d = np.array(f1) - np.array(f2)
    d = d - np.round(d)
    return np.linalg.norm(d @ cell)


def _expand_symmetry(sites, cell, tol=0.05):
    """Apply the listed symops to every site, deduplicating exact overlaps.

    A site sitting exactly on the z=0 (or z=1/2) mirror plane maps onto
    itself under two of the four operators - that is a lower-multiplicity
    special position (Wat2, O5 here), not a disorder, and is folded down
    to its true count rather than counted twice.
    """
    expanded = {}   # label -> list of dict(frac, elem, occ, nh)
    for s in sites:
        pts = []
        for op in SYMOPS:
            g = tuple(v % 1.0 for v in op(*s["frac"]))
            if not any(_mic_frac_dist(g, p["frac"], cell) < tol for p in pts):
                pts.append(dict(frac=g, elem=s["elem"], occ=s["occ"], nh=s["nh"]))
        expanded[s["label"]] = pts
    return expanded


def _resolve_disorder(expanded, cell, cluster_cutoff=2.5, verbose=True):
    """Resolve genuinely partial-occupancy labels into one explicit ordering.

    Each partial label's symmetry copies separate cleanly into local clusters
    (mirror-related pairs a bond length apart) that sit far from each other
    (related by the lattice centring - tens of angstrom away here, never
    confused with a real split pair). `occupancy * raw copies` gives how many
    of those clusters are genuinely occupied; the lowest-indexed clusters are
    activated ("ordering A", the same deliberate-and-labelled convention used
    for the tobermorite 9 A Ca2 split site), keeping one representative atom
    per activated cluster. This is one ordering among several valid ones.
    """
    resolved = {}
    for label, pts in expanded.items():
        occs = {round(p["occ"], 4) for p in pts}
        if occs == {1.0}:
            resolved[label] = pts
            continue
        occ = pts[0]["occ"]
        n = len(pts)
        target = round(occ * n)

        # union-find clustering by mutual distance
        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                i = parent[i]
            return i

        for i in range(n):
            for j in range(i + 1, n):
                if _mic_frac_dist(pts[i]["frac"], pts[j]["frac"], cell) < cluster_cutoff:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pi] = pj
        clusters: dict[int, list[int]] = {}
        for i in range(n):
            clusters.setdefault(find(i), []).append(i)
        ordered = sorted(clusters.values(), key=min)

        kept = [pts[c[0]] for c in ordered[:target]]
        resolved[label] = kept
        if verbose:
            sizes = [len(c) for c in ordered]
            print(f"  {label}: occ={occ:.3f} x {n} raw copies -> {len(ordered)} "
                  f"cluster(s) of size {sizes}, keeping {target} atom(s) "
                  f"(ordering A: lowest-index cluster(s) first)")
    return resolved


def build_structure(cif_path, verbose=True):
    cell, sites = read_tobermorite_cif(cif_path)
    expanded = _expand_symmetry(sites, cell)
    resolved = _resolve_disorder(expanded, cell, verbose=verbose)

    labels, positions, nh_map = [], [], {}
    for pts in resolved.values():
        for p in pts:
            nh_map[len(labels)] = p["nh"]
            labels.append(p["elem"])
            positions.append(np.array(p["frac"]) @ cell)
    info = {"source": str(cif_path),
            "attached_hydrogens": {i: n for i, n in nh_map.items() if n}}
    return Structure(labels, np.array(positions), cell, info)


# -- hydroxyl placement and charge closing (shared with the 9 A / anomalous
# 11 A pipeline; see build_csh067.py's history for the original derivation) --

def _place_hydroxyl_h(struct, o_index, cov_index, avoid_index=None,
                      hbond_range=(2.5, 3.4), n_trials=200, seed=0):
    """Place a hydroxide H on an oxygen that keeps one covalent bond.

    The Si-O bond is real and fixed, so H goes on the tetrahedral cone
    around it; the remaining rotational freedom is chosen to point away
    from the coordinating cation and toward a hydrogen-bond acceptor. See
    references/charge.md for why terminal, least-coordinated oxygens are the
    right protonation candidates.
    """
    rng = np.random.default_rng(seed)
    elems = struct.elements
    origin = struct.positions[o_index]
    v_cov = struct.mic(struct.positions[cov_index] - origin)[0]
    v_cov /= np.linalg.norm(v_cov)
    rot = _align_lone_pair(v_cov)
    theta = math.radians(TETRAHEDRAL_ANGLE)

    d = struct.distances_from(o_index)
    acc = [j for j in np.where((d > hbond_range[0]) & (d < hbond_range[1]))[0]
           if elems[j] == "O" and j != o_index]
    acc_dirs = [struct.mic(struct.positions[j] - origin)[0] for j in acc]
    acc_dirs = [v / np.linalg.norm(v) for v in acc_dirs if np.linalg.norm(v) > 1e-6]

    v_avoid = None
    if avoid_index is not None:
        v_avoid = struct.mic(struct.positions[avoid_index] - origin)[0]
        v_avoid /= np.linalg.norm(v_avoid)

    near = [j for j in np.where(d < 3.2)[0] if j not in (o_index, cov_index)]
    near_pos = struct.positions[near] if near else np.zeros((0, 3))
    near_h = np.array([elems[j] == "H" for j in near])

    best, best_score = None, -1e9
    for _ in range(n_trials):
        phi = rng.uniform(0, 2 * math.pi)
        local = np.array([math.sin(theta) * math.cos(phi),
                          math.sin(theta) * math.sin(phi),
                          -math.cos(theta)])
        h_dir = rot @ local
        h_pos = origin + h_dir * O_H

        rel = struct.mic(near_pos - h_pos) if len(near) else np.zeros((0, 3))
        dist = np.linalg.norm(rel, axis=1) if len(near) else np.array([])
        if dist.size and np.any(dist < np.where(near_h, 1.45, 1.55)):
            continue

        score = -float(np.sum(np.clip(2.0 - dist, 0, None))) if dist.size else 0.0
        if acc_dirs:
            score += 2.0 * max(float(np.dot(h_dir, a)) for a in acc_dirs)
        if v_avoid is not None:
            score -= 2.0 * float(np.dot(h_dir, v_avoid))
        if score > best_score:
            best, best_score = h_pos, score

    if best is None:
        raise RuntimeError(
            f"could not place hydroxyl hydrogen on oxygen {o_index} without "
            "a clash - the site is too crowded; check the structure")
    struct.extend(["H"], best[None, :])


def protonate_to_neutral(struct, verbose=True):
    """Close the remaining charge by protonating terminal oxygens.

    Only genuinely terminal oxygens (bonded to at most one Si) are
    candidates; a bridging Si-O-Si oxygen never takes a third bond. Among
    terminal candidates, the fewest additional (ionic) contacts identifies
    the site sticking furthest into the interlayer, which is where a real
    hydroxyl sits.
    """
    q = formal_charge(struct, ChargeModel())["total"]
    if q == 0:
        return 0
    if q > 0:
        raise RuntimeError(f"cell is {q:+d}; protonation only fixes a deficit")

    bonds = struct.neighbours()
    elems = struct.elements
    cand = []
    for i, e in enumerate(elems):
        if e != "O":
            continue
        if any(elems[j] == "H" for j in bonds[i]):
            continue
        n_si = sum(1 for j in bonds[i] if elems[j] == "Si")
        if n_si >= 2:
            continue
        heavy = [j for j in bonds[i] if elems[j] in FRAMEWORK_METALS or
                 elems[j] == "Si"]
        cand.append((len(heavy), i))
    cand.sort()

    added = 0
    for _, i in cand[: -q]:
        si_partners = [j for j in bonds[i] if elems[j] == "Si"]
        metal_partners = [j for j in bonds[i]
                          if elems[j] in FRAMEWORK_METALS and elems[j] != "Si"]
        if si_partners:
            nearest_metal = None
            if metal_partners:
                nearest_metal = min(
                    metal_partners,
                    key=lambda j: np.linalg.norm(
                        struct.mic(struct.positions[j] - struct.positions[i])[0]))
            _place_hydroxyl_h(struct, i, si_partners[0], avoid_index=nearest_metal)
        else:
            centre = (struct.positions[metal_partners].mean(axis=0)
                      if metal_partners else struct.positions[i])
            direction = struct.positions[i] - centre
            if np.linalg.norm(direction) < 1e-6:
                direction = np.array([0.0, 0.0, 1.0])
            _add_h(struct, i, direction, O_H)
        added += 1
    if verbose:
        print(f"  protonated {added} oxygen(s) to close {q:+d}")
    return added


def wrap_report(struct):
    """Assemble molecules for viewing (cosmetic only) and report what moved."""
    def split_members(s, groups):
        n = 0
        for g in groups:
            anchor = g.indices[0]
            for i in g.indices[1:]:
                raw = np.linalg.norm(s.positions[i] - s.positions[anchor])
                mic = np.linalg.norm(s.mic(s.positions[i] - s.positions[anchor])[0])
                if abs(raw - mic) > 0.01:
                    n += 1
        return n

    groups_before = detect_groups(struct)
    before_finite = split_members(struct, [g for g in groups_before if g.kind != "SiO4"])
    struct.wrap_molecules(groups_before)
    groups_after = detect_groups(struct)
    after_finite = split_members(struct, [g for g in groups_after if g.kind != "SiO4"])
    after_sio4 = split_members(struct, [g for g in groups_after if g.kind == "SiO4"])
    return {"finite_fixed": before_finite - after_finite,
            "finite_remaining": after_finite, "chain_seams": after_sio4}


def build(cif_name, form, out_stem, expect_ca_si, churakov_ref):
    print(f"\n=== tobermorite 11 A, {form} form ({cif_name}) ===")
    struct = build_structure(CIF_DIR / cif_name)
    print(f"  as read      : {len(struct):3d} atoms  {struct.composition()}  "
          f"q={formal_charge(struct)['total']:+d}")

    add_water_hydrogens(struct, verbose=False)
    print(f"  + declared H : {len(struct):3d} atoms  {struct.composition()}  "
          f"q={formal_charge(struct)['total']:+d}")

    protonate_to_neutral(struct)
    comp = struct.composition()
    ca_si = comp["Ca"] / comp["Si"]
    groups = {}
    for g in detect_groups(struct):
        groups[g.kind] = groups.get(g.kind, 0) + 1
    print(f"  neutral cell : {len(struct):3d} atoms  {comp}  "
          f"q={formal_charge(struct)['total']:+d}")
    print(f"  groups       : {groups}")
    print(f"  Ca/Si        : {ca_si:.4f}  (expected {expect_ca_si})")
    print(f"  close contacts <0.75 A: {len(close_contacts(struct, 0.75))}")
    assert abs(ca_si - expect_ca_si) < 0.01, \
        f"Ca/Si {ca_si:.4f} does not match the {form} formula ({expect_ca_si})"

    ortho = transform_cell(struct, [[2, 1, 0], [0, 2, 0], [0, 0, 1]], align=True)
    print(f"  orthogonal supercell: {len(ortho)} atoms  {ortho.composition()}  "
          f"q={formal_charge(ortho)['total']:+d}")
    print(f"  cell {np.round(ortho.lengths, 4).tolist()}  "
          f"angles {np.round(ortho.angles, 2).tolist()}")
    print(f"  Churakov reference : {churakov_ref}")

    bonds = ortho.neighbours()
    el = ortho.elements
    vecs = []
    for i, e in enumerate(el):
        if e != "Si":
            continue
        for o in [j for j in bonds[i] if el[j] == "O"]:
            vecs += [ortho.mic(ortho.positions[k] - ortho.positions[i])[0]
                     for k in bonds[o] if el[k] == "Si" and k != i]
    u = np.array([v / np.linalg.norm(v) for v in vecs])
    u = np.array([x if x[1] >= 0 else -x for x in u])
    along_y = sum(1 for x in u
                  if abs(x[1] - 1) < 0.05 and abs(x[0]) < 0.05 and abs(x[2]) < 0.05)
    print(f"  Si-O-Si links along +y: {along_y} of {len(u)}  (chain direction)")

    ortho.sort_by(2)
    report = wrap_report(ortho)
    print(f"  molecule assembly for viewing: {report['finite_fixed']} atom(s) "
          f"pulled together ({report['finite_remaining']} still split - should "
          f"be 0); {report['chain_seams']} SiO4 member(s) crossing the "
          "boundary (the chain continuing into the next cell, not a defect)")

    out_path = OUT_DIR / f"{out_stem}.xyz"
    write_xyz(out_path, ortho, wrap_molecules=True)
    print(f"  wrote {out_path}")
    return ortho


if __name__ == "__main__":
    build("tobermorite11A_anomalous_cod9005498.cif", "anomalous",
          "CSH_11A_anomalous", 0.667, "11.265 x 14.792 x 22.487, 90/90/90")
    build("tobermorite11A_normal_cod9005499.cif", "normal",
          "CSH_11A_normal", 0.75, "11.265 x 14.792 x 22.680, 90/90/90")
