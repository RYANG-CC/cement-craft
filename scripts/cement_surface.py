#!/usr/bin/env python3
"""
Surface construction: cut a slab along an axis, then heal and neutralise it.

Cutting a hydrated mineral leaves under-coordinated oxygens at both faces and,
almost always, a net charge. Neither can be ignored: a charged periodic cell is
not a physical model, and bare oxygens are not what the mineral exposes to
water. The reconstruction here follows the rule these phases actually obey:

    a surface oxygen bonded to ONE framework metal becomes a water molecule
    a surface oxygen bonded to TWO or more remains a hydroxyl

which is what returns every framework metal to its bulk coordination and, for a
stoichiometric cut, takes the slab to neutrality at the same time. Whatever
charge is left after that is reported and balanced explicitly, by adding or
removing protons, never absorbed silently.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cement_core import (  # noqa: E402
    Structure, detect_groups, element_of, close_contacts,
)
from cement_charge import ChargeModel, formal_charge  # noqa: E402
from cement_phases import O_H, WATER_OH, WATER_ANGLE, relax_overlaps  # noqa: E402

class SurfaceChargeError(RuntimeError):
    """Raised when a cut slab cannot be neutralised by surface protonation."""


FRAMEWORK_METALS = {
    "Mg", "Ca", "Al", "Fe", "Ni", "Zn", "Mn", "Co", "Cu", "Cr", "Ga",
    "Si", "Sc", "La", "Sr", "Ba", "Ti", "Zr",
}


def cut_slab(struct: Structure, axis: int = 2, thickness: float | None = None,
             origin: float = 0.0, vacuum: float = 15.0,
             keep_groups: bool = True) -> Structure:
    """Cut a slab of `thickness` starting at `origin` along `axis`.

    Molecular units are kept whole by default: a water or carbonate whose
    centre lies inside the slab is retained entirely, and one whose centre lies
    outside is removed entirely. Slicing on atom position alone would leave
    half-molecules at both faces, which is not a surface but a defect.
    """
    work = struct.copy().wrap()
    span = float(work.cell[axis, axis])
    thickness = thickness if thickness is not None else span

    lo, hi = origin, origin + thickness
    coords = work.positions[:, axis]

    if keep_groups:
        groups = detect_groups(work)
        in_group = {}
        for gi, g in enumerate(groups):
            for i in g.indices:
                in_group[i] = gi
        keep = np.zeros(len(work), dtype=bool)
        for i in range(len(work)):
            gi = in_group.get(i)
            ref = groups[gi].centre[axis] if gi is not None else coords[i]
            keep[i] = (lo <= ref < hi)
    else:
        keep = (coords >= lo) & (coords < hi)

    out = Structure(
        [l for l, k in zip(work.labels, keep) if k],
        work.positions[keep],
        work.cell.copy(),
        dict(work.info),
    )
    if len(out) == 0:
        raise ValueError("the slab window contains no atoms")

    # Re-zero along the cut axis and open the vacuum gap.
    out.positions[:, axis] -= out.positions[:, axis].min()
    new_len = out.positions[:, axis].max() + vacuum
    cell = out.cell.copy()
    cell[axis] = cell[axis] / np.linalg.norm(cell[axis]) * new_len
    out.cell = cell
    out.info.update({"slab_axis": axis, "slab_thickness": thickness,
                     "vacuum": vacuum})
    return out


def coordination(struct: Structure, bonds=None) -> dict[int, list[int]]:
    """Framework-metal neighbours of every oxygen."""
    bonds = bonds if bonds is not None else struct.neighbours()
    elems = struct.elements
    out = {}
    for i, e in enumerate(elems):
        if e != "O":
            continue
        out[i] = [j for j in bonds[i] if elems[j] in FRAMEWORK_METALS]
    return out


def _add_h(struct: Structure, o_index: int, direction: np.ndarray,
           length: float = O_H) -> None:
    d = direction / max(np.linalg.norm(direction), 1e-9)
    struct.extend(["H"], (struct.positions[o_index] + d * length)[None, :])


def protonate_surface(struct: Structure, axis: int = 2,
                      bulk_coordination: int | None = None,
                      verbose: bool = False) -> dict:
    """Heal a cut surface by protonating under-coordinated oxygens.

    Oxygens already carrying hydrogen keep it; the reconstruction only touches
    those the cut left bare. One framework metal neighbour means the oxygen was
    a bridging site that lost its partner, and it becomes a water molecule; two
    or more means it remains a hydroxyl.
    """
    bonds = struct.neighbours()
    elems = struct.elements
    coord = coordination(struct, bonds)

    counts = sorted(len(v) for v in coord.values())
    if bulk_coordination is None and counts:
        # The modal coordination in the interior is what "bulk" means here.
        bulk_coordination = max(set(counts), key=counts.count)

    added_h = 0
    to_water, to_hydroxyl = [], []
    for i, metals in coord.items():
        n_h = sum(1 for j in bonds[i] if elems[j] == "H")
        n_m = len(metals)
        if n_m == 0:
            continue                      # free water or a stray anion oxygen
        if n_m >= bulk_coordination:
            continue                      # fully coordinated: nothing to heal
        if n_h >= 2:
            continue                      # already water
        if n_m == 1 and n_h < 2:
            to_water.append((i, metals, n_h))
        elif n_h == 0:
            to_hydroxyl.append((i, metals))

    # Outward direction: away from the mean position of the bonded metals, and
    # biased along the surface normal so new hydrogens point into the vacuum.
    for i, metals, n_h in to_water:
        centre = struct.positions[list(metals)].mean(axis=0)
        out = struct.positions[i] - centre
        out[axis] += np.sign(out[axis] if out[axis] != 0 else 1.0) * 0.5
        need = 2 - n_h
        for k in range(need):
            perp = np.zeros(3)
            perp[(axis + 1) % 3] = (1.0 if k == 0 else -1.0)
            direction = out / max(np.linalg.norm(out), 1e-9) + 0.6 * perp
            _add_h(struct, i, direction, WATER_OH)
            added_h += 1

    for i, metals in to_hydroxyl:
        centre = struct.positions[list(metals)].mean(axis=0)
        out = struct.positions[i] - centre
        _add_h(struct, i, out, O_H)
        added_h += 1

    report = {
        "bulk_coordination": bulk_coordination,
        "to_water": len(to_water),
        "to_hydroxyl": len(to_hydroxyl),
        "hydrogens_added": added_h,
    }
    if verbose:
        print(f"  bulk O coordination taken as {bulk_coordination}")
        print(f"  singly coordinated O -> H2O: {len(to_water)}")
        print(f"  under-coordinated O  -> OH : {len(to_hydroxyl)}")
        print(f"  hydrogens added: {added_h}")
    return report


def neutralise(struct: Structure, model: ChargeModel | None = None,
               axis: int = 2, verbose: bool = False,
               counter_anion: str | None = None) -> dict:
    """Bring a slab to zero formal charge by adding or removing protons.

    Protons are the right currency: adding one converts a surface hydroxyl to
    water, removing one does the reverse, and both are chemistry the surface
    performs anyway in contact with water. Sites are chosen at the outermost
    faces, so the interior of the slab keeps its bulk protonation.
    """
    model = model or ChargeModel()
    q = formal_charge(struct, model)["total"]
    if q == 0:
        if verbose:
            print("  slab already neutral")
        return {"initial_charge": 0, "protons_added": 0, "protons_removed": 0,
                "counter_anions": 0, "final_charge": 0}

    coords = struct.positions[:, axis]
    mid = 0.5 * (coords.min() + coords.max())
    bonds = struct.neighbours()
    elems = struct.elements

    added = removed = 0
    if q < 0:
        # Too negative: protonate surface hydroxyls, outermost first.
        cand = []
        for i, e in enumerate(elems):
            if e != "O":
                continue
            hs = [j for j in bonds[i] if elems[j] == "H"]
            metals = [j for j in bonds[i] if elems[j] in FRAMEWORK_METALS]
            if len(hs) == 1 and metals:
                cand.append((abs(coords[i] - mid), i, metals))
        cand.sort(reverse=True)
        for _, i, metals in cand[: -q]:
            centre = struct.positions[list(metals)].mean(axis=0)
            direction = struct.positions[i] - centre
            perp = np.zeros(3)
            perp[(axis + 1) % 3] = 1.0
            _add_h(struct, i, direction / max(np.linalg.norm(direction), 1e-9)
                   + 0.6 * perp, WATER_OH)
            added += 1
    else:
        # Too positive: deprotonate surface water, outermost first.
        cand = []
        for g in detect_groups(struct):
            if g.kind != "H2O":
                continue
            o = [i for i in g.indices if elems[i] == "O"][0]
            metals = [j for j in bonds[o] if elems[j] in FRAMEWORK_METALS]
            if metals:
                h = [i for i in g.indices if elems[i] == "H"]
                cand.append((abs(coords[o] - mid), h[-1]))
        cand.sort(reverse=True)
        drop = [i for _, i in cand[:q]]
        struct.delete(drop)
        removed = len(drop)

    final = formal_charge(struct, model)["total"]
    if verbose:
        print(f"  charge {q:+d} -> {final:+d} "
              f"(+{added} H, -{removed} H)")
    counter_added = 0
    if final > 0 and counter_anion:
        # A positively charged hydroxide surface is balanced by anions, not by
        # stripping protons - which is what the mineral does in solution, and
        # what a cut through an LDH gallery leaves behind. Place them just
        # outside each face, where the interlayer species used to sit.
        counter_added = _add_counter_anions(struct, final, counter_anion, axis)
        final = formal_charge(struct, model)["total"]
        if verbose:
            print(f"  added {counter_added} {counter_anion} counter-ions "
                  f"-> {final:+d}")

    if final != 0:
        # Nothing left to try. For a layered phase this almost always means the
        # cut discarded interlayer anions that were balancing the sheets, and
        # no amount of surface chemistry replaces them - the cut itself has to
        # change, or counter-ions be supplied.
        raise SurfaceChargeError(
            f"slab still carries {final:+d} after healing "
            f"(started {q:+d}, +{added} H, -{removed} H, "
            f"{counter_added} counter-ions). For a layered phase cut normal "
            "to the stacking axis this usually means the window dropped "
            "interlayer anions: move `origin` so the cut falls inside a "
            "gallery, or pass counter_anion= to put them back.",
        )
    return {"initial_charge": q, "protons_added": added,
            "protons_removed": removed, "counter_anions": counter_added,
            "final_charge": final}


def _add_counter_anions(struct: Structure, charge: int, anion: str,
                        axis: int) -> int:
    """Place `charge`/|q| anions just beyond each face of the slab."""
    from cement_charge import ANION_CHARGE
    from cement_phases import FRAGMENTS, PACK_TARGET_HEAVY, PACK_TARGET_H

    q = abs(ANION_CHARGE[anion])
    n, residual = divmod(charge, q)
    if residual:
        raise SurfaceChargeError(
            f"residual charge {charge:+d} is not divisible by |{anion}| = {q}; "
            "choose a monovalent counter-anion or adjust the cut.")

    labels_f, geom = FRAGMENTS[anion]()
    coords = struct.positions[:, axis]
    lo, hi = coords.min(), coords.max()
    rng = np.random.default_rng(0)
    other = [k for k in range(3) if k != axis]
    placed = 0

    for k in range(n):
        face_hi = (k % 2 == 0)
        for _ in range(3000):
            centre = np.zeros(3)
            for ax in other:
                centre[ax] = rng.uniform(0, struct.cell[ax, ax])
            centre[axis] = (hi + rng.uniform(1.8, 3.4)) if face_hi \
                else (lo - rng.uniform(1.8, 3.4))
            trial = centre + geom
            d = np.linalg.norm(
                struct.mic(struct.positions[None, :, :] - trial[:, None, :]
                           ).reshape(-1, 3), axis=1)
            if d.min() >= PACK_TARGET_H + 0.2:
                struct.extend(labels_f, trial)
                placed += 1
                break
    if placed < n:
        raise SurfaceChargeError(
            f"could only place {placed}/{n} {anion} counter-ions at the surface")
    return placed





def build_surface(struct: Structure, axis: int = 2,
                  thickness: float | None = None, origin: float = 0.0,
                  vacuum: float = 15.0, oxidation: dict | None = None,
                  counter_anion: str | None = None,
                  verbose: bool = True) -> tuple[Structure, dict]:
    """Cut, protonate and neutralise in one call - the usual entry point."""
    model = ChargeModel(oxidation)
    slab = cut_slab(struct, axis=axis, thickness=thickness, origin=origin,
                    vacuum=vacuum)
    q0 = formal_charge(slab, ChargeModel(oxidation))["total"]
    if verbose:
        print(f"  cut slab: {len(slab)} atoms, charge {q0:+d} before healing")

    prot = protonate_surface(slab, axis=axis, verbose=verbose)
    neut = neutralise(slab, ChargeModel(oxidation), axis=axis,
                      verbose=verbose, counter_anion=counter_anion)

    # New hydrogens can land close to existing atoms; nudge them clear.
    h_added = [[i] for i in range(len(slab))
               if element_of(slab.labels[i]) == "H"]
    relax_overlaps(slab, h_added[-max(prot["hydrogens_added"], 1):],
                   min_heavy=1.9, min_h=1.4, n_iter=150)

    report = {"charge_after_cut": q0, **prot, **neut,
              "atoms": len(slab),
              "close_contacts": len(close_contacts(slab, 0.75))}
    slab.info.update({k: v for k, v in report.items()})
    slab.sort_by(axis)
    return slab, report


# ---------------------------------------------------------------------------
# Solid-solvent interface
# ---------------------------------------------------------------------------

# Bulk liquid water at ambient conditions, molecules per cubic angstrom.
WATER_NUMBER_DENSITY = 0.0334


def solvate(struct: Structure, axis: int = 2, density: float | None = None,
            n_water: int | None = None, gap: float = 2.6,
            counter_ion: str | None = None, n_counter_ion: int = 0,
            seed: int = 0, verbose: bool = False) -> dict:
    """Fill the vacuum above a slab with liquid water, making an interface.

    A slab in vacuum answers a different question from a slab in contact with
    water, and for cement phases it is nearly always the wetted surface that
    matters - the hydroxyls reorganise, the interlayer exchanges, and adsorbed
    species are solvated. This fills the empty region at bulk water density
    rather than by a molecule count, so the reservoir is physically meaningful
    whatever the slab area happens to be.

    `gap` is the clearance left between the slab surface and the first water,
    so the reservoir starts at a hydrogen-bond distance rather than on top of
    the surface hydroxyls.

    Counter-ions may be added alongside, for a slab whose surface chemistry
    leaves it charged in solution.
    """
    from cement_phases import FRAGMENTS, _fill_interlayer

    rng = np.random.default_rng(seed)
    coords = struct.positions[:, axis]
    top = float(coords.max())
    cell_len = float(struct.cell[axis, axis])
    free = cell_len - top - gap
    if free <= 2.0:
        raise ValueError(
            f"only {free:.1f} A of vacuum above the slab along axis {axis}; "
            "rebuild the slab with a larger `vacuum` before solvating")

    other = [k for k in range(3) if k != axis]
    area = float(struct.cell[other[0], other[0]] * struct.cell[other[1], other[1]])
    volume = area * free

    if n_water is None:
        density = density if density is not None else WATER_NUMBER_DENSITY
        n_water = int(round(volume * density))
    if n_water <= 0:
        raise ValueError("computed zero water molecules; check the vacuum size")

    inserted: list[list[int]] = []
    z_mid = top + gap + free / 2.0
    z_half = free / 2.0

    if counter_ion and n_counter_ion:
        if counter_ion not in FRAGMENTS:
            raise ValueError(f"no fragment template for {counter_ion!r}")
        _fill_interlayer(struct, counter_ion, n_counter_ion, z_mid, z_half,
                         rng, inserted=inserted)

    _fill_interlayer(struct, "H2O", n_water, z_mid, z_half, rng,
                     inserted=inserted)
    remaining = relax_overlaps(struct, inserted)

    report = {
        "reservoir_thickness": round(free, 2),
        "reservoir_volume": round(volume, 1),
        "n_water": n_water,
        "water_density": round(n_water / volume, 5),
        "counter_ions": n_counter_ion,
        "residual_clashes": remaining,
    }
    if verbose:
        print(f"  reservoir {free:.1f} A thick, {volume:.0f} A^3")
        print(f"  {n_water} H2O at {n_water / volume:.4f} /A^3 "
              f"(bulk water is {WATER_NUMBER_DENSITY})")
        if n_counter_ion:
            print(f"  {n_counter_ion} {counter_ion} counter-ions")
    struct.info.update(report)
    return report


def strip_interlayer(struct: Structure, axis: int = 2,
                     keep_water: int = 0, model: ChargeModel | None = None,
                     verbose: bool = False) -> dict:
    """Remove interlayer species from a cut slab, staying charge neutral.

    Cutting one layer out of a layered phase leaves the interlayer contents
    that happened to fall inside the window - loosely bound cations and water
    that belonged to the gallery, not to the sheet. Stripping them gives the
    "simplified basal surface" that is often what is wanted: a single layer
    presenting its own termination, with no gallery chemistry confusing the
    picture.

    Neutrality is maintained by construction. Free water is uncharged and can
    go in any number; every interlayer cation removed takes its charge with it,
    so an equal charge of hydroxide leaves with it. The routine refuses rather
    than returning a charged slab if it cannot pair them up.

    `keep_water` retains that many of the most strongly bound waters.
    """
    model = model or ChargeModel()
    bonds = struct.neighbours()
    elems = struct.elements
    groups = detect_groups(struct)

    # Framework = the connected silicate/metal-oxide network. An oxygen bonded
    # to Si or Al is framework; a cation bonded only to such oxygens through
    # fewer than three contacts is interlayer rather than structural.
    framework_o = {
        i for i, e in enumerate(elems)
        if e == "O" and any(elems[j] in ("Si", "Al", "S") for j in bonds[i])
    }

    free_water, free_oh = [], []
    for g in groups:
        o = next((i for i in g.indices if elems[i] == "O"), None)
        if o is None or o in framework_o:
            continue
        if g.kind == "H2O":
            free_water.append(g)
        elif g.kind == "OH":
            n_metal = sum(1 for j in bonds[o] if elems[j] in FRAMEWORK_METALS)
            if n_metal <= 1:
                free_oh.append(g)

    # Interlayer cations: bonded to few framework oxygens.
    interlayer_cations = []
    for i, e in enumerate(elems):
        if e not in ("Ca", "Na", "K", "Mg", "Sr", "Ba"):
            continue
        n_fw = sum(1 for j in bonds[i] if j in framework_o)
        if n_fw < 3:
            interlayer_cations.append(i)

    # Order water by how weakly it is held, so the retained ones are the bound ones.
    fw_pos = struct.positions[sorted(framework_o)] if framework_o else None
    def isolation(g):
        if fw_pos is None:
            return 0.0
        return float(np.linalg.norm(struct.mic(fw_pos - g.centre), axis=1).min())
    free_water.sort(key=isolation, reverse=True)
    drop_water = free_water[: max(len(free_water) - keep_water, 0)]

    # Pair each cation's charge against hydroxide so the removal is neutral.
    cation_charge = sum(model.state(elems[i]) for i in interlayer_cations)
    n_oh_needed = cation_charge
    # Short of hydroxide, balance by cation/proton exchange instead: an
    # interlayer Ca(2+) leaving is replaced by two H(+) on framework oxygens.
    # That is what a real basal surface does in contact with water, and it
    # keeps the sheet itself intact rather than eroding it.
    protonate_deficit = max(n_oh_needed - len(free_oh), 0)
    n_oh_needed = min(n_oh_needed, len(free_oh))
    drop_oh = free_oh[:n_oh_needed]
    protonate_n = 0

    doomed = set(interlayer_cations)
    for g in drop_water + drop_oh:
        doomed.update(g.indices)

    # Choose protonation sites before deleting, while indices are still valid:
    # framework oxygens nearest the departing cations, and least protonated.
    proton_sites = []
    if protonate_deficit:
        cand = []
        for i in sorted(framework_o):
            if i in doomed:
                continue
            if any(elems[j] == "H" for j in bonds[i]):
                continue
            if interlayer_cations:
                d = min(float(np.linalg.norm(
                    struct.mic(struct.positions[c] - struct.positions[i])[0]))
                    for c in interlayer_cations)
            else:
                d = 0.0
            cand.append((d, i))
        cand.sort()
        proton_sites = [i for _, i in cand[:protonate_deficit]]
        if len(proton_sites) < protonate_deficit:
            raise SurfaceChargeError(
                f"need {protonate_deficit} framework oxygens to protonate but "
                f"only {len(proton_sites)} are free; cut a different window.")
        for i in proton_sites:
            metals = [j for j in bonds[i] if elems[j] in FRAMEWORK_METALS]
            centre = (struct.positions[metals].mean(axis=0) if metals
                      else struct.positions[i])
            direction = struct.positions[i] - centre
            if np.linalg.norm(direction) < 1e-6:
                direction = np.array([0.0, 0.0, 1.0])
            _add_h(struct, i, direction, O_H)
        protonate_n = len(proton_sites)

    struct.delete(sorted(doomed))

    q = formal_charge(struct, model)["total"]
    report = {
        "cations_removed": len(interlayer_cations),
        "water_removed": len(drop_water),
        "hydroxide_removed": len(drop_oh),
        "framework_protonated": protonate_n,
        "atoms": len(struct),
        "final_charge": q,
    }
    if verbose:
        print(f"  removed {len(interlayer_cations)} interlayer cation(s), "
              f"{len(drop_water)} H2O, {len(drop_oh)} OH; protonated "
              f"{protonate_n} framework O -> {len(struct)} atoms, charge {q:+d}")
    if q != 0:
        raise SurfaceChargeError(
            f"stripping the interlayer left charge {q:+d}; the cation/hydroxide "
            "pairing did not balance. Inspect the slab before using it.")
    return report
