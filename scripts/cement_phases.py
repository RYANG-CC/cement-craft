#!/usr/bin/env python3
"""
Builders for the cement phases: LDH, AFm, AFt, hydrogarnet, C-S-H.

Two construction routes, chosen per phase for good reason:

* **LDH is built procedurally.** Its structure is a brucite-like M(OH)2 sheet
  with some M(II) replaced by M(III), and the interlayer filled with whatever
  anion balances the resulting charge. Generating that from geometry rather
  than editing a fixed cell is what makes arbitrary cation pairs, ratios,
  anions and water contents reachable.

* **AFm, AFt, hydrogarnet and C-S-H start from a reference cell.** Their
  frameworks are rigid and crystallographically specific; rebuilding them from
  scratch would invent coordinates. These are read from `assets/structures/`
  and then substituted, hydrated or expanded.

Every builder returns a `Structure` and records how it was made in `.info`.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cement_core import (  # noqa: E402
    Structure, cell_from_parameters, detect_groups, element_of,
    read_cif, read_xyz, close_contacts,
)
from cement_charge import ANION_CHARGE, ChargeModel, formal_charge  # noqa: E402

ASSETS = Path(__file__).resolve().parent.parent / "assets" / "structures"

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------

# In-plane metal-metal distance of the brucite-like sheet, angstrom. Measured
# 3.07-3.09 for Mg-Al LDH; it shifts with the cation pair, so it is a parameter.
A_HEX_DEFAULT = 3.07

# Basal spacing (metal sheet to metal sheet) by interlayer anion, angstrom.
# These set the interlayer volume and therefore how much water fits.
BASAL_SPACING = {
    "CO3": 7.78,
    "SO4": 8.90,
    "Cl": 7.80,
    "OH": 7.60,
    "NO3": 8.80,
}

M_O_DISTANCE = {"Mg": 2.07, "Ca": 2.35, "Fe": 2.10, "Ni": 2.05, "Zn": 2.08,
                "Mn": 2.15, "Co": 2.08, "Cu": 2.10,
                "Al": 1.90, "Cr": 1.98, "Ga": 1.99, "Sc": 2.10, "La": 2.40}
O_H = 0.97
WATER_OH = 0.96
WATER_ANGLE = 104.5

# Packing separations, angstrom. These structures are starting points for DFT
# or MD relaxation, not final geometries, so the bar is "no unphysical
# overlap", not "equilibrium distances". TARGET is what packing aims for -
# roughly a hydrogen-bonded contact; FLOOR is the hard limit below which a
# structure is rejected, because anything closer will not relax out cleanly.
PACK_TARGET_HEAVY, PACK_TARGET_H = 2.45, 1.55
PACK_FLOOR_HEAVY, PACK_FLOOR_H = 2.20, 1.30


# ---------------------------------------------------------------------------
# Molecular fragment templates
# ---------------------------------------------------------------------------

def _carbonate() -> tuple[list[str], np.ndarray]:
    """Planar CO3 2-, C-O 1.29 A at 120 degrees."""
    r = 1.29
    pos = [[0, 0, 0]]
    labels = ["C"]
    for k in range(3):
        t = math.radians(90 + 120 * k)
        pos.append([r * math.cos(t), r * math.sin(t), 0.0])
        labels.append("O")
    return labels, np.array(pos)


def _sulfate() -> tuple[list[str], np.ndarray]:
    """Tetrahedral SO4 2-, S-O 1.48 A."""
    r = 1.48 / math.sqrt(3)
    pos = [[0, 0, 0]]
    labels = ["S"]
    for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
        pos.append([r * v[0], r * v[1], r * v[2]])
        labels.append("O")
    return labels, np.array(pos)


def _water() -> tuple[list[str], np.ndarray]:
    half = math.radians(WATER_ANGLE / 2)
    return (
        ["O", "H", "H"],
        np.array([
            [0.0, 0.0, 0.0],
            [WATER_OH * math.sin(half), 0.0, WATER_OH * math.cos(half)],
            [-WATER_OH * math.sin(half), 0.0, WATER_OH * math.cos(half)],
        ]),
    )


def _hydroxide() -> tuple[list[str], np.ndarray]:
    return ["O", "H"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, O_H]])


FRAGMENTS = {
    "CO3": _carbonate, "SO4": _sulfate, "H2O": _water, "OH": _hydroxide,
    "Cl": lambda: (["Cl"], np.array([[0.0, 0.0, 0.0]])),
}


def _random_rotation(rng) -> np.ndarray:
    """Uniform random rotation matrix via QR of a Gaussian matrix."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.diag(r))


def _planar_rotation(rng, max_tilt_deg: float = 20.0) -> np.ndarray:
    """Free spin about z, plus a small tilt.

    Planar interlayer anions (carbonate, nitrate) sit with their plane roughly
    parallel to the hydroxide sheets - the gallery is barely wider than the
    anion, so a freely oriented one simply does not fit, and the flat-lying
    arrangement is what diffraction shows in any case."""
    a = rng.uniform(0, 2 * math.pi)
    cz, sz = math.cos(a), math.sin(a)
    spin = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    t = math.radians(rng.uniform(-max_tilt_deg, max_tilt_deg))
    b = rng.uniform(0, 2 * math.pi)
    axis = np.array([math.cos(b), math.sin(b), 0.0])
    K = np.array([[0.0, -axis[2], axis[1]],
                  [axis[2], 0.0, -axis[0]],
                  [-axis[1], axis[0], 0.0]])
    tilt = np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * (K @ K)
    return tilt @ spin


# Anions whose shape forces a flat-lying orientation in a narrow gallery.
PLANAR_SPECIES = {"CO3", "NO3"}


# ---------------------------------------------------------------------------
# LDH
# ---------------------------------------------------------------------------

def build_ldh(
    divalent: str = "Mg",
    trivalent: str = "Al",
    ratio: float = 2.0,
    anion: str = "CO3",
    n_a: int = 4,
    n_b: int = 6,
    n_layers: int = 2,
    water_per_anion: float = 4.0,
    a_hex: float | None = None,
    basal: float | None = None,
    seed: int = 0,
    label_trivalent: str | None = None,
    auto_basal: bool = True,
) -> Structure:
    """Build a layered double hydroxide.

    `ratio` is M(II):M(III) on the metal sublattice, so 2.0 means Mg2Al1. The
    sheet charge is one per trivalent cation, and the interlayer anion count
    follows from it - the structure is neutral by construction, not by
    adjustment afterwards.

    The orthogonal cell is a n_a x n_b supercell of the hexagonal sheet, giving
    2*n_a*n_b metals per layer:
        a = n_a * a_hex * sqrt(3)
        b = n_b * a_hex

    With `auto_basal`, the gallery is allowed to expand when the requested
    water does not fit at the tabulated spacing. That is the real behaviour of
    these materials - basal spacing is set by hydration state, and a swollen
    gallery is a physical structure rather than a fudge - but the amount of
    swelling is reported so it can be checked against measurement.
    """
    if trivalent == divalent:
        raise ValueError("divalent and trivalent cations must differ")
    if anion not in ANION_CHARGE:
        raise ValueError(f"unknown anion {anion!r}; known: {sorted(ANION_CHARGE)}")

    a_hex = a_hex or A_HEX_DEFAULT
    basal_requested = basal or BASAL_SPACING.get(anion, 7.8)
    max_swell = 3.0 if auto_basal else 0.0
    last_error: Exception | None = None

    swell = 0.0
    while swell <= max_swell + 1e-9:
        try:
            return _build_ldh_once(
                divalent, trivalent, ratio, anion, n_a, n_b, n_layers,
                water_per_anion, a_hex, basal_requested + swell, seed,
                label_trivalent, basal_requested)
        except RuntimeError as exc:
            last_error = exc
            swell += 0.3
    raise RuntimeError(
        f"{last_error} Gallery could not be filled even after swelling the "
        f"basal spacing by {max_swell:.1f} A. Reduce water_per_anion.")


def _build_ldh_once(divalent, trivalent, ratio, anion, n_a, n_b, n_layers,
                    water_per_anion, a_hex, basal, seed, label_trivalent,
                    basal_requested):
    rng = np.random.default_rng(seed)

    per_layer = 2 * n_a * n_b
    n_metal = per_layer * n_layers
    n_tri_total = int(round(n_metal / (ratio + 1.0)))
    if n_tri_total == 0:
        raise ValueError("ratio too high for this cell: no trivalent cations fit")
    # Distribute trivalent cations evenly between layers.
    per_layer_tri = [n_tri_total // n_layers] * n_layers
    for k in range(n_tri_total % n_layers):
        per_layer_tri[k] += 1

    a = n_a * a_hex * math.sqrt(3.0)
    b = n_b * a_hex
    c = n_layers * basal
    cell = cell_from_parameters(a, b, c, 90.0, 90.0, 90.0)

    # Brucite sheet geometry. The metals form a triangular lattice; the
    # hydroxyl oxygens sit above and below the centroids of alternating
    # triangles, which is what makes every metal six-coordinate octahedral.
    # Distance from a metal to each of its six oxygens projected in-plane is
    # a_hex/sqrt(3), so the out-of-plane offset follows from the M-O bond.
    root3 = math.sqrt(3.0)
    d_in = a_hex / root3
    m_o = M_O_DISTANCE.get(divalent, 2.07)
    dz = math.sqrt(max(m_o ** 2 - d_in ** 2, 0.25))

    # Within one orthogonal subcell (a_hex*sqrt3 by a_hex) there are two metal
    # sites, two upper-hollow oxygens and two lower-hollow oxygens.
    metal_basis = [
        np.array([0.0, 0.0, 0.0]),
        np.array([a_hex * root3 / 2.0, a_hex / 2.0, 0.0]),
    ]
    up_offset = np.array([a_hex * root3 / 6.0, a_hex / 2.0, 0.0])
    down_offset = np.array([a_hex * root3 / 3.0, 0.0, 0.0])

    labels: list[str] = []
    pos: list[list[float]] = []
    tri_label = label_trivalent or trivalent

    for layer in range(n_layers):
        z0 = layer * basal
        sites = []
        for i in range(n_a):
            for j in range(n_b):
                base = np.array([i * a_hex * root3, j * a_hex, 0.0])
                for m in metal_basis:
                    sites.append(base + m)
        sites = np.array(sites)

        order = rng.permutation(len(sites))
        tri_idx = set(int(k) for k in order[: per_layer_tri[layer]])

        for k, s in enumerate(sites):
            labels.append(tri_label if k in tri_idx else divalent)
            pos.append([s[0], s[1], z0])

        # Hydroxyls: upper sheet on one hollow set, lower sheet on the other.
        for sign, offset in ((+1.0, up_offset), (-1.0, down_offset)):
            for s in sites:
                ox = s + offset + np.array([0.0, 0.0, sign * dz])
                labels.append("O")
                pos.append([ox[0], ox[1], z0 + sign * dz])
                # H points out of the sheet, into the gallery.
                labels.append("H")
                pos.append([ox[0], ox[1], z0 + sign * (dz + O_H)])

    struct = Structure(labels, np.array(pos), cell)
    inserted: list[list[int]] = []

    # Interlayer contents: anions to balance, then water.
    q_anion = abs(ANION_CHARGE[anion])
    n_anion_total, residual = divmod(n_tri_total, q_anion)
    if residual:
        raise ValueError(
            f"sheet charge {n_tri_total}+ is not divisible by |{anion}| = {q_anion}. "
            f"Adjust the ratio or cell size so the trivalent count is a multiple "
            f"of {q_anion}, or use a monovalent anion."
        )
    n_water_total = int(round(n_anion_total * water_per_anion))

    per_gap_anion = [n_anion_total // n_layers] * n_layers
    for k in range(n_anion_total % n_layers):
        per_gap_anion[k] += 1
    per_gap_water = [n_water_total // n_layers] * n_layers
    for k in range(n_water_total % n_layers):
        per_gap_water[k] += 1

    for layer in range(n_layers):
        # Interlayer midplane sits between this sheet and the next.
        mid = layer * basal + basal / 2.0
        # The usable gallery runs between the two hydroxyl-H planes, less
        # room for the inserted molecule itself. Seeding anything below that
        # boundary is unrecoverable: framework atoms never move, so a molecule
        # started inside the sheet stays clashing however long it is relaxed.
        half = (basal / 2.0) - (dz + O_H) - 0.20
        half = max(half, 0.4)
        _fill_interlayer(struct, anion, per_gap_anion[layer], mid, half, rng,
                         inserted=inserted)
        _fill_interlayer(struct, "H2O", per_gap_water[layer], mid, half, rng,
                         inserted=inserted)

    relax_overlaps(struct, inserted)
    # Accept at the physical floor rather than the packing target: a contact a
    # little short of ideal relaxes out in the first few DFT steps, whereas a
    # genuine overlap does not.
    hard = relax_overlaps(struct, inserted, min_heavy=PACK_FLOOR_HEAVY,
                          min_h=PACK_FLOOR_H, n_iter=0)
    if hard:
        raise RuntimeError(
            f"{hard} contacts remain below the physical floor "
            f"({PACK_FLOOR_HEAVY} A heavy / {PACK_FLOOR_H} A with H). The "
            "interlayer is over-filled for this basal spacing - reduce "
            "water_per_anion or increase basal.")

    struct.wrap()
    struct.info.update({
        "phase": "LDH",
        "divalent": divalent, "trivalent": tri_label, "ratio": ratio,
        "anion": anion, "n_anion": n_anion_total, "n_water": n_water_total,
        "sheet_charge": n_tri_total, "layers": n_layers,
        "a_hex": a_hex, "basal": round(basal, 3),
        "basal_requested": round(basal_requested, 3),
        "basal_swelling": round(basal - basal_requested, 3),
    })
    return struct


def _fill_interlayer(struct: Structure, species: str, count: int,
                     z_mid: float, z_half: float, rng, tries: int = 3000,
                     inserted: list | None = None) -> None:
    """Insert `count` copies of a fragment into an interlayer gallery.

    Candidate sites come from a jittered grid spanning the gallery rather than
    from uniform random sampling. At the densities a real LDH interlayer
    reaches - roughly four waters per carbonate - pure rejection sampling
    stalls well short of the target, because late insertions almost always
    land on top of something. A shuffled grid keeps the disorder that matters
    (position jitter and free orientation) while guaranteeing the candidates
    are spread out to begin with.
    """
    if count <= 0:
        return
    make = FRAGMENTS[species]
    labels_f, geom = make()
    cell_a, cell_b = struct.cell[0, 0], struct.cell[1, 1]
    frag_h = np.array([element_of(l) == "H" for l in labels_f])
    # Insert at the packing target. A permissive threshold reaches the count
    # quickly but seeds molecules against the immobile sheet, where no amount
    # of relaxation frees them; density is instead reached by alternating
    # insertion with relaxation, which opens genuine space.
    insert_heavy, insert_h = PACK_TARGET_HEAVY, PACK_TARGET_H
    inv_cell = np.linalg.inv(struct.cell)
    existing_h = np.array([element_of(l) == "H" for l in struct.labels])

    # Grid sized so there are comfortably more candidate sites than fragments.
    target = max(count * 4, 8)
    nz = max(1, int(round(2 * z_half / 2.8)))
    per_plane = math.ceil(target / nz)
    aspect = cell_a / cell_b
    nx = max(1, int(round(math.sqrt(per_plane * aspect))))
    ny = max(1, math.ceil(per_plane / nx))

    candidates = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                candidates.append((
                    (ix + 0.5) * cell_a / nx,
                    (iy + 0.5) * cell_b / ny,
                    z_mid - z_half + (iz + 0.5) * (2 * z_half) / nz,
                ))
    candidates = np.array(candidates)
    rng.shuffle(candidates)

    jitter = np.array([cell_a / nx, cell_b / ny, 2 * z_half / nz]) * 0.35
    placed = 0
    attempts_per_site = max(4, tries // max(len(candidates), 1))
    mine: list[list[int]] = []

    for round_ in range(10):
      if placed >= count:
          break
      if round_:
          # Nothing more fits at this arrangement; let what is already in the
          # gallery spread out, then try the remaining molecules again.
          relax_overlaps(struct, (inserted if inserted is not None else mine),
                         n_iter=120, step=0.3)
          rng.shuffle(candidates)
      for base in candidates:
        if placed >= count:
            break
        for _ in range(attempts_per_site):
            centre = base + rng.uniform(-1, 1, size=3) * jitter
            centre[2] = np.clip(centre[2], z_mid - z_half, z_mid + z_half)
            if len(geom) == 1:
                rot = np.eye(3)
            elif species in PLANAR_SPECIES:
                rot = _planar_rotation(rng)
            else:
                rot = _random_rotation(rng)
            trial = centre + geom @ rot.T

            # One vectorised test over the whole fragment against the whole
            # structure. Rebuilding an element list per trial, as an earlier
            # version did, dominated the runtime once the gallery filled up.
            delta = struct.positions[None, :, :] - trial[:, None, :]
            frac = delta @ inv_cell
            frac -= np.round(frac)
            d = np.linalg.norm(frac @ struct.cell, axis=2)
            # Heavy-heavy contacts must clear the packing target; anything
            # involving a hydrogen may approach to a hydrogen-bond distance.
            limit = np.where(existing_h[None, :] | frag_h[:, None],
                             insert_h, insert_heavy)
            ok = not np.any(d < limit)
            if ok:
                start = len(struct)
                struct.extend(labels_f, trial)
                group = list(range(start, len(struct)))
                existing_h = np.concatenate([existing_h, frag_h])
                mine.append(group)
                if inserted is not None:
                    inserted.append(group)
                placed += 1
                break

    if placed < count:
        raise RuntimeError(
            f"could only place {placed}/{count} {species} in the interlayer. "
            "The gallery is too tight - increase the basal spacing, reduce the "
            "water content, or enlarge the cell."
        )


def relax_overlaps(struct: Structure, mobile: list[list[int]],
                   min_heavy: float = PACK_TARGET_HEAVY,
                   min_h: float = PACK_TARGET_H,
                   n_iter: int = 600, step: float = 0.20) -> int:
    """Push overlapping inserted molecules apart, treating each as rigid.

    Filling an interlayer to the density a real LDH reaches means accepting
    some overlap at insertion and removing it afterwards - the same trade a
    dedicated packing code makes. Molecules translate but never distort, so
    bond lengths and angles stay exactly as the fragment templates defined
    them. Framework atoms are immobile, so the sheet is never disturbed.

    Returns the number of remaining clashes."""
    if not mobile:
        return 0

    heavy = np.array([element_of(l) != "H" for l in struct.labels])
    n_atoms = len(struct)
    # Flat arrays instead of nested loops: the whole mobile-vs-all distance
    # matrix is small enough to evaluate at once, and doing so turns a
    # minutes-long relaxation into a fraction of a second.
    mob_idx = np.array([i for g in mobile for i in g], dtype=int)
    mob_grp = np.array([gi for gi, g in enumerate(mobile) for _ in g], dtype=int)
    n_groups = len(mobile)

    group_of = np.full(n_atoms, -1, dtype=int)
    group_of[mob_idx] = mob_grp

    limits = np.where(heavy[mob_idx][:, None] & heavy[None, :], min_heavy, min_h)
    # A pair inside one molecule is a bond, not a clash.
    same_mol = group_of[None, :] == mob_grp[:, None]
    inv_cell = np.linalg.inv(struct.cell)

    def _clash_state():
        delta = struct.positions[None, :, :] - struct.positions[mob_idx][:, None, :]
        frac = delta @ inv_cell
        frac -= np.round(frac)
        delta = frac @ struct.cell
        r = np.linalg.norm(delta, axis=2)
        bad = (r < limits) & (r > 1e-6) & (~same_mol)
        return delta, r, bad

    for _ in range(n_iter):
        delta, r, bad = _clash_state()
        if not bad.any():
            return 0
        overlap = np.where(bad, limits - r, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            unit = np.where(r[:, :, None] > 1e-6, -delta / r[:, :, None], 0.0)
        push = (unit * overlap[:, :, None]).sum(axis=1)

        shifts = np.zeros((n_groups, 3))
        np.add.at(shifts, mob_grp, push)
        shifts *= step
        mag = np.linalg.norm(shifts, axis=1, keepdims=True)
        big = (mag > 0.3).ravel()
        shifts[big] *= (0.3 / mag[big])
        struct.positions[mob_idx] += shifts[mob_grp]

    _, _, bad = _clash_state()
    return int(bad.sum())

# ---------------------------------------------------------------------------
# Reference-cell phases
# ---------------------------------------------------------------------------

REFERENCE = {
    "AFm": {
        "file": "AFm_SO4.xyz",
        "cell": (5.7586002350, 9.9741878510, 26.7945995331, 90.0, 90.0, 90.0),
        "formula": "Ca4Al2(SO4)(OH)12.4H2O x3",
        "note": "monosulfate AFm; interlayer anion SO4",
    },
    "AFt": {
        "file": "AFt_SO4.xyz",
        "cell": (11.2290, 11.2290, 21.4780, 90.0, 90.0, 120.0),
        "formula": "Ca6Al2(SO4)3(OH)12.24H2O x2",
        "note": "ettringite; columnar Ca/Al with SO4 and water in channels",
    },
    "CSH": {
        "file": "CSH_14A.xyz",
        "cell": (14.85, 13.47, 27.987, 90.0, 90.0, 90.0),
        "formula": "Ca5Si6O16(OH)2.6H2O x8",
        "note": "14 A tobermorite-like C-S-H, Ca/Si = 0.83",
    },
    "hydrogarnet": {
        "file": "hydrogarnet_katoite.cif",
        "cell": None,          # taken from the CIF
        "formula": "Ca3Al2(OH)12 x8",
        "note": "katoite, cubic Ia-3d; Si-free hydrogarnet end member",
    },
}


# ---------------------------------------------------------------------------
# C-S-H reference series
# ---------------------------------------------------------------------------
#
# C-S-H is not one structure but a compositional family, and Ca/Si is the
# variable that matters most: it sets silicate chain length, and 0.67, 0.83 and
# ~1.25 are all common. Rather than derive one ratio from another by removing
# bridging tetrahedra - a modelling choice with several published conventions -
# each ratio is served by its own reference cell.
#
# The tobermorite series these sit on (Churakov 2009a, Eur. J. Mineral. 21,
# 261-271; Churakov 2009b, Am. Mineral. 94, 156-165):
#
#   14 A plombierite    Ca5Si6O16(OH)2.7H2O    Ca/Si 0.833
#   11 A normal         Ca5-xSi6O17-2x(OH)2x.5H2O  Ca/Si 0.75 as simulated
#   11 A anomalous      Ca4Si6O15(OH)2.5H2O    Ca/Si 0.667
#    9 A riversideite   Ca5Si6O16(OH)2         Ca/Si 0.833
#
# Churakov's validated orthorhombic AIMD supercells, for reference when
# building an 11 A cell from crystallography rather than from these files:
#   anomalous 11 A   11.265 x 14.792 x 22.487 A   (90/90/90)
#   normal    11 A   11.265 x 14.792 x 22.680 A   (90/90/90)
# both being 2 x 2 x 1 of the monoclinic B11m cell (13.47 x 14.74 x 22.487,
# gamma 123.25).

CSH_REFERENCES = {
    0.667: {
        "file": "CSH_11A_anomalous.xyz",
        "cell": (11.2648, 14.7700, 22.4870, 90.0, 90.0, 90.0),
        "formula": "Ca32Si48O176H96 (2x2x1 of COD 9005498, Wessels mine)",
        "note": "11 A anomalous tobermorite; Ca2 interlayer site absent, not "
                "disordered - see references/phases.md",
        "cell_verified": True,
    },
    0.75: {
        "file": "CSH_11A_normal.xyz",
        "cell": (11.2688, 14.7380, 22.6800, 90.0, 90.0, 90.0),
        "formula": "Ca36Si48O176H88 (2x2x1 of COD 9005499, Basenov, ordering A)",
        "note": "11 A normal tobermorite; Ca2/Wat1/Wat3 disorder resolved as "
                "one explicit ordering ('A') - see references/phases.md",
        "cell_verified": True,
    },
    0.833: {
        "file": "CSH_14A.xyz",
        "cell": (14.85, 13.47, 27.987, 90.0, 90.0, 90.0),
        "formula": "Ca5Si6O16(OH)2.6H2O x8",
        "note": "14 A tobermorite (plombierite); 2x2x1 of a=6.735 b=7.425",
        "cell_verified": True,
    },
    1.25: {
        "file": "CSH_v_CaSi125.xyz",
        "cell": (11.4769, 15.0660, 19.2663, 90.0, 90.0, 90.0),
        "formula": "Ca40Si32O144H80",
        "note": "defect C-S-H with bridging-tetrahedra vacancies",
        "cell_verified": True,
    },
}

# `CSH_CaSi067.xyz` in the archive is deliberately NOT in this registry, even
# though it also has Ca/Si 0.667. It is a single layer cut from the 0.833
# phase with its interlayer stripped - a simplified basal surface, not a bulk
# phase - and its Ca/Si is a consequence of removing interlayer calcium, not
# a third bulk stoichiometry. The bulk 0.667 phase is the anomalous 11 A
# tobermorite above. Reproduce the slab with:
#
#     slab = cut_slab(build_csh(ca_si=0.833), axis=2, thickness=c/2)
#     strip_interlayer(slab, axis=2, keep_water=10)
#
# which returns Ca16 Si24 O80 H32, neutral - the archive structure exactly.
# Its in-plane cell is the parent's; the c axis is chosen for the water and
# vacuum wanted.


def csh_chain_report(struct: Structure) -> dict:
    """Ca/Si and silicate-chain connectivity - check for every CSH structure.

    Qn = how many bridging (Si-O-Si) oxygens a tetrahedron has. In a
    periodic, undefected chain there are no chain *ends* (the repeat
    continues into the next cell), so an intact single dreierkette is
    all-Q2, not Q1/Q2 - Q1 is what a genuinely truncated chain (a removed
    bridging tetrahedron, an amorphous surface) looks like. A single chain
    stays 2-connected throughout; a **double** dreierkette (two chains
    condensed at the bridging tetrahedron) shows Q3 at exactly the bridging
    site, 1 Si in 3 - the same distinction the tobermorite 9 A validation
    drew qualitatively (all-Q2 single chains there; see
    csh-tobermorite9a-validated.md), now a routine quantitative check
    instead of a one-off one.

    Mean chain length uses the standard silicate-NMR convention
    MCL = 2*(Q1+Q2)/Q1, meaningful only when Q1 is real chain truncation
    (Q3=Q4=0); an infinite periodic chain has Q1=0 and MCL is reported as
    inf rather than dividing by zero. When Q3/Q4 are present, `verdict`
    describes the connectivity instead of forcing it through that formula.
    """
    bonds = struct.neighbours()
    elems = struct.elements
    si = [i for i, e in enumerate(elems) if e == "Si"]
    if not si:
        raise ValueError("no Si atoms - not a silicate structure")

    qn = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    for i in si:
        o_neighbours = [j for j in bonds[i] if elems[j] == "O"]
        n_bridge = sum(1 for o in o_neighbours
                       if sum(1 for k in bonds[o] if elems[k] == "Si") == 2)
        qn[min(n_bridge, 4)] += 1

    comp = struct.composition()
    ca_si = comp.get("Ca", 0) / comp["Si"]
    n = len(si)
    branched = qn[3] > 0 or qn[4] > 0
    mcl = float("inf") if qn[1] == 0 else 2 * (qn[1] + qn[2]) / qn[1]

    if branched and abs(qn[3] - n / 3) < 0.6 and qn[0] == qn[1] == qn[4] == 0:
        verdict = "double chain (Q3 at the bridging site, Q2 elsewhere)"
    elif qn[2] == n and qn[0] == qn[1] == qn[3] == qn[4] == 0:
        verdict = "single chain, fully connected (all Q2)"
    elif qn[1] == n and qn[0] == qn[2] == qn[3] == qn[4] == 0:
        verdict = "isolated dimers (all Q1, bridging tetrahedra removed)"
    elif branched:
        verdict = "cross-linked / sheet-like (irregular Q3+)"
    else:
        verdict = "mixed chain lengths (Q0-Q2 present)"

    return {
        "ca_si": round(ca_si, 4),
        "n_si": n,
        "Qn": qn,
        "mean_chain_length": None if branched else mcl,
        "verdict": verdict,
    }


def load_reference(phase: str) -> Structure:
    """Read a bundled reference cell and attach its crystallographic cell."""
    if phase not in REFERENCE:
        raise ValueError(f"unknown phase {phase!r}; known: {sorted(REFERENCE)}")
    spec = REFERENCE[phase]
    path = ASSETS / spec["file"]
    if not path.exists():
        raise FileNotFoundError(f"reference structure missing: {path}")
    struct = read_cif(path) if path.suffix == ".cif" else read_xyz(path)
    if spec["cell"]:
        struct.cell = cell_from_parameters(*spec["cell"])
    struct.info.update({"phase": phase, "formula": spec["formula"],
                        "note": spec["note"], "source": spec["file"]})
    return struct


def substitute(struct: Structure, frm: str, to: str, count: int | None = None,
               fraction: float | None = None, seed: int = 0,
               label: str | None = None) -> Structure:
    """Replace some `frm` atoms with `to`, chosen at random.

    Used for cation substitution across every phase: Fe(III) for Al(III) in AFm
    and AFt, Fe(II) for Mg(II) in LDH. Random selection rather than a fixed
    pattern, because the site preference is generally what a study is trying to
    determine - imposing one would prejudge it."""
    rng = np.random.default_rng(seed)
    idx = [i for i, l in enumerate(struct.labels) if element_of(l) == frm]
    if not idx:
        raise ValueError(f"no {frm} atoms present")
    if fraction is not None:
        count = int(round(len(idx) * fraction))
    if count is None:
        raise ValueError("give either count or fraction")
    if count > len(idx):
        raise ValueError(f"asked for {count} substitutions but only {len(idx)} {frm}")
    chosen = rng.choice(idx, size=count, replace=False)
    out = struct.copy()
    for i in chosen:
        out.labels[int(i)] = label or to
    out.info = dict(struct.info)
    out.info["substitution"] = f"{count} {frm} -> {label or to}"
    return out


def set_interlayer_water(struct: Structure, n_water: int, seed: int = 0,
                         axis: int = 2) -> Structure:
    """Adjust the interlayer water count of a layered phase.

    Removing water picks the molecules furthest from the framework first, so
    the ones that remain are the more strongly bound; adding places new ones in
    the gallery midplane with the same clash test used when building."""
    rng = np.random.default_rng(seed)
    out = struct.copy()
    waters = [g for g in detect_groups(out) if g.kind == "H2O"]
    current = len(waters)

    if n_water == current:
        return out
    if n_water < current:
        framework = [i for i, l in enumerate(out.labels)
                     if element_of(l) not in ("H", "O")]
        if not framework:
            drop = waters[n_water:]
        else:
            fpos = out.positions[framework]
            def isolation(g):
                d = np.linalg.norm(out.mic(fpos - g.centre), axis=1)
                return float(d.min())
            waters.sort(key=isolation, reverse=True)
            drop = waters[: current - n_water]
        out.delete([i for g in drop for i in g.indices])
        out.info = dict(struct.info)
        out.info["n_water"] = n_water
        return out

    zs = out.positions[:, axis]
    z_mid = float(np.median(zs))
    span = float(zs.max() - zs.min()) * 0.15
    _fill_interlayer(out, "H2O", n_water - current, z_mid, max(span, 1.0), rng)
    out.info = dict(struct.info)
    out.info["n_water"] = n_water
    return out


def build_csh(ca_si: float | None = None, al_si: float = 0.0, seed: int = 0,
              supercell=(1, 1, 1), cell=None) -> Structure:
    """C-S-H at a chosen Ca/Si, from the reference series.

    `ca_si` picks the nearest bundled reference and reports the exact ratio it
    provides; omit it for the 14 A tobermorite default. Ratios between the
    references are NOT interpolated - reaching one means removing bridging
    silicate tetrahedra and charge-compensating, and there is more than one
    published convention for which sites to remove.

    `al_si` substitutes Al for Si to give C-A-S-H. Al(III) on a Si(IV) site
    leaves one negative charge per substitution, which the caller must balance;
    the returned structure reports the imbalance rather than hiding it.
    """
    if ca_si is None:
        key = 0.833
    else:
        key = min(CSH_REFERENCES, key=lambda k: abs(k - ca_si))
        if abs(key - ca_si) > 0.03:
            available = ", ".join(f"{k:.3f}" for k in sorted(CSH_REFERENCES))
            raise NotImplementedError(
                f"no C-S-H reference at Ca/Si = {ca_si:.3f}. Available: "
                f"{available}. Reaching an intermediate ratio means removing "
                "bridging silicate tetrahedra and charge-compensating, which "
                "is a modelling decision with several conventions - supply a "
                "reference structure at the target ratio instead.")

    spec = CSH_REFERENCES[key]
    path = ASSETS / spec["file"]
    if not path.exists():
        raise FileNotFoundError(f"C-S-H reference missing: {path}")
    struct = read_xyz(path)

    if cell is not None:
        struct.cell = cell_from_parameters(*cell)
    elif spec["cell"]:
        struct.cell = cell_from_parameters(*spec["cell"])
    else:
        raise ValueError(
            f"the Ca/Si = {key} reference ({spec['file']}) has no cell on "
            f"record - {spec['note']}. Pass cell=(a,b,c,alpha,beta,gamma).")

    struct.info.update({"phase": "CSH", "formula": spec["formula"],
                        "note": spec["note"], "source": spec["file"],
                        "cell_verified": spec["cell_verified"]})
    if supercell != (1, 1, 1):
        struct = struct.repeat(*supercell)

    comp = struct.composition()
    n_si, n_ca = comp.get("Si", 0), comp.get("Ca", 0)
    actual = n_ca / n_si if n_si else float("nan")

    if al_si > 0:
        n_al = int(round(n_si * al_si))
        struct = substitute(struct, "Si", "Al", count=n_al, seed=seed)
        struct.info["charge_note"] = (
            f"{n_al} Al(III) on Si(IV) sites leaves {n_al}- to balance "
            "(add Ca2+ or protonate bridging O)")

    struct.info.update({"ca_si": round(actual, 4), "al_si": al_si})
    return struct


def summary(struct: Structure, model: ChargeModel | None = None) -> str:
    q = formal_charge(struct, model or ChargeModel())
    groups: dict[str, int] = {}
    for g in detect_groups(struct):
        groups[g.kind] = groups.get(g.kind, 0) + 1
    lines = [
        f"  atoms        {len(struct)}",
        f"  composition  {struct.composition()}",
        f"  cell         {np.round(struct.lengths, 4).tolist()}  "
        f"angles {np.round(struct.angles, 2).tolist()}",
        f"  groups       {groups}",
        f"  formal charge {q['total']:+d}",
    ]
    if q["assumed"]:
        lines.append(f"  assumed states {q['assumed']} (not derived)")
    bad = close_contacts(struct)
    if bad:
        lines.append(f"  CLOSE CONTACTS {len(bad)}, worst {bad[0][2]:.2f} A")
    return "\n".join(lines)
