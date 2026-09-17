#!/usr/bin/env python3
"""
Core structure handling for cement_craft.

A deliberately small structure model: labels, Cartesian positions, and a 3x3
cell. Everything else in the toolkit is built on this plus the geometry helpers
here - periodic distances, neighbour lists, molecular-group detection and layer
clustering - because those are the operations every cement phase needs.

Only NumPy is required. ASE is used opportunistically for CIF files with
symmetry operators, and its absence is reported rather than guessed around.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Covalent radii (angstrom) for the elements that appear in cement systems.
# Used only for bond perception, so approximate values are fine.
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Ti": 1.60, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Sr": 1.95, "Zr": 1.75, "Ba": 2.15, "La": 2.07, "Ce": 2.04,
}

BOND_TOLERANCE = 1.25   # a pair bonds if d < tol * (r_i + r_j)


def element_of(label: str) -> str:
    """Strip a CP2K kind label (Fe_A, Ow, Mg1) down to its element symbol."""
    if not label:
        return label
    base = label.split("_")[0]
    base = re.sub(r"\d+$", "", base)
    if base in COVALENT_RADII:
        return base
    if len(base) > 1 and base[:1] in COVALENT_RADII:
        return base[:1]
    return base


@dataclass
class Structure:
    """Labels, positions (N,3) in angstrom, and a (3,3) cell matrix."""

    labels: list[str]
    positions: np.ndarray
    cell: np.ndarray
    info: dict = field(default_factory=dict)

    def __post_init__(self):
        self.positions = np.asarray(self.positions, dtype=float).reshape(-1, 3)
        self.cell = np.asarray(self.cell, dtype=float).reshape(3, 3)
        if len(self.labels) != len(self.positions):
            raise ValueError(
                f"{len(self.labels)} labels but {len(self.positions)} positions"
            )

    # -- basics -----------------------------------------------------------
    def __len__(self):
        return len(self.labels)

    @property
    def elements(self) -> list[str]:
        return [element_of(l) for l in self.labels]

    @property
    def lengths(self) -> np.ndarray:
        return np.linalg.norm(self.cell, axis=1)

    @property
    def angles(self) -> np.ndarray:
        a, b, c = self.cell
        def ang(u, v):
            return math.degrees(
                math.acos(
                    float(np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)),
                                  -1.0, 1.0))
                )
            )
        return np.array([ang(b, c), ang(a, c), ang(a, b)])

    @property
    def volume(self) -> float:
        return abs(float(np.linalg.det(self.cell)))

    def composition(self) -> dict:
        out: dict = {}
        for e in self.elements:
            out[e] = out.get(e, 0) + 1
        return dict(sorted(out.items()))

    def label_counts(self) -> dict:
        out: dict = {}
        for l in self.labels:
            out[l] = out.get(l, 0) + 1
        return dict(sorted(out.items()))

    def copy(self) -> "Structure":
        return Structure(list(self.labels), self.positions.copy(), self.cell.copy(),
                         dict(self.info))

    # -- geometry ---------------------------------------------------------
    def fractional(self) -> np.ndarray:
        return self.positions @ np.linalg.inv(self.cell)

    def wrap(self) -> "Structure":
        """Fold every atom into the primary cell."""
        frac = self.fractional() % 1.0
        self.positions = frac @ self.cell
        return self

    def wrap_molecules(self, groups=None) -> "Structure":
        """Fold the structure into the primary cell one molecule at a time.

        A plain per-atom wrap() (or the raw coordinates from a builder that
        never wraps at all) is correct for any PBC-aware calculation, which
        only ever sees minimum-image distances - but it can leave one H of a
        water molecule on the far side of the cell from its own O, which
        renders as a bond stretched across the whole box in VMD, VESTA, or any
        viewer that plots raw coordinates. The physics is identical either
        way; only the picture changes.

        Every atom in a FINITE, non-overlapping group (H2O, OH, an isolated
        CO3 or SO4) is shifted by the same whole number of lattice vectors,
        chosen so the group's first atom lands in the primary cell - a rigid
        shift, so bond lengths and internal geometry are exactly preserved.
        A molecule that straddles a face can still poke slightly outside the
        cell after this by about a bond length - that is normal and is what
        every MD viewer's "wrap by molecule" does too.

        A polymerised species - SiO4 tetrahedra sharing bridging oxygens along
        a silicate chain - is left alone entirely: detect_groups() lists each
        tetrahedron separately, so a bridging oxygen is a member of two groups
        at once, and a chain spanning the whole cell is not a finite molecule
        to begin with (it genuinely continues into the next periodic image).
        Assembling it is not what this method is for - it gets the same plain
        per-atom wrap as any other framework atom, split-looking seam and
        all, which is the correct picture for something that really is
        periodic."""
        groups = groups if groups is not None else detect_groups(self)
        inv = np.linalg.inv(self.cell)

        membership: dict[int, int] = {}
        for g in groups:
            for i in g.indices:
                membership[i] = membership.get(i, 0) + 1

        grouped = set()
        for g in groups:
            if any(membership[i] > 1 for i in g.indices):
                continue                  # shares an atom with another group
            anchor = g.indices[0]
            # .copy() is essential: self.positions[anchor] is a VIEW into the
            # array. Without it, reassigning self.positions[anchor] below (for
            # i == anchor) silently mutates anchor_pos too, since they alias
            # the same memory - and every later member in this group would
            # then have the wrap shift applied to it a second time.
            anchor_pos = self.positions[anchor].copy()
            anchor_frac = anchor_pos @ inv
            shift = np.floor(anchor_frac) @ self.cell
            for i in g.indices:
                if i == anchor:
                    self.positions[i] = anchor_pos - shift
                else:
                    delta = self.mic(self.positions[i] - anchor_pos)[0]
                    self.positions[i] = (anchor_pos - shift) + delta
                grouped.add(i)

        # Everything else - framework atoms, and any polymerised-species
        # member left unassembled above - gets an ordinary independent wrap.
        # Deliberately no attempt to keep a silicate chain looking bonded
        # across the cell: it is not a finite molecule, and reconstructing
        # its connectivity is a different, larger job than the one this
        # method is for.
        ungrouped = [i for i in range(len(self)) if i not in grouped]
        if ungrouped:
            frac = self.positions[ungrouped] @ inv
            self.positions[ungrouped] = (frac % 1.0) @ self.cell
        return self

    def mic(self, delta: np.ndarray) -> np.ndarray:
        """Minimum-image convention for displacement vectors."""
        delta = np.atleast_2d(delta)
        frac = delta @ np.linalg.inv(self.cell)
        frac -= np.round(frac)
        return frac @ self.cell

    def distance(self, i: int, j: int) -> float:
        return float(np.linalg.norm(self.mic(self.positions[j] - self.positions[i])[0]))

    def distances_from(self, i: int) -> np.ndarray:
        return np.linalg.norm(self.mic(self.positions - self.positions[i]), axis=1)

    def neighbours(self, tolerance: float = BOND_TOLERANCE) -> list[list[int]]:
        """Bond list by covalent-radius overlap, under the minimum image.

        O(N^2) but vectorised per atom, which is fine for the few-thousand-atom
        cells these phases produce."""
        elems = self.elements
        radii = np.array([COVALENT_RADII.get(e, 1.2) for e in elems])
        out: list[list[int]] = [[] for _ in range(len(self))]
        for i in range(len(self)):
            d = self.distances_from(i)
            cut = tolerance * (radii + radii[i])
            hits = np.where((d < cut) & (d > 1e-6))[0]
            out[i] = [int(j) for j in hits]
        return out

    def layers(self, axis: int = 2, tolerance: float = 1.0,
               indices=None) -> list[list[int]]:
        """Cluster atoms into layers along an axis, ordered low to high.

        Used for slab construction and for locating metal sheets in layered
        phases, where the meaningful structural unit is a plane, not a plane
        index in some idealised lattice."""
        idx = list(range(len(self))) if indices is None else list(indices)
        if not idx:
            return []
        idx.sort(key=lambda i: self.positions[i, axis])
        groups = [[idx[0]]]
        ref = self.positions[idx[0], axis]
        for i in idx[1:]:
            z = self.positions[i, axis]
            if abs(z - ref) <= tolerance:
                groups[-1].append(i)
                ref = float(np.mean([self.positions[k, axis] for k in groups[-1]]))
            else:
                groups.append([i])
                ref = z
        return groups

    def extent(self, axis: int) -> float:
        col = self.positions[:, axis]
        return float(col.max() - col.min())

    # -- editing ----------------------------------------------------------
    def extend(self, labels, positions) -> "Structure":
        self.labels.extend(list(labels))
        self.positions = np.vstack([self.positions, np.atleast_2d(positions)])
        return self

    def delete(self, indices) -> "Structure":
        keep = [i for i in range(len(self)) if i not in set(indices)]
        self.labels = [self.labels[i] for i in keep]
        self.positions = self.positions[keep]
        return self

    def repeat(self, na: int, nb: int, nc: int) -> "Structure":
        labels, pos = [], []
        for i in range(na):
            for j in range(nb):
                for k in range(nc):
                    shift = i * self.cell[0] + j * self.cell[1] + k * self.cell[2]
                    labels.extend(self.labels)
                    pos.append(self.positions + shift)
        cell = self.cell * np.array([[na], [nb], [nc]])
        return Structure(labels, np.vstack(pos), cell, dict(self.info))

    def sort_by(self, axis: int = 2) -> "Structure":
        """Sort atoms along an axis.

        Index-based constraints (CP2K &FIXED_ATOMS) refer to file order, so
        sorting before writing - never after - is what keeps a fixed-layer
        selection meaningful."""
        order = np.argsort(self.positions[:, axis], kind="stable")
        self.labels = [self.labels[i] for i in order]
        self.positions = self.positions[order]
        return self


# --------------------------------------------------------------------------
# Molecular group detection
# --------------------------------------------------------------------------

@dataclass
class Group:
    kind: str                  # 'H2O', 'OH', 'CO3', 'SO4', 'SiO4'
    indices: list[int]
    centre: np.ndarray


def detect_groups(struct: Structure, bonds=None) -> list[Group]:
    """Find water, hydroxyl, carbonate, sulfate and silicate units.

    Perception is by connectivity rather than by label, so a structure that
    arrived from someone else's naming convention still resolves correctly.
    An oxygen with one H is hydroxyl, with two is water; oxygens bonded to a
    central C/S/Si belong to that oxyanion."""
    bonds = bonds if bonds is not None else struct.neighbours()
    elems = struct.elements
    groups: list[Group] = []
    claimed: set[int] = set()

    # Oxyanions first, so their oxygens are not mistaken for hydroxyl.
    for centre_el, kind, n_o in (("S", "SO4", 4), ("C", "CO3", 3), ("Si", "SiO4", 4)):
        for i, e in enumerate(elems):
            if e != centre_el or i in claimed:
                continue
            oxy = [j for j in bonds[i] if elems[j] == "O"]
            if len(oxy) < n_o - 1:
                continue
            members = [i] + oxy
            # Protonated oxyanion oxygens keep their H with the group.
            for j in list(oxy):
                members += [k for k in bonds[j] if elems[k] == "H"]
            groups.append(Group(kind, sorted(set(members)),
                                _centre(struct, members)))
            claimed.update(members)

    for i, e in enumerate(elems):
        if e != "O" or i in claimed:
            continue
        hs = [j for j in bonds[i] if elems[j] == "H"]
        if len(hs) >= 2:
            members = [i] + hs[:2]
            groups.append(Group("H2O", sorted(members), _centre(struct, members)))
            claimed.update(members)
        elif len(hs) == 1:
            members = [i, hs[0]]
            groups.append(Group("OH", sorted(members), _centre(struct, members)))
            claimed.update(members)

    return groups


def _centre(struct: Structure, members) -> np.ndarray:
    """Centroid under the minimum image, anchored on the first member so a
    group straddling a periodic boundary does not average to the cell middle."""
    anchor = struct.positions[members[0]]
    rel = struct.mic(struct.positions[list(members)] - anchor)
    return anchor + rel.mean(axis=0)


# --------------------------------------------------------------------------
# File I/O
# --------------------------------------------------------------------------

def read_xyz(path) -> Structure:
    """Read xyz or extended-xyz. A Lattice="..." comment is used when present."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    n = int(lines[0].split()[0])
    comment = lines[1] if len(lines) > 1 else ""
    cell = np.zeros((3, 3))
    m = re.search(r'Lattice="([^"]+)"', comment)
    if m:
        v = [float(x) for x in m.group(1).split()]
        if len(v) == 9:
            cell = np.array(v).reshape(3, 3)
    else:
        # `CELL lx ly lz | notes` - an orthorhombic shorthand used by some
        # interface-building workflows. Anything after the pipe is commentary.
        m = re.match(r"\s*CELL\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)",
                     comment)
        if m:
            cell = np.diag([float(m.group(k)) for k in (1, 2, 3)])
    labels, pos = [], []
    for line in lines[2 : 2 + n]:
        p = line.split()
        if len(p) < 4:
            break
        labels.append(p[0])
        pos.append([float(x) for x in p[1:4]])
    if len(labels) != n:
        raise ValueError(f"{path}: header says {n} atoms, found {len(labels)}")
    return Structure(labels, np.array(pos), cell, {"comment": comment})


def write_xyz(path, struct: Structure, comment: str | None = None,
              lattice: bool = True, wrap_molecules: bool = False) -> None:
    """Write extended xyz. The Lattice tag keeps the cell with the structure,
    which matters because a cement slab is meaningless without it.

    `wrap_molecules=True` assembles each water/hydroxyl/oxyanion into a single
    contiguous unit before writing (see Structure.wrap_molecules) - purely
    cosmetic, so a copy is wrapped rather than the caller's own structure."""
    if wrap_molecules:
        struct = struct.copy().wrap_molecules()
    if comment is None:
        if lattice and struct.cell.any():
            flat = " ".join(f"{x:.10f}" for x in struct.cell.reshape(-1))
            comment = f'Lattice="{flat}" Properties=species:S:1:pos:R:3 pbc="T T T"'
        else:
            comment = struct.info.get("comment", "")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{len(struct)}\n{comment}\n")
        for lab, (x, y, z) in zip(struct.labels, struct.positions):
            fh.write(f"{lab:<5}{x:20.10f}{y:20.10f}{z:20.10f}\n")


def cell_from_parameters(a, b, c, alpha=90.0, beta=90.0, gamma=90.0) -> np.ndarray:
    """Standard crystallographic cell matrix: a along x, b in the xy plane."""
    al, be, ga = (math.radians(x) for x in (alpha, beta, gamma))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(ga), b * math.sin(ga), 0.0])
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz2 = c * c - cx * cx - cy * cy
    vc = np.array([cx, cy, math.sqrt(max(cz2, 0.0))])
    return np.vstack([va, vb, vc])


class OccupancyError(ValueError):
    """Raised when a CIF's partial occupancies cannot be used as they stand."""


def _supercell_for_occupancies(occs, max_n: int = 8):
    """Smallest supercell multiplier making every occupancy a whole number.

    An occupancy of 1/3 needs three of that site before it can be represented
    by whole atoms; 0.25 needs four. The least common multiple over all partial
    sites is the smallest cell in which an ordered approximant exists at all."""
    from fractions import Fraction
    lcm = 1
    for o in occs:
        if abs(o - round(o)) < 1e-6:
            continue
        d = Fraction(o).limit_denominator(12).denominator
        lcm = lcm * d // math.gcd(lcm, d)
        if lcm > max_n ** 3:
            return None
    return lcm


def read_cif(path, occupancy: str = "strict", threshold: float = 0.5,
             seed: int = 0) -> Structure:
    """Read a CIF, expanding symmetry operators when present.

    **Partial occupancies are not silently ignored.** A crystallographic site
    at occupancy 0.5 means the site is filled half the time across the crystal,
    which is a statistical statement; an atomistic simulation needs whole atoms
    in definite places, and a disordered CIF does not define a unique ordered
    structure. Different valid orderings give different energies, so the choice
    has to be made deliberately.

    `occupancy` selects what to do:

    * ``"strict"`` (default) - refuse, reporting which sites are partial and
      the smallest supercell in which they could be made integral.
    * ``"round"`` - keep sites with occupancy >= `threshold`, drop the rest.
      Fast, but it changes stoichiometry and usually breaks charge balance;
      check the composition afterwards.
    * ``"ignore"`` - take every site at face value. Only correct when the
      partial sites are split positions for the *same* atom, and even then it
      duplicates atoms.

    For a proper ordered approximant, build the supercell this reports and
    enumerate orderings, comparing several - one configuration is a sample of
    the disorder, not the answer.

    Falls back to ASE when installed, since its CIF handling covers more of the
    format than is worth reimplementing; the built-in path handles the common
    case of explicit symmetry operation lists.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")

    def _num(key):
        m = re.search(rf"^{key}\s+([-\d.]+)", text, re.MULTILINE)
        return float(m.group(1)) if m else None

    a, b, c = (_num(f"_cell_length_{x}") for x in "abc")
    al, be, ga = (_num(f"_cell_angle_{x}") for x in ("alpha", "beta", "gamma"))
    if None in (a, b, c):
        try:
            from ase.io import read as ase_read  # noqa: F401
        except ImportError:
            raise ValueError(
                f"{path}: could not read cell parameters, and ASE is not "
                "installed to fall back on (pip install ase)"
            )
        from ase.io import read as ase_read
        atoms = ase_read(str(path))
        return Structure(list(atoms.get_chemical_symbols()),
                         atoms.get_positions(), np.array(atoms.get_cell()))
    cell = cell_from_parameters(a, b, c, al or 90, be or 90, ga or 90)

    # Symmetry operations. They may be quoted or bare, and may sit in a loop
    # with a leading index column - all three forms occur in the wild. Reading
    # only the quoted form silently drops every operator but the identity,
    # which yields a structure with a fraction of the atoms and a density to
    # match, so parse the loop properly and fall back to quotes only after.
    ops: list[str] = []
    sym_header, sym_rows = _cif_loop(text, "_symmetry_equiv_pos_as_xyz")
    if not sym_header:
        sym_header, sym_rows = _cif_loop(text, "_space_group_symop_operation_xyz")
    if sym_header:
        for row in sym_rows:
            for token in row:
                tok = token.strip().strip("'\"")
                if tok.count(",") == 2 and re.fullmatch(r"[-+0-9xyzXYZ/,. ]+", tok):
                    ops.append(tok)
                    break
    if not ops:
        ops = [o for o in re.findall(r"['\"]([-+0-9xyzXYZ/,. ]+)['\"]", text)
               if o.count(",") == 2]
    if not ops:
        ops = ["x,y,z"]

    # atom site loop
    header, rows = _cif_loop(text, "_atom_site_")
    if not header:
        raise ValueError(f"{path}: no _atom_site_ loop found")
    def col(*names):
        for nm in names:
            if nm in header:
                return header.index(nm)
        return None
    ci = col("_atom_site_type_symbol", "_atom_site_label")
    cx = col("_atom_site_fract_x")
    cy = col("_atom_site_fract_y")
    cz = col("_atom_site_fract_z")
    co = col("_atom_site_occupancy")
    ch = col("_atom_site_attached_hydrogens")
    if None in (ci, cx, cy, cz):
        raise ValueError(f"{path}: atom site loop lacks fractional coordinates")

    labels, frac, occs, attached = [], [], [], []
    for row in rows:
        if len(row) <= max(ci, cx, cy, cz):
            continue
        sym = re.sub(r"[^A-Za-z]", "", row[ci])
        sym = element_of(sym) if element_of(sym) in COVALENT_RADII else sym[:2]
        try:
            f = [float(re.sub(r"\(.*\)", "", row[k])) for k in (cx, cy, cz)]
        except ValueError:
            continue
        o = 1.0
        if co is not None and len(row) > co:
            try:
                o = float(re.sub(r"\(.*\)", "", row[co]))
            except ValueError:
                o = 1.0
        # Refinements that cannot locate hydrogen often still record how many
        # belong on each site. That is the crystallographer stating where the
        # protons go, and it is far better than inferring it later.
        nh = 0
        if ch is not None and len(row) > ch:
            try:
                nh = int(float(row[ch]))
            except ValueError:
                nh = 0
        labels.append(sym)
        frac.append(f)
        occs.append(o)
        attached.append(nh)

    partial = [(l, f, o) for l, f, o in zip(labels, frac, occs) if o < 0.999]
    if partial:
        needed = _supercell_for_occupancies([o for _, _, o in partial])
        report = ", ".join(
            f"{l} @ {o:.3f}" for l, _, o in sorted(partial, key=lambda t: -t[2])[:8]
        )
        if occupancy == "strict":
            if needed and needed <= 8:
                hint = (f"a {needed}x supercell would make them integral, so an "
                        "ordered approximant exists at reasonable size")
            else:
                hint = (f"the smallest integral supercell is {needed}x"
                        if needed else "no supercell under 8^3 makes them integral")
                hint += (", i.e. these occupancies are not simple fractions and no "
                         "practical ordered approximant reproduces them exactly")
            raise OccupancyError(
                f"{Path(path).name}: {len(partial)} of {len(labels)} sites have "
                f"partial occupancy ({report}). A disordered CIF does not define "
                f"a unique ordered structure - {hint}. Choose explicitly: "
                "read_cif(..., occupancy='round') to keep sites above a "
                "threshold, or build the supercell and enumerate orderings, "
                "comparing several rather than trusting one.")
        if occupancy == "round":
            keep = [i for i, o in enumerate(occs) if o >= threshold]
            dropped = len(labels) - len(keep)
            labels = [labels[i] for i in keep]
            frac = [frac[i] for i in keep]
            occs = [occs[i] for i in keep]
            _occ_note = (f"rounded at {threshold}: dropped {dropped} partial "
                         f"site(s); stoichiometry and charge will have changed")
        elif occupancy == "ignore":
            _occ_note = (f"{len(partial)} partial site(s) taken at face value; "
                         "atoms may be duplicated at split positions")
        else:
            raise ValueError(f"unknown occupancy mode {occupancy!r}")
    else:
        _occ_note = None

    out_lab, out_frac, out_nh = [], [], []
    for lab, f, nh in zip(labels, frac, attached):
        for op in ops:
            g = _apply_symop(op, f)
            g = [x % 1.0 for x in g]
            if not any(
                lab == ol and np.allclose(_wrap_delta(np.array(g) - np.array(of)),
                                          0, atol=1e-4)
                for ol, of in zip(out_lab, out_frac)
            ):
                out_lab.append(lab)
                out_frac.append(g)
                out_nh.append(nh)

    positions = np.array(out_frac) @ cell
    info = {"source": str(path)}
    if _occ_note:
        info["occupancy_note"] = _occ_note
    if any(out_nh):
        info["attached_hydrogens"] = {i: n for i, n in enumerate(out_nh) if n}
    return Structure(out_lab, positions, cell, info)


def _wrap_delta(d):
    return d - np.round(d)


def _cif_loop(text, prefix):
    """Return (header, rows) for the first loop_ containing prefix columns."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        header, j = [], i + 1
        while j < len(lines) and lines[j].strip().startswith("_"):
            header.append(lines[j].strip())
            j += 1
        if not any(h.startswith(prefix) for h in header):
            continue
        rows = []
        while j < len(lines):
            s = lines[j].strip()
            if not s or s.startswith("_") or s.startswith("loop_") or s.startswith("#"):
                break
            rows.append(s.split())
            j += 1
        return header, rows
    return [], []


def _apply_symop(op, frac):
    x, y, z = frac
    out = []
    for part in op.lower().split(","):
        part = part.strip().replace(" ", "")
        expr = re.sub(r"(?<=[\d])(?=[xyz])", "*", part)
        expr = re.sub(r"(\d)/(\d)", r"(\1/\2)", expr)
        out.append(eval(expr, {"__builtins__": {}}, {"x": x, "y": y, "z": z}))
    return out


def close_contacts(struct: Structure, cutoff: float = 0.75) -> list[tuple]:
    """Atom pairs closer than cutoff - always a structure error, never an SCF one."""
    bad = []
    for i in range(len(struct)):
        d = struct.distances_from(i)
        for j in np.where((d < cutoff) & (d > 1e-9))[0]:
            if j > i:
                bad.append((i, int(j), float(d[j])))
    return bad


# --------------------------------------------------------------------------
# Completing water molecules
# --------------------------------------------------------------------------

WATER_OH_LENGTH = 0.96
WATER_HOH_ANGLE = 104.5

# Labels that conventionally mark a water oxygen in a refinement.
WATER_OXYGEN_LABELS = {"Wa", "WA", "Ow", "OW", "OW1", "OW2", "W", "Owat"}

# Cations that bind water as an ionic ligand through the oxygen lone pair,
# as opposed to a hydrogen-bond acceptor, which is approached by a hydrogen
# instead. Si is deliberately excluded - it never binds water directly, only
# through covalent Si-O-H chemistry, which is a different site entirely.
WATER_COORDINATING_CATIONS = {
    "Na", "K", "Mg", "Ca", "Al", "Sc", "Ti", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Sr", "Zr", "Ba", "La",
}


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues' rotation matrix for `angle` radians about unit `axis`."""
    x, y, z = axis
    k = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * (k @ k)


def _align_lone_pair(target: np.ndarray) -> np.ndarray:
    """A rotation matrix that sends the template's lone-pair axis (0,0,-1) to
    unit vector `target`, leaving free rotation about `target` unspecified -
    the caller composes that spin on top."""
    z = np.array([0.0, 0.0, -1.0])
    c = float(np.dot(z, target))
    if c > 1 - 1e-8:
        return np.eye(3)
    if c < -1 + 1e-8:
        perp = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(z, perp)
        return _axis_rotation(axis / np.linalg.norm(axis), math.pi)
    axis = np.cross(z, target)
    return _axis_rotation(axis / np.linalg.norm(axis), math.acos(np.clip(c, -1, 1)))


def add_water_hydrogens(struct: Structure, labels=None, hbond_range=(2.5, 3.4),
                        n_trials: int = 400, seed: int = 0,
                        relabel: str = "O", verbose: bool = False) -> dict:
    """Add the missing hydrogens to bare water oxygens.

    Crystallographic refinements against X-ray data rarely locate hydrogen -
    it scatters too weakly - so interlayer and channel water often arrives as a
    lone oxygen on a site labelled ``Wa`` or ``Ow``. That is fine for a
    diffraction model and useless for a simulation, where the missing protons
    carry charge and make the hydrogen bonds.

    Each bare oxygen gets two hydrogens at proper water geometry (0.96 A,
    104.5 degrees). Orientation is chosen by trial: the rigid H-O-H unit is
    rotated at random and scored on how well its hydrogens point at nearby
    acceptor oxygens without clashing into anything. Water in these phases is
    held in place by its hydrogen bonds, so pointing the protons at acceptors
    puts them far closer to their relaxed positions than an arbitrary
    orientation would.

    When the oxygen is also bonded to a cation (Ca, Mg, ...) - an aqua ligand,
    not just a hydrogen-bonded water - the search is constrained so the
    oxygen's lone pair points at that cation, since that is the bond the
    ligand is actually making; only the free spin about that direction is
    then chosen by the same acceptor/clash scoring.

    `relabel` renames the oxygen once it is a real water (the ``Wa`` label
    means nothing to a simulation code); pass None to keep it.
    """
    rng = np.random.default_rng(seed)
    labels_wanted = set(labels) if labels else WATER_OXYGEN_LABELS
    bonds = struct.neighbours()
    elems = struct.elements

    declared = struct.info.get("attached_hydrogens") or {}

    targets = []
    for i, lab in enumerate(struct.labels):
        # A CIF that states _atom_site_attached_hydrogens has already told us
        # which sites carry protons; trust it over any label convention.
        if declared:
            want = declared.get(i, 0)
            if want:
                n_h_now = sum(1 for j in bonds[i] if elems[j] == "H")
                if n_h_now < want:
                    targets.append((i, want - n_h_now))
            continue
        is_marked = lab in labels_wanted
        if not (is_marked or elems[i] == "O"):
            continue
        n_h = sum(1 for j in bonds[i] if elems[j] == "H")
        if n_h >= 2:
            continue
        # An unmarked oxygen is only treated as water when nothing else claims
        # it: no metal, no silicon or sulfur, no hydrogen at all.
        if not is_marked:
            continue
        if n_h == 1:
            targets.append((i, 1))
        else:
            targets.append((i, 2))

    if not targets:
        if verbose:
            print("  no bare water oxygens found")
        return {"waters_completed": 0, "hydrogens_added": 0}

    half = math.radians(WATER_HOH_ANGLE / 2)
    template = np.array([
        [WATER_OH_LENGTH * math.sin(half), 0.0, WATER_OH_LENGTH * math.cos(half)],
        [-WATER_OH_LENGTH * math.sin(half), 0.0, WATER_OH_LENGTH * math.cos(half)],
    ])

    lo, hi = hbond_range
    added = 0
    for o_index, n_needed in targets:
        # Recomputed each pass: the structure grows as hydrogens are added.
        elems = struct.elements
        origin = struct.positions[o_index]
        d = struct.distances_from(o_index)
        # Acceptors: oxygens at hydrogen-bonding distance.
        acc = [j for j in np.where((d > lo) & (d < hi))[0]
               if elems[j] == "O" and j != o_index]
        acc_dirs = [struct.mic(struct.positions[j] - origin)[0] for j in acc]
        acc_dirs = [v / np.linalg.norm(v) for v in acc_dirs if np.linalg.norm(v) > 1e-6]

        # An aqua ligand: the nearest bonded cation is who the lone pair is
        # actually pointing at, not just another acceptor to weight amongst
        # others. bonds[] was built before any H was added, but cation
        # positions never move, so it is still valid here.
        cation_j = [j for j in bonds[o_index] if elems[j] in WATER_COORDINATING_CATIONS]
        lone_pair_rot, cation_unit = None, None
        if cation_j:
            nearest = min(cation_j, key=lambda j: d[j])
            cation_unit = struct.mic(struct.positions[nearest] - origin)[0]
            cation_unit = cation_unit / np.linalg.norm(cation_unit)
            lone_pair_rot = _align_lone_pair(cation_unit)

        # Anything nearby that a hydrogen must not run into.
        near = [j for j in np.where(d < 3.2)[0] if j != o_index]
        near_pos = struct.positions[near] if near else np.zeros((0, 3))
        near_h = np.array([elems[j] == "H" for j in near])

        best, best_score = None, -1e9
        for _ in range(n_trials):
            if lone_pair_rot is not None:
                # The Ca-O direction is fixed; only the spin about it is free.
                theta = rng.uniform(0, 2 * math.pi)
                spin = _axis_rotation(cation_unit, theta)
                rot = spin @ lone_pair_rot
            else:
                q, r = np.linalg.qr(rng.normal(size=(3, 3)))
                rot = q * np.sign(np.diag(r))
            trial = origin + template @ rot.T
            if n_needed == 1:
                trial = trial[:1]

            score = 0.0
            ok = True
            for h in trial:
                rel = struct.mic(near_pos - h) if len(near) else np.zeros((0, 3))
                dist = np.linalg.norm(rel, axis=1) if len(near) else np.array([])
                if dist.size:
                    limit = np.where(near_h, 1.45, 1.55)
                    if np.any(dist < limit):
                        ok = False
                        break
                    score -= float(np.sum(np.clip(2.0 - dist, 0, None)))
                # Reward alignment with an acceptor direction.
                hdir = (h - origin) / WATER_OH_LENGTH
                if acc_dirs:
                    score += 2.0 * max(float(np.dot(hdir, a)) for a in acc_dirs)
            if ok and score > best_score:
                best, best_score = trial, score

        if best is None:
            raise RuntimeError(
                f"could not place hydrogens on water oxygen {o_index} without "
                "a clash - the site is too crowded; check the structure")
        struct.extend(["H"] * len(best), best)
        added += len(best)
        if relabel:
            struct.labels[o_index] = relabel

    if verbose:
        print(f"  completed {len(targets)} water molecules, added {added} H")
    return {"waters_completed": len(targets), "hydrogens_added": added}


def transform_cell(struct: Structure, matrix, align: bool = True,
                   tol: float = 1e-4) -> Structure:
    """Re-express a structure on a new lattice basis defined by an integer matrix.

    `matrix` gives the new cell vectors as integer combinations of the old, so
    its determinant is the number of primitive cells the new cell contains.
    The classic use is turning a monoclinic setting into an orthogonal
    supercell: a layered silicate whose chains run along a non-orthogonal `b`
    comes out with those chains diagonal in Cartesian space, which is awkward
    to constrain, to slice, and to look at.

    With `align`, the result is additionally rotated so the cell vectors lie on
    the Cartesian axes where the new cell is orthogonal. That is a rigid
    rotation - no coordinate is distorted - and it is what puts a chain
    direction along x, y or z.
    """
    matrix = np.asarray(matrix, dtype=int).reshape(3, 3)
    det = int(round(abs(np.linalg.det(matrix))))
    if det == 0:
        raise ValueError("transformation matrix is singular")

    new_cell = matrix @ struct.cell
    inv_new = np.linalg.inv(new_cell)

    # Search far enough in the old basis to cover the new cell in every
    # direction; the bound is loose on purpose, duplicates are filtered below.
    reach = int(abs(matrix).sum()) + 1
    labels: list[str] = []
    positions: list[np.ndarray] = []
    seen: list[np.ndarray] = []

    for i in range(-reach, reach + 1):
        for j in range(-reach, reach + 1):
            for k in range(-reach, reach + 1):
                shift = i * struct.cell[0] + j * struct.cell[1] + k * struct.cell[2]
                trial = struct.positions + shift
                frac = trial @ inv_new
                inside = np.all((frac >= -tol) & (frac < 1.0 - tol), axis=1)
                for idx in np.where(inside)[0]:
                    f = frac[idx] % 1.0
                    if any(
                        struct.labels[idx] == lab
                        and np.allclose(((f - g + 0.5) % 1.0) - 0.5, 0, atol=1e-3)
                        for lab, g in zip(labels, seen)
                    ):
                        continue
                    labels.append(struct.labels[idx])
                    seen.append(f)
                    positions.append(trial[idx])

    expected = len(struct) * det
    if len(labels) != expected:
        raise RuntimeError(
            f"transformation produced {len(labels)} atoms, expected "
            f"{expected} ({len(struct)} x |det| {det}). The matrix may not map "
            "the lattice onto itself.")

    out = Structure(labels, np.array(positions), new_cell, dict(struct.info))

    if align:
        lengths = np.linalg.norm(new_cell, axis=1)
        gram = new_cell @ new_cell.T
        off = max(abs(gram[0, 1]), abs(gram[0, 2]), abs(gram[1, 2]))
        if off / (lengths[0] * lengths[1]) > 1e-3:
            raise ValueError(
                "align=True requires an orthogonal new cell; this one is not "
                f"(largest off-diagonal {off:.4f}). Use align=False.")
        # Fractional coordinates are invariant under the rotation, so rebuilding
        # positions from them against a diagonal cell IS the rotation.
        frac = out.positions @ np.linalg.inv(new_cell)
        out.cell = np.diag(lengths)
        out.positions = (frac % 1.0) @ out.cell

    out.info["transformed_by"] = matrix.tolist()
    return out
