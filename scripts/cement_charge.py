#!/usr/bin/env python3
"""
Formal charge bookkeeping for cement phases.

Charge balance is what decides whether a constructed cement structure is
physical. A layered double hydroxide sheet carries a positive charge set by its
trivalent fraction, and that charge must be met by interlayer anions; a cut
surface leaves dangling oxygens that must be protonated back to neutrality.
Getting it wrong does not crash anything downstream - it produces a charged cell
that a DFT code will happily converge.

Main-group oxidation states are treated as known. Transition metals are not:
their state is asked for, or defaulted to the most common one and reported as
an assumption, never folded silently into a total.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from settings import MAIN_GROUP

from cement_core import Structure, detect_groups, element_of  # noqa: E402

# Main-group and alkaline-earth states relevant to cement chemistry. These are
# unambiguous in this context and are applied without asking.
MAIN_GROUP = {
    "H": 1, "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
    "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2,
    "B": 3, "Al": 3, "Ga": 3, "In": 3,
    "C": 4, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 2,
    "N": -3, "P": 5, "As": 5, "Sb": 3, "Bi": 3,
    "O": -2, "S": 6, "Se": 4, "Te": 4,
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
}

# Most common state, then other states seen often enough to matter. Applied
# only as a labelled default; the caller is expected to override when the
# chemistry says otherwise.
TRANSITION_METAL = {
    "Sc": (3, []), "Ti": (4, [3]), "V": (3, [4, 5]), "Cr": (3, [6]),
    "Mn": (2, [3, 4]), "Fe": (3, [2]), "Co": (2, [3]), "Ni": (2, [3]),
    "Cu": (2, [1]), "Zn": (2, []),
    "Y": (3, []), "Zr": (4, []), "Nb": (5, []), "Mo": (6, [4]),
    "Ag": (1, []), "Cd": (2, []),
    "La": (3, []), "Ce": (3, [4]), "Hf": (4, []), "Ta": (5, []), "W": (6, [4]),
}

# Oxyanion charges, used when a group is perceived rather than summed atomwise.
#
# Only ISOLATED oxyanions belong here. Silicate is deliberately absent: in
# C-S-H the tetrahedra are polymerised into chains, so their bridging oxygens
# are shared between two silicons. Charging each SiO4 as a discrete -4 unit
# counts those oxygens twice and puts a tobermorite cell around -100 when it is
# in fact neutral. Silicate therefore falls through to the atomwise sum, where
# every oxygen is counted exactly once whether it bridges or not.
GROUP_CHARGE = {"H2O": 0, "OH": -1, "CO3": -2, "SO4": -2}


class ChargeModel:
    """Assigns an oxidation state to every element in a structure.

    `overrides` wins over everything; transition metals absent from it fall back
    to the most common state and are recorded in `assumed` so the caller can
    report them."""

    def __init__(self, overrides: dict | None = None):
        self.overrides = {k: int(v) for k, v in (overrides or {}).items()}
        self.assumed: dict[str, int] = {}
        self.unknown: set[str] = set()

    def state(self, element: str) -> int:
        if element in self.overrides:
            return self.overrides[element]
        if element in MAIN_GROUP:
            return MAIN_GROUP[element]
        if element in TRANSITION_METAL:
            primary, _alts = TRANSITION_METAL[element]
            self.assumed[element] = primary
            return primary
        self.unknown.add(element)
        return 0

    def alternatives(self, element: str) -> list[int]:
        return TRANSITION_METAL.get(element, (None, []))[1]

    def needs_input(self, struct: Structure) -> list[str]:
        """Transition metals present with no explicit override.

        The caller should ask about these rather than accept the default
        silently - which metal state is right changes the whole charge balance
        of a substituted phase."""
        present = {element_of(l) for l in struct.labels}
        return sorted(
            e for e in present
            if e in TRANSITION_METAL and e not in self.overrides
        )


def formal_charge(struct: Structure, model: ChargeModel | None = None,
                  use_groups: bool = True) -> dict:
    """Total formal charge, by group where groups are perceived.

    Summing atomwise states double-counts nothing but mis-handles peroxo-like
    bonding; for cement phases the group route (OH-, CO3 2-, SO4 2-, SiO4 4-)
    is both more robust and closer to how the chemistry is actually reasoned
    about."""
    model = model or ChargeModel()
    elems = struct.elements
    total = 0
    breakdown: dict[str, int] = {}
    claimed: set[int] = set()

    if use_groups:
        for g in detect_groups(struct):
            q = GROUP_CHARGE.get(g.kind)
            if q is None:
                continue
            # A protonated oxyanion (HCO3-, HSO4-) carries its H's charge too.
            n_h = sum(1 for i in g.indices if elems[i] == "H")
            if g.kind in ("CO3", "SO4", "SiO4"):
                q += n_h
            total += q
            breakdown[g.kind] = breakdown.get(g.kind, 0) + 1
            claimed.update(g.indices)

    for i, e in enumerate(elems):
        if i in claimed:
            continue
        s = model.state(e)
        total += s
        key = f"{e}({s:+d})"
        breakdown[key] = breakdown.get(key, 0) + 1

    return {
        "total": total,
        "breakdown": breakdown,
        "assumed": dict(model.assumed),
        "unknown": sorted(model.unknown),
    }


def layer_charge(n_divalent: int, n_trivalent: int) -> int:
    """Net charge of an LDH metal sheet.

    A brucite sheet M(II)(OH)2 is neutral; every M(III) substituted onto it
    adds one positive charge, which is exactly what the interlayer anions have
    to cancel."""
    return n_trivalent


def anions_for_balance(layer_q: int, anion: str) -> tuple[int, int]:
    """How many anions of a given type neutralise a layer charge.

    Returns (count, residual). A residual is non-zero when the charge does not
    divide evenly - for carbonate on an odd layer charge, say - and the caller
    must resolve it (adjust the M(III) count, or mix in a monovalent anion)
    rather than round it away."""
    q = ANION_CHARGE[anion]
    count = layer_q // abs(q)
    residual = layer_q - count * abs(q)
    return count, residual


ANION_CHARGE = {"CO3": -2, "SO4": -2, "Cl": -1, "OH": -1, "NO3": -1, "CLO4": -1}


def describe(result: dict) -> str:
    lines = [f"formal charge: {result['total']:+d}"]
    for k, v in sorted(result["breakdown"].items()):
        lines.append(f"    {k:<12} x {v}")
    if result["assumed"]:
        lines.append("  ASSUMED oxidation states (not derived):")
        for e, s in sorted(result["assumed"].items()):
            alts = TRANSITION_METAL.get(e, (None, []))[1]
            extra = f"   also common: {', '.join(f'{a:+d}' for a in alts)}" if alts else ""
            lines.append(f"    {e:<3} {s:+d}{extra}")
    if result["unknown"]:
        lines.append(f"  UNKNOWN elements, counted as 0: {', '.join(result['unknown'])}")
    return "\n".join(lines)


def main(argv=None):
    import argparse

    from cement_core import read_xyz

    ap = argparse.ArgumentParser(description="Formal charge of a structure.")
    ap.add_argument("structure")
    ap.add_argument(
        "--oxidation", action="append", default=[], metavar="EL=STATE",
        help="override an oxidation state, e.g. --oxidation Fe=2 (repeatable)",
    )
    args = ap.parse_args(argv)

    overrides = {}
    for item in args.oxidation:
        el, _, val = item.partition("=")
        overrides[el.strip()] = int(val)

    s = read_xyz(args.structure)
    model = ChargeModel(overrides)
    need = model.needs_input(s)
    res = formal_charge(s, model)
    print(f"\n{args.structure}: {len(s)} atoms  {s.composition()}")
    print(describe(res))
    if need:
        print(f"\n  Transition metals with no explicit state: {', '.join(need)}")
        print("  Pass --oxidation EL=STATE to set them.")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
