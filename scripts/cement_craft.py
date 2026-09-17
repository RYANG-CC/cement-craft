#!/usr/bin/env python3
"""
cement_craft command line: build a cement phase, optionally cut a surface.

  python cement_craft.py ldh --divalent Mg --trivalent Al --ratio 2 \
      --anion CO3 --water-per-anion 4 -o mg2al1_co3.xyz

  python cement_craft.py afm --substitute Al=Fe:2 -o afm_fe.xyz
  python cement_craft.py aft --water 40 -o aft_dry.xyz
  python cement_craft.py hydrogarnet --surface-axis a -o katoite_100.xyz
  python cement_craft.py csh --al-si 0.1 -o cash.xyz

  python cement_craft.py ldh --ratio 2 --surface-axis c \
      --surface-fraction 0.6 --counter-anion CO3 -o ldh_slab.xyz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cement_core import write_xyz  # noqa: E402
from cement_charge import ChargeModel, describe, formal_charge  # noqa: E402
from cement_phases import (  # noqa: E402
    build_csh, build_ldh, load_reference, set_interlayer_water, substitute,
    summary,
)
from cement_surface import build_surface, solvate  # noqa: E402

AXES = {"a": 0, "b": 1, "c": 2, "x": 0, "y": 1, "z": 2}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build cement phase structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("phase", choices=["ldh", "afm", "aft", "hydrogarnet", "csh"])
    ap.add_argument("-o", "--output", help="write the structure here (.xyz)")

    g = ap.add_argument_group("LDH")
    g.add_argument("--divalent", default="Mg")
    g.add_argument("--trivalent", default="Al")
    g.add_argument("--ratio", type=float, default=2.0,
                   help="M(II):M(III) on the metal sublattice")
    g.add_argument("--anion", default="CO3",
                   help="interlayer anion: CO3, SO4, Cl, OH, NO3")
    g.add_argument("--water-per-anion", type=float, default=4.0)
    g.add_argument("--na", type=int, default=4, help="sheet repeats along a")
    g.add_argument("--nb", type=int, default=6, help="sheet repeats along b")
    g.add_argument("--layers", type=int, default=2)
    g.add_argument("--basal", type=float, help="override basal spacing (A)")
    g.add_argument("--no-swell", action="store_true",
                   help="do not expand the gallery to fit the requested water")

    g2 = ap.add_argument_group("all phases")
    g2.add_argument("--supercell", default="1,1,1")
    g2.add_argument("--substitute", action="append", default=[],
                    metavar="FROM=TO:N",
                    help="cation substitution, e.g. Al=Fe:2 (repeatable)")
    g2.add_argument("--water", type=int,
                    help="set the interlayer water count")
    g2.add_argument("--al-si", type=float, default=0.0,
                    help="C-S-H only: Al/Si substitution giving C-A-S-H")
    g2.add_argument("--ca-si", type=float,
                    help="C-S-H only: target Ca/Si ratio")
    g2.add_argument("--oxidation", action="append", default=[],
                    metavar="EL=STATE",
                    help="fix a transition-metal oxidation state, e.g. Fe=2")
    g2.add_argument("--seed", type=int, default=0)

    g3 = ap.add_argument_group("surface")
    g3.add_argument("--surface-axis", choices=list(AXES),
                    help="cut a slab normal to this axis")
    g3.add_argument("--surface-fraction", type=float, default=0.6,
                    help="fraction of the cell kept as slab")
    g3.add_argument("--surface-origin", type=float, default=0.0)
    g3.add_argument("--vacuum", type=float, default=15.0)
    g3.add_argument("--counter-anion",
                    help="anion used to balance residual surface charge")
    g3.add_argument("--solvate", action="store_true",
                    help="fill the vacuum with water at bulk density, "
                         "making a solid-solvent interface")
    g3.add_argument("--n-water", type=int,
                    help="explicit reservoir water count (default: bulk density)")
    g3.add_argument("--counter-ion", help="ion added to the reservoir")
    g3.add_argument("--n-counter-ion", type=int, default=0)

    ap.add_argument(
        "--no-wrap-molecules", action="store_true",
        help="write raw coordinates instead of assembling each water/hydroxyl "
             "into one contiguous unit. Cosmetic only - a PBC calculation sees "
             "an identical structure either way; the default assembles them "
             "so they don't render as bonds stretched across the cell.",
    )

    args = ap.parse_args(argv)

    overrides = {}
    for item in args.oxidation:
        el, _, val = item.partition("=")
        overrides[el.strip()] = int(val)
    model = ChargeModel(overrides)

    if args.phase == "ldh":
        struct = build_ldh(
            divalent=args.divalent, trivalent=args.trivalent, ratio=args.ratio,
            anion=args.anion, n_a=args.na, n_b=args.nb, n_layers=args.layers,
            water_per_anion=args.water_per_anion, basal=args.basal,
            seed=args.seed, auto_basal=not args.no_swell,
        )
    elif args.phase == "csh":
        sc = tuple(int(x) for x in args.supercell.split(","))
        struct = build_csh(ca_si=args.ca_si, al_si=args.al_si, seed=args.seed,
                           supercell=sc)
    else:
        name = {"afm": "AFm", "aft": "AFt",
                "hydrogarnet": "hydrogarnet"}[args.phase]
        struct = load_reference(name)
        sc = tuple(int(x) for x in args.supercell.split(","))
        if sc != (1, 1, 1):
            struct = struct.repeat(*sc)

    for spec in args.substitute:
        frm, _, rest = spec.partition("=")
        to, _, n = rest.partition(":")
        struct = substitute(struct, frm.strip(), to.strip(),
                            count=int(n), seed=args.seed)

    if args.water is not None:
        struct = set_interlayer_water(struct, args.water, seed=args.seed)

    print(f"\n{args.phase.upper()}")
    print(summary(struct, model))

    need = model.needs_input(struct)
    if need:
        print(f"\n  Transition metals present with no stated oxidation state: "
              f"{', '.join(need)}")
        print("  The most common state was assumed. Pass --oxidation EL=STATE "
              "to set it explicitly.")

    if args.surface_axis:
        axis = AXES[args.surface_axis]
        print(f"\nSURFACE normal to {args.surface_axis}")
        struct, report = build_surface(
            struct, axis=axis,
            thickness=struct.cell[axis, axis] * args.surface_fraction,
            origin=args.surface_origin, vacuum=args.vacuum,
            oxidation=overrides, counter_anion=args.counter_anion,
            verbose=True,
        )
        print(summary(struct, model))

        if args.solvate:
            print("\nSOLVATION")
            solvate(struct, axis=axis, n_water=args.n_water,
                    counter_ion=args.counter_ion,
                    n_counter_ion=args.n_counter_ion, seed=args.seed,
                    verbose=True)
            print(summary(struct, model))

    if args.output:
        write_xyz(args.output, struct, wrap_molecules=not args.no_wrap_molecules)
        print(f"\n  wrote {args.output}")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
