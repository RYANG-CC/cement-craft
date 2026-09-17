#!/usr/bin/env python3
"""End-to-end check: every phase builds neutral and clash-free."""
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

from cement_charge import formal_charge                      # noqa: E402
from cement_core import close_contacts, detect_groups        # noqa: E402
from cement_phases import (                                  # noqa: E402
    CSH_REFERENCES, build_csh, build_ldh, csh_chain_report, load_reference,
    set_interlayer_water, substitute,
)
from cement_surface import build_surface                     # noqa: E402

ok = fail = 0


def check(name, fn):
    global ok, fail
    try:
        s = fn()
        q = formal_charge(s)["total"]
        cc = len(close_contacts(s, 0.75))
        g = collections.Counter(x.kind for x in detect_groups(s))
        good = (q == 0 and cc == 0)
        ok, fail = (ok + 1, fail) if good else (ok, fail + 1)
        print(f"  {'OK ' if good else 'BAD'} {name:26} {len(s):4d} at  "
              f"q={q:+d}  clash={cc}  {dict(g)}")
    except Exception as exc:
        fail += 1
        print(f"  ERR {name:26} {type(exc).__name__}: {str(exc)[:55]}")


print("BULK PHASES")
check("LDH Mg2Al1-CO3", lambda: build_ldh(ratio=2.0, anion="CO3",
                                          water_per_anion=4.0))
check("LDH Mg3Al1-SO4", lambda: build_ldh(ratio=3.0, anion="SO4",
                                          water_per_anion=4.0))
check("LDH Ni/Fe 3:1 CO3", lambda: build_ldh(divalent="Ni", trivalent="Fe",
                                             ratio=3.0, anion="CO3"))
check("LDH Mg2Al1-Cl", lambda: build_ldh(ratio=2.0, anion="Cl",
                                         water_per_anion=2.0))
check("AFm", lambda: load_reference("AFm"))
check("AFt", lambda: load_reference("AFt"))
check("hydrogarnet", lambda: load_reference("hydrogarnet"))
check("C-S-H 14A", lambda: build_csh())

print("MODIFICATIONS")
check("AFm + 2 Fe", lambda: substitute(load_reference("AFm"), "Al", "Fe",
                                       count=2))
check("AFt water -> 40", lambda: set_interlayer_water(load_reference("AFt"), 40))

print("SURFACES")
_b = build_ldh(ratio=2.0, anion="CO3", water_per_anion=2.0)
check("LDH slab perp c", lambda: build_surface(
    _b, axis=2, thickness=_b.cell[2, 2] * 0.6, counter_anion="CO3",
    verbose=False)[0])
_h = load_reference("hydrogarnet")
check("katoite slab perp a", lambda: build_surface(
    _h, axis=0, thickness=_h.cell[0, 0] * 0.6, verbose=False)[0])
_a = load_reference("AFt")
check("AFt slab perp c", lambda: build_surface(
    _a, axis=2, thickness=_a.cell[2, 2] * 0.6, verbose=False)[0])

print("CSH CHAIN CHECKS")
for ca_si in sorted(CSH_REFERENCES):
    try:
        s = build_csh(ca_si=ca_si)
        r = csh_chain_report(s)
        good = abs(r["ca_si"] - ca_si) < 0.01
        ok, fail = (ok + 1, fail) if good else (ok, fail + 1)
        print(f"  {'OK ' if good else 'BAD'} Ca/Si {ca_si:<5} -> {r['ca_si']:<7} "
              f"Qn={r['Qn']}  {r['verdict']}")
    except Exception as exc:
        fail += 1
        print(f"  ERR Ca/Si {ca_si:<5} {type(exc).__name__}: {str(exc)[:55]}")

print(f"\n{ok} passed, {fail} failed")
sys.exit(1 if fail else 0)
