# Charge bookkeeping

## Why it is not an afterthought

A periodic cell with a net charge converges in any DFT code and returns a
number that means nothing. Charge balance is therefore built into construction:
the LDH anion count is derived from the sheet charge, and a surface that cannot
be neutralised raises rather than returning.

## Main group: known

Applied without asking, because they are unambiguous in cement chemistry:

H +1 · Li/Na/K +1 · Mg/Ca/Sr/Ba +2 · Al/B/Ga +3 · Si/C +4 · P +5 · S +6 ·
O −2 · F/Cl/Br/I −1 · N −3

## Transition metals: asked, or defaulted and labelled

Most common state used as a default, alternatives recorded:

| | default | also common |
|---|---|---|
| Ti | +4 | +3 |
| Cr | +3 | +6 |
| Mn | +2 | +3, +4 |
| Fe | +3 | +2 |
| Co, Ni | +2 | +3 |
| Cu | +2 | +1 |
| Zn, Cd | +2 | |
| Zr, Hf | +4 | |

`ChargeModel.needs_input(struct)` lists transition metals present with no
explicit override. **Ask about these.** Fe(II) versus Fe(III) on an LDH metal
sublattice changes the sheet charge and therefore the entire interlayer
composition — it is not a detail.

```python
model = ChargeModel({"Fe": 2})       # explicit, no assumption recorded
formal_charge(struct, model)
```

## Groups, not atoms — except when the anion is polymerised

Isolated oxyanions are charged as units: OH −1, CO3 −2, SO4 −2. A protonated
oxyanion carries its hydrogen's charge too (HCO3 −1).

**Silicate is deliberately excluded.** In C-S-H the tetrahedra are polymerised
into chains, so bridging oxygens are shared between two silicons. Charging each
SiO4 as a discrete −4 unit counts those oxygens twice and puts a neutral
tobermorite cell near −100. Silicate falls through to the atomwise sum, where
every oxygen is counted exactly once whether it bridges or not.

The same reasoning applies to any condensed anion: if the polyhedra share
corners, count atomwise.

## LDH balance

Sheet charge = number of trivalent cations. Anion count = sheet charge divided
by anion charge, and the division must be exact:

| Cell | Metals | 2:1 → M(III) | CO3 needed |
|---|---|---|---|
| 4×6, 2 layers | 96 | 32 | 16 ✓ |
| 4×4, 2 layers | 64 | 21 | refused — odd |

When it does not divide, adjust the cell or ratio, or use a monovalent anion.

## Which oxygen carries the proton: bond valence

Formal charge tells you *how many* protons a cell needs. It does not tell you
*where* they go, and the builders decide that with a coordination heuristic —
`protonate_to_neutral` in `build_tobermorite11.py` takes the terminal oxygens, those
bonded to at most one Si, with the fewest additional cation contacts.

The rigorous check on that choice is a **bond-valence sum**. Compute each
oxygen's BVS with hydrogen excluded; a site that is already saturated cannot
take a proton, and a site that is short by ~1 valence unit must have one.

Worked example, tobermorite 9 Å (ICSD 87689), H excluded:

```
OH6 (= ICSD O9)   V = 1.03          -> short by ~1.0, needs exactly one H
O1..O8            V = 1.92 - 2.26   -> saturated, no H
```

Exactly one oxygen site in the whole structure asks for a proton, which settles
the assignment without appealing to the file's own header or to any heuristic.
The same site turns out to be a **silanol** — 1.669 Å from Si2 — not a free
hydroxide, which matches the 11 Å anomalous case where the protonated O6 is also
a terminal chain oxygen rather than an interlayer OH⁻.

Two things this method is good for:

- **Confirming a protonation heuristic picked the right site.** The heuristic
  reasons from coordination; BVS reasons from bond lengths. Agreement between
  two independent arguments is worth much more than either alone.
- **Catching a structure that cannot be neutralised sensibly.** If no oxygen is
  meaningfully under-bonded but the cell still carries charge, the composition
  or the occupancy resolution is wrong, not the protonation.

**Where the H points is a separate question, and geometry alone can get it
wrong.** In that same structure the H makes a clean near-linear bond to O7
(H···O 2.081 Å, 168.3°) — but O7 is the most *over*-bonded oxygen (BVS 2.26),
while the most under-bonded one, O8 (1.64), sits at O···O 2.767 Å and is
reachable only through an 81° angle. Relaxation may well rotate the proton
toward O8. Check H orientation after optimization rather than trusting the
placement.
