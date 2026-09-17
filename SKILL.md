---
name: cement_craft
description: Construct atomic structures of cement and cementitious phases - layered double hydroxides (LDH/hydrotalcite) with any divalent/trivalent cation pair, adjustable interlayer water and interlayer anion (carbonate, sulfate, chloride, hydroxide, nitrate); AFm (monosulfate) and AFt (ettringite); hydrogarnet/katoite; and C-S-H and C-A-S-H from tobermorite-like references. Also cuts surfaces normal to any axis, balances the resulting charge, and reconstructs the exposed faces by protonation, adding OH and water. Use this skill whenever the user wants to build, substitute into, hydrate, or cut a surface from a cement phase, an LDH or hydrotalcite, ettringite, monosulfate, katoite, tobermorite or calcium silicate hydrate - and whenever a structure needs its formal charge balanced or its cut surface healed.
---

# cement_craft

Builds cement phase structures ready for DFT or classical simulation. Output is
a structure file plus a report of composition, groups and formal charge — the
structure is neutral by construction, not by adjustment afterwards.

## The rule that matters most

**A structure that is not charge-balanced is not a model.** A charged periodic
cell will converge in any DFT code and give a number that means nothing.

Charge is therefore never an afterthought here:

- an LDH sheet carries one positive charge per trivalent cation, and the
  interlayer anion count is *derived* from that, never chosen independently;
- a cut surface is healed by protonation, and any residual charge is closed
  explicitly with protons or counter-ions;
- when a target cannot be reached honestly, the builder **refuses with an
  explanation** rather than returning something plausible-looking.

Main-group oxidation states are treated as known. **Transition metals are not.**
If a structure contains one and the user has not stated its oxidation state, ask.
The most common state is used as a labelled default and reported as an
assumption — never folded silently into a total. Fe(II) versus Fe(III) in an LDH
changes the sheet charge, and therefore the entire interlayer composition.

## Quick start

```bash
python scripts/cement_craft.py ldh --divalent Mg --trivalent Al --ratio 2 \
    --anion CO3 --water-per-anion 4 -o mg2al1.xyz

python scripts/cement_craft.py aft --substitute Al=Fe:2 -o aft_fe.xyz
python scripts/cement_craft.py csh --al-si 0.125 -o cash.xyz
python scripts/cement_craft.py hydrogarnet --surface-axis a -o katoite_slab.xyz
```

`--help` lists every option. All scripts are importable too:

```python
from cement_phases import build_ldh, load_reference, substitute
from cement_surface import build_surface
```

## Phases

| Phase | Route | Notes |
|---|---|---|
| **LDH** | procedural | Any M(II)/M(III) pair, any ratio, anion CO3/SO4/Cl/OH/NO3, tunable water |
| **AFm** | reference cell | Ca4Al2(SO4)(OH)12·4H2O ×3; substitute or re-hydrate |
| **AFt** | reference cell | Ca6Al2(SO4)3(OH)12·24H2O ×2 (ettringite) |
| **hydrogarnet** | reference CIF | Ca3Al2(OH)12 ×8 (katoite), cubic Ia-3d |
| **C-S-H** | reference series | Ca/Si 0.667 / 0.75 / 0.833 / 1.25, 11 Å and 14 Å tobermorite; Al→Si gives C-A-S-H |

**LDH is generated from geometry** — a brucite-like M(OH)2 sheet with M(III)
substituted onto it — because that is what makes arbitrary cation pairs,
ratios, anions and hydration states reachable. The others start from a
crystallographic reference, because their frameworks are specific and building
them from scratch would mean inventing coordinates.

Details of every phase, its cell and what can be varied:
**`references/phases.md`**.

## Building an LDH

```python
s = build_ldh(divalent="Mg", trivalent="Al", ratio=2.0,
              anion="CO3", n_a=4, n_b=6, n_layers=2, water_per_anion=4.0)
```

The cell is an `n_a × n_b` orthogonal supercell of the hexagonal sheet, giving
`2·n_a·n_b` metals per layer:
`a = n_a·a_hex·√3`, `b = n_b·a_hex`, with `a_hex ≈ 3.07 Å` for Mg-Al.

Two things follow from charge balance and are worth knowing before choosing a
cell:

- The trivalent count must be divisible by the anion charge. A 4×4 cell at 2:1
  gives 21 trivalent cations, which carbonate cannot balance; the builder says
  so and suggests adjusting the cell or using a monovalent anion.
- If the requested water does not fit at the tabulated basal spacing, the
  gallery **expands** rather than silently dropping molecules. That is the real
  behaviour of these materials — basal spacing is set by hydration state — and
  the swelling is reported in `info["basal_swelling"]` so it can be checked.

## Cutting a surface

```python
slab, report = build_surface(struct, axis=2, thickness=..., vacuum=15.0,
                             counter_anion="CO3")
```

Three steps, each of which can be run alone:

1. **Cut** — molecular units are kept whole. A water or carbonate whose centre
   is inside the window is retained entirely; one outside is removed entirely.
   Slicing on atom position alone leaves half-molecules, which is a defect and
   not a surface.
2. **Protonate** — a surface oxygen bonded to **one** framework metal becomes a
   **water molecule**; one bonded to **two or more** remains a **hydroxyl**.
   That is what returns every framework metal to its bulk coordination.
3. **Neutralise** — residual charge is closed with protons where possible
   (a surface hydroxyl and a surface water differ by exactly one), and with
   counter-anions when the slab is left positive, which is what happens when a
   cut through a layered phase discards interlayer anions.

If neutrality cannot be reached, `SurfaceChargeError` is raised with what was
tried. Do not work around it by returning the charged slab.

## Solid-solvent interfaces

A slab in vacuum answers a different question from a slab in contact with
water, and for cement phases it is usually the wetted surface that matters.

```python
solvate(slab, axis=2, counter_ion="Cl", n_counter_ion=2)
```

The vacuum region is filled at **bulk water density** (0.0334 molecules/A^3)
rather than by a molecule count, so the reservoir is physically meaningful
whatever the slab area is. A clearance (`gap`, default 2.6 A) is left between
the surface and the first water so the reservoir starts at a hydrogen-bond
distance. Counter-ions can be dissolved in it for a surface that is charged in
solution.

Build the slab with enough vacuum first - `--vacuum 22` or more - or there is
nothing to fill.

```bash
python scripts/cement_craft.py hydrogarnet --surface-axis a     --vacuum 22 --solvate -o katoite_water.xyz
```

Surface rules and the reasoning behind them: **`references/surfaces.md`**.

## Assembling molecules for viewing

```python
struct.wrap_molecules()
```

Purely cosmetic - a PBC calculation sees an identical structure either way -
but it stops a water or hydroxyl from rendering as a bond stretched across the
whole cell when its atoms land on opposite faces. Every finite molecule is
shifted as one rigid unit; bond lengths are exactly preserved. `cement_craft.py`
applies this by default (`--no-wrap-molecules` to disable). Polymerised species
(silicate chains) are left to the ordinary per-atom wrap, since a bridging
oxygen belongs to two tetrahedra at once and a chain spanning the cell is
genuinely continuing into the next periodic image, not broken. See
`references/surfaces.md`.

## Reporting a structure

Always report composition, perceived groups and formal charge — the same three
things the builders print. A structure whose group count looks wrong (sheet
hydroxyls perceived as water, say) is usually a geometry error, and the group
census catches it before a calculation does.

```python
from cement_phases import summary
print(summary(struct))
```

## Honest limits

- **C-S-H comes from a reference series, not interpolation.** Ca/Si 0.667,
  0.833 and 1.25 each have their own cell. A ratio between them raises with the
  list of what is available, because reaching one means removing bridging
  silicate tetrahedra and charge-compensating — a modelling decision with
  several published conventions. The 0.667 reference has no cell on record and
  refuses until one is supplied. Al→Si substitution for C-A-S-H *is* supported,
  and the resulting charge imbalance is reported for the caller to balance.
- **Structures are starting points, not equilibrium geometries.** Packing aims
  for hydrogen-bond separations and rejects anything below a physical floor,
  but interlayer water is placed by a rigid-body packer, not equilibrated.
  Relax with DFT or MD before drawing conclusions.
- **Water from a CIF usually has no hydrogen.** X-rays barely see it, so
  `Wa`/`Ow` sites arrive as lone oxygens. `add_water_hydrogens()` completes
  them at 0.96 Å / 104.5°, orienting each toward nearby acceptors, and
  relabels the oxygen. When the oxygen is also bonded to a cation (Ca, Mg,
  ...) — an aqua ligand, not just a hydrogen-bonded water — the lone pair is
  pointed at that cation instead, since that is the bond the ligand is
  actually making; only the free spin about that direction is then chosen by
  acceptor/clash scoring. Re-check the charge afterwards — a bare oxygen
  counts −2 and a completed water 0.
- **Disordered CIFs are refused, not guessed at.** A CIF with partial
  occupancies does not define a unique ordered structure, so `read_cif` raises
  by default and reports the smallest supercell that would make the
  occupancies integral. Resolve it deliberately — see `references/occupancy.md`
  — and compare several orderings rather than trusting one.
- **Cation substitution is random.** Site preference is usually the thing a
  study is trying to determine, so imposing an ordering would prejudge it. Pass
  a `seed` for reproducibility, and build several configurations.

## Reference files

| File | Read when |
|---|---|
| `references/phases.md` | Building or substituting into any phase; cells, formulas, what varies |
| `references/surfaces.md` | Cutting slabs, protonation rules, charge balancing, solvation |
| `references/charge.md` | Oxidation states, group charges, why silicate is handled atomwise, bond valence for protonation sites |
| `references/occupancy.md` | Reading CIFs with partial occupancy; ordering, split sites, supercells |

## Scripts

| Script | Purpose |
|---|---|
| `scripts/cement_craft.py` | Command line front end |
| `scripts/cement_phases.py` | Phase builders, substitution, hydration, packing |
| `scripts/cement_surface.py` | Slab cutting, protonation, neutralisation |
| `scripts/cement_charge.py` | Oxidation states and formal charge |
| `scripts/cement_core.py` | Structure, I/O (xyz/extxyz/CIF), PBC, groups, layers |
| `scripts/build_tobermorite11.py` | Builds the bundled 11 Å normal/anomalous C-S-H references from their COD CIFs |

NumPy is the only requirement. ASE is used, if installed, as a fallback for
CIF files the built-in reader cannot handle.
