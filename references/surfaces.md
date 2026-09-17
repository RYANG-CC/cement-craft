# Surfaces: cutting, protonation, charge balance

Cutting a hydrated mineral leaves two problems that must both be solved before
the slab is a model of anything: under-coordinated oxygens at the exposed
faces, and a net charge.

## 1. Cut

```python
slab = cut_slab(struct, axis=2, thickness=..., origin=0.0, vacuum=15.0)
```

Molecular units are kept whole. A water, carbonate or sulfate whose **centre**
lies inside the window is retained entirely; one whose centre lies outside is
removed entirely. Slicing on atom position alone leaves half-molecules at both
faces, which is a defect, not a surface.

The centroid is computed under the minimum image and anchored on the first
atom of the group, so a molecule straddling a periodic boundary does not
average to the middle of the cell.

Vacuum is opened along the cut axis after re-zeroing the slab. 15 Å is a
reasonable default; a slab with a dipole needs more.

**Choose the cut plane deliberately.** For a layered phase, where the window
falls decides whether the interlayer anions come with their sheets. A cut
through the middle of a gallery keeps them; one that clips the gallery edge
discards them and leaves the sheets strongly positive.

## 2. Protonate

```python
protonate_surface(slab, axis=2)
```

The rule these phases obey:

| Surface oxygen bonded to | Becomes |
|---|---|
| **one** framework metal | **H2O** — it was a bridging site that lost its partner |
| **two or more** framework metals | **OH** — still bridging, just needs its proton |
| already carrying 2 H | left alone |
| no framework metal | left alone (free water, or an anion oxygen) |

Bulk coordination is taken as the modal metal-coordination of oxygen in the
structure, so the check adapts to the phase rather than assuming a number.
Oxygens already at bulk coordination are untouched.

New hydrogens point away from the mean position of the bonded metals, biased
along the surface normal so they project into the vacuum rather than back into
the slab.

## 3. Neutralise

```python
neutralise(slab, axis=2, counter_anion="CO3")
```

Protons are the natural currency: a surface hydroxyl and a surface water differ
by exactly one, and both interconvert in contact with water. Sites are chosen
at the outermost faces first, so the interior keeps its bulk protonation.

- **Slab too negative** → protonate surface hydroxyls to water.
- **Slab too positive** → deprotonate surface water to hydroxyl.
- **Still positive with no water left** → add counter-anions just beyond each
  face. This is the normal outcome for an LDH: the sheets are intrinsically
  positive and it is anions, not protons, that balance them.

If neutrality still cannot be reached, `SurfaceChargeError` is raised with what
was tried and why. **Do not work around it** by accepting the charged slab —
a charged periodic cell converges to a meaningless number.

A residual charge that is not divisible by the counter-anion charge also
raises: +3 cannot be balanced by carbonate. Use a monovalent anion or move the
cut.

## Worked outcomes

| System | Cut | Result |
|---|---|---|
| Katoite, normal to a | 60% of cell | 8 singly-coordinated O → H2O, then −8 H → **neutral** |
| Mg2Al1-CO3 LDH, normal to c | 60% | cut leaves +64; 32 CO3 counter-ions → **neutral** |
| Mg2Al1-CO3 LDH, normal to b | 60% | cut leaves +4; 2 CO3 → **neutral** |
| Mg2Al1-CO3 LDH, normal to a | 60% | leaves +3 — **refused**, not divisible by carbonate |

All produced zero close contacts below 0.75 Å.

## After building

New hydrogens can land close to existing atoms, so a short rigid-body
relaxation nudges them clear. The slab is then sorted along the cut axis, which
is what makes an index-based constraint (`&FIXED_ATOMS` in CP2K, say) meaningful:
the bottom layers are the first indices. Sort before writing, never after.


## 4. Solvate — the solid-solvent interface

```python
solvate(slab, axis=2, counter_ion="Na", n_counter_ion=1)
```

Fills the vacuum above the slab with water at **bulk density**
(0.0334 molecules/A^3, about 1 g/cm^3). Specifying a density rather than a
count means the reservoir stays physical when the slab area changes; pass
`n_water` to override when matching an existing structure.

- `gap` (default 2.6 A) is the clearance between the topmost slab atom and the
  first water, so the reservoir begins at a hydrogen-bond distance instead of
  on top of the surface hydroxyls.
- Counter-ions are placed in the reservoir before the water, so they end up
  solvated rather than buried.
- The same grid packing and rigid-body relaxation used for interlayers applies,
  so the result is clash-free but **not equilibrated**. Run MD on the water
  before production: the reservoir is a starting configuration, not a liquid.

Worked example: a katoite slab cut normal to a with 22 A of vacuum gives a
19.4 A reservoir of 3050 A^3, filled with 102 water molecules at exactly
0.0334 /A^3. Result 450 atoms, neutral, no close contacts.

## File format note

Interface structures are sometimes written with an orthorhombic shorthand on
the comment line instead of an extended-xyz `Lattice=` block:

```
990
CELL 34.000000000 19.948375702 15.000000000 | notes | formal charge 0
```

`read_xyz` understands both. Recording the formal charge in the comment, as
that example does, is a habit worth keeping - it makes a structure's charge
state checkable long after the build.

## Assembling molecules for viewing (wrap_molecules)

A PBC calculation only ever sees minimum-image distances, so it makes no
difference whether a water molecule's H sits at fractional coordinate 0.98 or
-0.02 - both describe the same physical structure. But a plain viewer that
plots raw coordinates does not know that, and will draw a bond stretched clean
across the cell if the O and its H happen to land on opposite sides.

```python
struct.wrap_molecules()             # in place
write_xyz(path, struct, wrap_molecules=True)   # or via the CLI: --wrap-molecules
```

Each finite molecule (H2O, OH, an isolated CO3 or SO4) is shifted as a rigid
body by whole lattice vectors, chosen so it lands together in the primary
cell. Bond lengths and internal geometry are exactly preserved - this changes
nothing a calculation would see, only where a human looking at the file sees
the atoms sitting. `cement_craft.py` does this by default; pass
`--no-wrap-molecules` for raw builder coordinates.

**Polymerised species are deliberately left alone.** A silicate chain's
bridging oxygen belongs to two SiO4 tetrahedra at once, and a chain spanning
the whole cell is not a finite molecule in the first place - it genuinely
continues into the next periodic image. Trying to rigidly assemble two
overlapping groups around a shared atom just moves it to please whichever
group is processed last, leaving it many angstroms from the tetrahedron it
was supposed to also belong to. Those atoms are left to the ordinary per-atom
wrap instead, which keeps every atom inside the primary cell but does not try
to make the chain look unbroken - because within one unit cell, it isn't.
