# Fractional occupancy in CIFs

Experimental CIFs routinely carry `_atom_site_occupancy` below 1. An atomistic
simulation cannot use that directly, and the reason is worth being precise
about: **an occupancy is a statistical statement about a crystal, not a
description of one cell.** A site at 0.5 is filled in half the cells of the
sample; a simulation cell has to decide.

So a disordered CIF **does not define a unique ordered structure**. Different
orderings are all consistent with the same diffraction data and give different
energies. Picking one silently hides a choice that changes the answer.

## What the reader does

```python
read_cif(path)                              # strict: refuses, and explains
read_cif(path, occupancy="round")           # keep occ >= threshold
read_cif(path, occupancy="round", threshold=0.7)
read_cif(path, occupancy="ignore")          # face value; duplicates split sites
```

**`strict` is the default** and raises `OccupancyError`, listing the partial
sites and the smallest supercell in which the occupancies become whole numbers.
That number is the useful diagnostic:

- **small (≤ 8)** — an ordered approximant exists at a practical size. Build
  that supercell, enumerate distinct orderings, and compare several.
- **large** — the occupancies are not simple fractions and no practical
  supercell reproduces them exactly. Round, or go back to the refinement.

Real example, a stratlingite CIF:

```
6 of 11 sites have partial occupancy (O @ 0.500, Al @ 0.285, Si @ 0.285,
Al @ 0.270, Si @ 0.270, Wa @ 0.250). A disordered CIF does not define a unique
ordered structure - the smallest integral supercell is 308x, i.e. these
occupancies are not simple fractions and no practical ordered approximant
reproduces them exactly.
```

Those 0.27/0.285 pairs are a mixed Al/Si site — the classic case where the
refinement reports an average and the modeller has to choose an arrangement.

## Two different problems that look alike

**Split sites** — two positions 0.3–0.8 Å apart, occupancies summing to ~1.
This is *one* atom whose position is uncertain or dynamically averaged. Keep
one; taking both duplicates the atom and creates an impossible contact.

**Genuine partial occupancy** — a site at 0.25 with no partner nearby. The site
is empty three quarters of the time. This needs a supercell and an ordering
decision, or a vacancy.

`occupancy="ignore"` is only defensible for the first case, and only after
checking the geometry.

### A split site the CIF does not admit to

Worked example, tobermorite 9 Å (ICSD 87689). The CIF lists **occupancy 1.0 for
every Ca**, which reads as a normal fully occupied structure and implies Ca12
per parent cell. It is wrong, and only the geometry shows it: Ca2 and its
inversion image sit **0.593 Å apart**, far too close to both be occupied. Ca2 is
a split site at effective occupancy ½, and the correct content is **Ca10**.

The lesson generalises: **occupancy 1.0 everywhere is not proof of an ordered
structure.** Run a close-contact check on any newly read CIF — a heavy-atom pair
under ~0.8 Å is a split site being reported as two real atoms, whatever the
occupancy column says.

```python
close_contacts(struct, cutoff=0.8)   # catches split sites the CIF hides
```

Resolving it is an ordering decision like any other. Dropping one of each pair
consistently (the "A" arrangement) makes the ordering periodic with the *parent*
cell, so a supercell built on it adds no ordering freedom and the alternative
arrangements go unsampled — worth knowing before concluding anything from a
supercell. Consistent dropping also removes the inversion centre, so the result
is genuinely P1.

## Recommended workflow

1. Read strict. Look at what is partial and why.
2. If the partial sites are a **mixed cation site** (Al/Si, Mg/Al, Fe/Al),
   decide the ratio you want and substitute deliberately with
   `substitute(...)`, which places them randomly with a seed.
3. If they are **vacancies or interlayer disorder**, build the supercell and
   enumerate.
4. **Generate several configurations and compare.** One ordering is a sample of
   the disorder, not the structure. Relax each; the spread tells you how much
   the disorder matters.
5. Check composition and formal charge afterwards — rounding or ordering almost
   always changes both.

## Charge is the check that catches mistakes

After resolving occupancy, run the charge analysis. Dropping a partial anion
site or over-filling a mixed cation site shows up immediately as a non-zero
formal charge, which is far easier to see than a subtly wrong stoichiometry.


## Missing hydrogen on water sites

X-rays scatter weakly off hydrogen, so refinements routinely locate the water
oxygen and not its protons. Such sites arrive labelled `Wa`, `Ow`, `OW` or `W`
with nothing attached. That is a complete diffraction model and an unusable
simulation input: the missing protons carry charge and make the hydrogen bonds
that hold the water in place.

```python
s = read_cif(path, occupancy="round")
add_water_hydrogens(s, verbose=True)
```

Each bare oxygen gets two hydrogens at proper geometry — 0.96 Å, 104.5° — and
the rigid H–O–H unit is rotated by trial, scored on how well its protons point
at nearby acceptor oxygens (2.5–3.4 Å) without clashing. Water in these phases
is held by its hydrogen bonds, so aiming the protons at acceptors starts them
far closer to their relaxed positions than an arbitrary orientation would.

If the oxygen is also bonded to a cation (Ca, Mg, Al, ...), it is being held
as an aqua ligand, not just hydrogen-bonded, and the lone pair is what forms
that bond — so the search is constrained to point the lone pair at the cation
first, then only the remaining free spin about that direction is chosen by
the acceptor/clash scoring above. A water with no nearby cation is unaffected.

The oxygen is relabelled (`Wa` → `O`) once it is a real water, since the
refinement label means nothing to a simulation code. Pass `relabel=None` to
keep it, for instance when the label is doing duty as a CP2K `&KIND` name.

Worked example on a stratlingite CIF: 6 `Wa` sites become 6 water molecules,
12 hydrogens added, O–H 0.960 Å and H–O–H 104.5° exactly, every proton with an
acceptor within 3.0 Å (shortest H···O 2.09 Å), no close contacts.

**Check the charge afterwards.** Adding protons changes it: a bare oxygen
counts −2, a completed water 0. If the CIF was charge-balanced with the water
oxygens bare, it will not be once they are completed — which usually means the
original bookkeeping was treating them as water all along.
