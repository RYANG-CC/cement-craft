# Cement phases: structures and what can be varied

## Layered double hydroxide (LDH / hydrotalcite)

`[M(II)_{1-x} M(III)_x (OH)2]^{x+} · (A^{n-})_{x/n} · mH2O`

A brucite-like sheet of edge-sharing M(OH)6 octahedra. Replacing M(II) with
M(III) leaves the sheet positively charged, and that charge is balanced by
anions in the gallery between sheets, along with water.

**Built procedurally.** Metals sit on a triangular lattice; hydroxyl oxygens
occupy the centroids of alternating triangles above and below, which is what
makes every metal six-coordinate. In-plane metal-to-oxygen distance is
`a_hex/√3`, so the out-of-plane offset follows from the M–O bond length.

| Parameter | Default | Notes |
|---|---|---|
| `divalent` / `trivalent` | Mg / Al | Any pair; M–O distances tabulated for common cations |
| `ratio` | 2.0 | M(II):M(III). 2, 3 and 5 are the usual synthesis targets |
| `anion` | CO3 | CO3, SO4, Cl, OH, NO3 |
| `water_per_anion` | 4.0 | ~4 is typical for carbonate LDH |
| `n_a`, `n_b` | 4, 6 | Sheet repeats; 2·n_a·n_b metals per layer |
| `n_layers` | 2 | Stacked sheets |
| `a_hex` | 3.07 Å | Measured 3.07–3.09 for Mg-Al |

Basal spacing by anion (Å): CO3 7.78, SO4 8.90, Cl 7.80, OH 7.60, NO3 8.80.
These set gallery volume and therefore how much water fits. With `auto_basal`
(default on), the gallery expands when the requested water will not fit, and
the swelling is reported — real LDH basal spacing does depend on hydration.

**Charge divisibility.** The trivalent count must be a multiple of the anion
charge. A 4×6 cell over two layers has 96 metals, so 2:1 gives 32 trivalent
cations and 16 carbonates — exact. A 4×4 cell gives 21, which carbonate cannot
balance; the builder refuses and says so.

**Validation.** A Mg2Al1-CO3 build at 4×6×2 layers reproduces the composition of
an independently prepared reference structure atom for atom: Mg64 Al32 C16 O306
H324, 742 atoms.

Planar anions (carbonate, nitrate) are placed lying flat between the sheets
with free rotation about the stacking axis and a small tilt. A freely oriented
carbonate is 2.6 Å across and simply does not fit a 3.7 Å gallery — and
flat-lying is what diffraction shows in any case.

## AFm — monosulfate

`Ca4Al2(SO4)(OH)12·4H2O` (bundled reference is ×3)

Cell 5.759 × 9.974 × 26.795 Å, all angles 90°. 141 atoms:
Ca12 Al6 S3 O60 H60. Groups: 3 SO4, 36 OH, 12 H2O. Neutral.

A layered calcium aluminate hydrate: positively charged
`[Ca2Al(OH)6]^+` layers with sulfate and water between them. Structurally an
LDH with Ca and Al on the metal sublattice, which is why the same substitution
and hydration tools apply.

Vary by substituting Al → Fe(III) (Fe-AFm) and by resetting the interlayer
water count.

## AFt — ettringite

`Ca6Al2(SO4)3(OH)12·24H2O` (bundled reference is ×2)

Cell 11.229 × 11.229 × 21.478 Å, γ = 120°. 238 atoms:
Ca12 Al4 S6 O96 H120. Groups: 6 SO4, 24 OH, 48 H2O. Neutral.

Columns of `[Ca6Al2(OH)12·24H2O]^{6+}` running along c, with sulfate and water
in the channels between them. The hexagonal cell matters: cutting a surface
normal to a or b is not equivalent to cutting normal to c.

## Hydrogarnet — katoite

`Ca3Al2(OH)12` (bundled reference is ×8, the full cubic cell)

Cubic Ia-3d, a = 12.5389 Å. 232 atoms: Ca24 Al16 O96 H96, all oxygen as
hydroxyl. Neutral.

The silicon-free end member of the hydrogarnet series
`Ca3Al2(SiO4)_x(OH)_{12-4x}`. Each SiO4 replaces four OH, so introducing
silicon means removing four hydroxyls per silicon added — not a substitution
the builders perform automatically, because the site choice matters.

Read from CIF with symmetry operators expanded, so the asymmetric unit becomes
the full cell.

## C-S-H and C-A-S-H

C-S-H is a compositional family, not one structure, and **Ca/Si is the variable
that matters most** — it sets the silicate chain length and how many bridging
tetrahedra are missing. Values of 0.67, 0.83 and ~1.25 are all common, so each
is served by its own reference cell rather than derived from another.

### The tobermorite series

From Churakov 2009a (*Eur. J. Mineral.* **21**, 261–271, normal 11 Å) and
Churakov 2009b (*Am. Mineral.* **94**, 156–165, anomalous 11 Å):

| Phase | Formula | Ca/Si |
|---|---|---|
| 14 Å plombierite | Ca5Si6O16(OH)2·7H2O | 0.833 |
| 11 Å normal | Ca5−xSi6O17−2x(OH)2x·5H2O | 0.75 as simulated |
| 11 Å anomalous | Ca4Si6O15(OH)2·5H2O | 0.667 |
| 9 Å riversideite | Ca5Si6O16(OH)2 | 0.833 |

14 Å converts to 11 Å on heating to 353–373 K; normal 11 Å loses interlayer
water at 573 K and goes to 9 Å, while anomalous 11 Å does not collapse.

**Chain topology differs across the series, and it is the thing most often got
wrong.** Tobermorite 9 Å (riversideite) has **single dreierketten — every Si is
Q², O:Si exactly 3.000**. Double chains belong to *clinotobermorite* and to
tobermorite 11 Å, not to the 9 Å phase. In a model validated against ICSD 87689
the chains run along **b** with a 7.303 Å Si1→Si2→Si3 repeat, and Si2···Si2
across the interlayer is 4.131 Å with **no bridging oxygen** — that missing
bridge is what makes it single-chain. Do not "repair" it into a double chain.

Same model, for reference when checking one: Ca1/Ca3 sevenfold in the Ca sheet,
Ca2 sixfold in the interlayer, Ca BVS 1.90–2.01, density 2.865 g/cm³ against
2.86 reported, basal spacing 9.376 Å along c*, anhydrous (largest void 1.95 Å).
Unrelaxed X-ray coordinates give an Si–O spread of 1.554–1.677 Å; that is
refinement noise, and relaxation tightens it to ~1.61–1.64. Optimize before
computing any property.

Churakov's validated orthorhombic AIMD supercells, useful when building an
11 Å cell from crystallography rather than from a bundled file:

| | supercell (Å) | angles |
|---|---|---|
| anomalous 11 Å | 11.265 × 14.792 × 22.487 | 90/90/90 |
| normal 11 Å | 11.265 × 14.792 × 22.680 | 90/90/90 |

Both are 2 × 2 × 1 of the monoclinic **B11m** cell
(13.47 × 14.74 × 22.487, γ = 123.25°), the doubling being along *a*.

### Bundled references

| Ca/Si | File | Cell (Å) | Atoms |
|---|---|---|---|
| 0.667 | `CSH_11A_anomalous.xyz` | 11.2648 × 14.7700 × 22.4870 | 352 |
| 0.75 | `CSH_11A_normal.xyz` | 11.2688 × 14.7380 × 22.6800 | 348 |
| 0.833 | `CSH_14A.xyz` | 14.85 × 13.47 × 27.987 | 392 |
| 1.25 | `CSH_v_CaSi125.xyz` | 11.4769 × 15.0660 × 19.2663 | 296 |

### The two 11 Å bulk references, and the disorder behind the normal one

Both are built by `scripts/build_tobermorite11.py` from Merlino, Bonaccorsi &
Armbruster (2001, *Eur. J. Mineral.* **13**, 577–590), the same paper that
supplies both polytypes in one refinement — CC0 on COD, sources kept under
`assets/cif_sources/`:

- **anomalous**, COD 9005498 (Wessels mine) — the Ca2 interlayer site is
  simply **absent**, not disordered. Every occupancy is 1; nothing to resolve.
- **normal**, COD 9005499 (Bašenov, Urals) — Ca2 (occ. 0.25) and the Wat1/Wat3
  waters (occ. 0.5) *are* genuinely disordered. Each expands under the listed
  symmetry operators into two locations related by the lattice centring
  (≈10–12 Å apart — not a split site), and each of those two locations is
  itself a mirror-related split pair (≈0.9–2.0 Å apart). `occupancy × raw
  copies` says how many of the two locations are real (1 of 2 for Ca2, both
  for Wat1/Wat3); the build keeps the lowest-indexed candidate at each real
  location ("ordering A" — the same explicit, labelled, one-of-several-valid
  convention used for the tobermorite 9 Å Ca2 split site) and drops the rest.
  This is one ordering, not the only one; re-run with a different tie-break
  to sample another.

Both close a real gap in the *published* formula: the CIF's own
`_chemical_formula_sum` (Ca2Si3H5O11 for anomalous, Ca2.25Si3H7O11 for normal)
is not charge-neutral even before any disorder resolution — the X-ray
refinement locates some water hydrogens explicitly but not every proton the
mineral actually needs. `protonate_to_neutral` closes the rest the same way
`add_water_hydrogens`/CIF-hydrogen completion does elsewhere in this project:
by charge, not by trusting the deposited formula. Check the resulting Ca/Si
against the target after building; that is the number the formula's own
stoichiometry constrains, independent of how the extra protons are placed.

**Both are confirmed double dreierketten**, unlike the single-chain 9 Å phase
(see `csh-tobermorite9a-validated.md`): `csh_chain_report()` shows Q2 at the
two "paired" Si sites and **Q3 at the bridging Si** (1 of every 3, not the 2 a
naive single-chain picture would predict), the quantitative signature of the
two flanking chains condensing at the bridging tetrahedron. Confirmed
identical between both forms, as expected — normal and anomalous differ only
in interlayer Ca/water content, not in the silicate framework itself.

```python
build_csh(ca_si=0.667)   # anomalous 11 A, bulk, exact
build_csh(ca_si=0.75)    # normal 11 A, bulk, ordering A
csh_chain_report(build_csh(ca_si=0.667))   # Ca/Si, Qn, chain verdict
```

### `CSH_CaSi067.xyz` is a surface, not this bulk phase

A *second*, unrelated 152-atom structure at the same nominal Ca/Si 0.667 also
lives in the archive, and it is tempting to confuse it with the bulk anomalous
phase above. It is not the same thing. It is a **single layer cut from the
0.833 phase with the interlayer stripped** — a simplified basal surface — and
its Ca/Si follows from removing interlayer calcium, not from the anomalous
form's actual missing-Ca2 stoichiometry. `build_csh(ca_si=0.667)` returns the
real bulk phase above; this surface is deliberately not in `CSH_REFERENCES`
and has to be rebuilt:

```python
par  = build_csh(ca_si=0.833)
slab = cut_slab(par, axis=2, thickness=par.cell[2, 2] / 2, vacuum=15.0)
strip_interlayer(slab, axis=2, keep_water=10)
# -> 152 atoms, Ca16 Si24 O80 H32, Ca/Si 0.667, neutral
```

That reproduces the archive structure's composition exactly. The in-plane cell
is the parent's; the c axis is chosen for the water and vacuum wanted, which is
why no fixed cell belongs with it. The arithmetic: against half the parent
cell, the cut removes 4 Ca, 16 O and 24 H, which resolves as
4 Ca²⁺ + 8 H₂O + 8 OH⁻ — charge-neutral, since +8 from the calcium cancels −8
from the hydroxide. Silicon is exactly halved, one layer of two.

```python
build_csh(ca_si=0.833)                 # nearest reference, exact ratio reported
build_csh(ca_si=1.25, al_si=0.125)     # C-A-S-H
```

**Ratios between the references are not interpolated.** Asking for one raises
with the list of what is available, because reaching an intermediate Ca/Si
means removing bridging silicate tetrahedra and charge-compensating, and there
is more than one published convention for which sites to remove.

### Check Ca/Si and chain connectivity on every CSH build

`csh_chain_report(struct)` (in `cement_phases.py`) reports Ca/Si, the Qn
histogram and a plain-text verdict (single chain / double chain / isolated
dimers / cross-linked), and is run against every `CSH_REFERENCES` entry in
`verify.py`. Run it on any new or modified CSH structure too — a fully
connected single chain is all-Q2, a double chain adds Q3 at exactly the
bridging site, and defect C-S-H with bridging tetrahedra removed goes to
all-Q1 (`build_csh(ca_si=1.25)` is exactly this: MCL 2, isolated dimers). A
result that does not match one of those clean patterns is worth looking at
before trusting the structure.

**The 0.667 reference has no cell on record** and no `cp2k.inp` alongside it,
so it refuses until a cell is passed. Its comment line also carries an energy
within 0.008 Ha of a 392-atom structure, which cannot be its own — treat the
provenance as unconfirmed.

**Al substitution (C-A-S-H)** replaces that fraction of Si with Al. Al(III) on
a Si(IV) site leaves one negative charge per substitution, and the structure
reports the imbalance for the caller to close — extra Ca(2+), or a protonated
bridging oxygen.

Silicate charge is computed **atomwise**, not per tetrahedron: the chains share
bridging oxygens, and counting each SiO4 as a discrete −4 unit double-counts
them and puts a neutral cell near −100.

## Substitution and hydration, all phases

```python
substitute(struct, "Al", "Fe", count=2)        # or fraction=0.25
substitute(struct, "Mg", "Fe", fraction=0.5, label="Fe_A")
set_interlayer_water(struct, n_water=40)
```

Substitution sites are chosen at random with a seed, because site preference is
generally what a study sets out to determine. `label` writes a distinct kind
name (`Fe_A`, `Fe_B`), which is how an antiferromagnetic arrangement is set up
downstream in a DFT input.

Removing water drops the molecules furthest from the framework first, so the
ones that remain are the more strongly bound.
