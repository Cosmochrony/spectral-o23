# Reproduction code for O23, subsection "A proved obstruction on Bridge 1"

Two scripts.

`o23_filtration_stabiliser.py` computes the exact stabiliser of each O12 shell-span **level** inside the
Weil image of $\mathrm{SL}(2,\mathbb{Z}/q\mathbb{Z})$, reproducing the pure-Fourier lemma, the exact
level stabiliser $U \rtimes M_n \subseteq B(L)$, the exclusion of the **exceptional** binary polyhedral
groups $2T$, $2O$, $2I$, and the measured multiplier groups.

`o23_dicyclic_decomposition.py` reproduces the two dicyclic witnesses. Borel confinement does **not**
exclude dicyclic groups, which are metacyclic and carry faithful two-dimensional spinor representations;
this script exhibits levels whose stabiliser is $\mathrm{Dic}_q$ and decomposes the level under it. It
verifies the dicyclic presentation, the invariance of the level, the homomorphism defect of the lift over
**all** ordered pairs of the stabiliser, the character norm, the multiplicities, and the isotypic
projectors (idempotence, ranks, completeness), exiting non-zero if any check fails.

```bash
python3 o23_dicyclic_decomposition.py                                  # both deposited witnesses
python3 o23_dicyclic_decomposition.py --q 53 --block 47,21,32 --level 1
```

Expected output for each witness: $\dim W_{<1} = 9$, $|\mathrm{Stab}| = 4q$, dicyclic confirmed,
$\langle\chi,\chi\rangle = 5$, one one-dimensional and four two-dimensional constituents each of
multiplicity one, and $1 + 4\times 2 = 9$. Evidence is written to
`checkpoints/dicyclic_decomposition.json`.

**Status of the lift.** A genuine linear Weil representation exists for $q$ an odd prime (Weil 1964;
Gérardin 1977). That this Bruhat-word implementation realises that genuine lift rather than a projective
twist is verified numerically (homomorphism defect $\sim 10^{-14}$ over all ordered pairs) and is not
proved. The conclusion is insensitive to the gap: two lifts differ by a linear character, and twisting by
a linear character permutes the constituents while preserving dimensions and multiplicity-freeness, so the
shape $1 + 4\times 2$ is lift-independent; only the component labels are not.

Everything the script reads or writes lies inside this repository. There are no external data
dependencies: all inputs are generated from the prime $q$ and the documented seed.

## Environment

Python 3.14.5 and the pinned `requirements.txt` (numpy 2.4.0) were used to produce the deposited
numbers. Numerical reproduction holds over a wider range; the pin exists so that the recorded
checkpoints are byte-reproducible.

```bash
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

## Exact commands

Primary runs (Remark 4.4, sampled blocks):

```bash
python3 o23_filtration_stabiliser.py --q 53  --blocks 3 --n-max 8
python3 o23_filtration_stabiliser.py --q 101 --blocks 3 --n-max 8
python3 o23_filtration_stabiliser.py --q 211 --blocks 2 --n-max 6
```

Controls on the ambient central character, at $q = 53$:

```bash
python3 o23_filtration_stabiliser.py --q 53 --blocks 2 --n-max 6 --central c1
python3 o23_filtration_stabiliser.py --q 53 --blocks 2 --n-max 6 --central one
```

Multiplier sweeps (the $25/2000$, $7/2000$, $1/2000$ figures of Remark 4.4):

```bash
python3 o23_filtration_stabiliser.py --q 53  --blocks 1 --n-max 6 --sweep 2000
python3 o23_filtration_stabiliser.py --q 101 --blocks 1 --n-max 6 --sweep 2000
python3 o23_filtration_stabiliser.py --q 211 --blocks 1 --n-max 6 --sweep 2000
```

Independent sampling stream (the $26/2000$, $4/2000$, $0/2000$ figures). The sweep RNG is seeded at
`seed + 1` and that derived value is recorded in each checkpoint manifest, so re-running with
`--seed 12345` replays the same stream; a genuine resample requires a different `--seed`:

```bash
python3 o23_filtration_stabiliser.py --q 53  --blocks 1 --n-max 6 --sweep 2000 --seed 777
python3 o23_filtration_stabiliser.py --q 101 --blocks 1 --n-max 6 --sweep 2000 --seed 777
python3 o23_filtration_stabiliser.py --q 211 --blocks 1 --n-max 6 --sweep 2000 --seed 777
```

Default seed is 12345 throughout. Runs are deterministic.

## What the script guarantees

- **Licence precheck.** It verifies $\rho_c(g)W_c(v)\rho_c(g)^{-1} = W_c(gv)$ on the exact phase
  (not merely up to proportionality) for the generators and for random Bruhat words, and aborts if
  the identity fails. The Weyl-coefficient reduction used by Theorem 4.2 is not valid otherwise.
- **Two independent tracks.** An integer-only combinatorial track computes $F_n$, $M_n$ and
  $U \rtimes M_n$ with no floating threshold; a floating-point track computes the projectors and
  their stabilisers directly. The script exits with status 3 on any disagreement between them. The
  thresholds `TOL_RANK`, `TOL_SUPP`, `TOL_MEMBER` govern only the second track.
- **No target dimension.** There is no `HEFF_DIM`, no fixed rank, and no branch specific to $2I$ or
  to any binary polyhedral group. Ranks and stabilisers are outputs; the classification is done
  after the fact against Dickson's list.
- **Checkpoint integrity.** Each run writes a JSON checkpoint whose filename and internal manifest
  carry every result-affecting parameter ($q$, seed, `--central`, `--blocks`, `--n-max`,
  `--bfs-frac`, `--max-nodes`, `--chunk`, `--sweep`) plus the derived sweep seed and the thresholds.
  `--resume` exits with status 4 rather than continuing from a checkpoint written under different
  parameters.

## Stored evidence

`checkpoints/` holds the eleven checkpoints backing Remark 4.4: five for the primary and control
runs (84 sampled levels, 168 dimension/stabiliser comparisons, all agreements), three carrying the
seed-12345 sweep histograms, and three carrying the independent seed-777 histograms. Each file
records the per-level dimensions, exact frequency-set sizes, multiplier groups, stabiliser orders,
pointwise stabilisers, projector cross-check norms and the classification.
