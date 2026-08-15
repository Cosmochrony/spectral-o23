This repository contains the source of the **O23 Cosmochrony paper**
*Three Neutral Directions from a Spinor Carrier:
Conditional Status of the Threefold Admissibility Threshold*.

This work belongs to the **spectral admissibility sub-programme** and determines the exact status of the
value three in the threshold condition isolated by **O21** and **O22**:

> A spinorial neutral sector carries exactly three independent directions;
> the threshold $\Sigma_c(n_3) = 3$ is a supplied selection rule compatible with that structure, not a derived constant.

## Context

**O21** established that:

- the physically relevant observable is the canonical fibre-level quantity
  $\sigma_{\mathrm{pair}}^{\mathrm{can}}(n)$
- the admissibility criterion can be reformulated intrinsically through the
  observable rank $n_3^{\mathrm{obs}}$
- the threshold condition
  $\Sigma_c(n_3) = 3$
  defines the physically relevant saturation shell

**O22** proved that:

- saturation necessarily occurs on a BFS shell (projection locking)
- shell-alignment is a theorem, not a conjecture

The value **3** entering the threshold is the object of the present paper.

## Core Results

### 1. Adjoint-dimension theorem (proved, conditional on the carrier)

If the admissible neutral sector is carried by an irreducible two-dimensional $\mathrm{SU}(2)$-valued
representation $V_\rho \cong \mathbb{C}^2$, then its traceless anti-Hermitian sector is the real Lie algebra

$\mathfrak{su}(V_\rho) \cong \mathfrak{su}(2) \cong \mathrm{Im}\,\mathbb{H}$,

of real dimension exactly **3**.
Every neutral generator image $\rho(s) = \mathrm{i}(\vec{u}_s \cdot \vec{\sigma})$ lies on the unit sphere of
this three-dimensional space, for the inner product $\langle A, B\rangle = \tfrac12\mathrm{Tr}(A^\dagger B)$.
This explains why a spinor-carried admissible sector naturally carries three independent directions.

### 2. Non-abelian support (proved)

An irreducible two-dimensional complex representation forces the group to be non-abelian (Schur).
The neutrality condition alone does **not** exclude abelian groups: $\mathbb{Z}/4$ admits a faithful
reducible two-dimensional representation with a neutral symmetric generating set.
The exclusion is carried by the irreducibility of the supplied carrier, not by admissibility.

### 3. Associativity does not select $\mathbb{H}$ (delimitation)

The Weil cocycle structure excludes octonionic realisations (associativity), but the representation
algebra $M_q(\mathbb{C})$ is not a division algebra, so the Hurwitz classification does not apply to it
and cannot force a quaternionic structure.
The quaternions enter only through the supplied $\mathrm{SU}(2)$ carrier.

### 4. Generator axes need not span the three directions (no-go)

In the faithful two-dimensional representation of $Q_8$, the symmetric neutral set
$\{\pm\mathbf{i}, \pm\mathbf{j}\}$ generates the group with only **two** generator axes; the third
direction arises only under commutator (Lie) closure.
Counting admissible generator axes therefore cannot replace the module statement of Result 1.

### 5. The O12 filtration cannot supply the carrier (proved, unconditional)

Every $k=3$ fingerprint vector of the O12 shell-span construction is a **pure Fourier mode**, of frequency
$f = \sum_i c_i b_i \bmod q$ (Lemma 4.1). Each shell span $W_{<n}$ is therefore a coordinate subspace in the
Fourier basis, its projector is Fourier-diagonal, and its Weyl support is a **single line** $L$. Consequently,
for every odd prime $q$, every generic block and every proper nonzero level,

$\mathrm{Stab}(W_{<n}) = U \rtimes M_n \subseteq B(L)$,   $M_n = \{s \in \mathbb{F}_q^\times : sF_n = F_n\}$

where $B(L)$ is the Borel subgroup of $\mathrm{SL}(2,\mathbb{Z}/q\mathbb{Z})$ stabilising $L$ and $U$ its
unipotent radical (Theorem 4.2). Since $B(L)$ is metacyclic and no binary polyhedral group is, **none of
$2T \cong \mathrm{SL}(2,3)$, $2O$, $2I \cong \mathrm{SL}(2,5)$ is selected**, so this filtration supplies no
two-dimensional spinor carrier (Corollary 4.3).

The exclusion is sharpest where the groups are available: by Dickson's classification $2I$ does sit inside
$\mathrm{SL}(2,\mathbb{Z}/q\mathbb{Z})$ whenever $q \equiv \pm 1 \pmod 5$, and the filtration still does not
select it. This closes one **source** of the carrier, not the carrier itself, and is **not** a no-go on the
value three.

Measured separately and kept apart from the theorem (Remark 4.4): on the documented sampled blocks at
$q = 53, 101, 211$ the multiplier groups are $M_n = \{\pm 1\}$, giving $\mathrm{Stab}(W_{<n}) \cong C_{2q}$
with normaliser $B(L)$ — but genericity does **not** force this, and generic blocks with $|M_n| \in \{4,6\}$
at some depth exist. That the intersection over levels is $\{\pm 1\}$ for *all* generic blocks is **not
established**.

## The Two Open Bridges

The conversion of the adjoint-dimension theorem into a derivation of $\Sigma_c(n_3) = 3$ requires two
bridges, both **open**:

1. **Carrier selection**: no result derives $V_\rho \cong \mathbb{C}^2$ from Born–Infeld parity.
   O18 establishes parity as a covariance of the response family only (even-order responses preserved,
   odd-order reversed); the fibre identification is a typed open classification problem (O18 Problem 2.8).
   One candidate source is now **closed exactly**: the O12 shell-span filtration selects no spinor carrier
   (Result 5). The bridge itself remains open — nothing proves that no admissible construction selects one.
   Any replacement must first satisfy the invariant criterion that a proper nonzero projector's Weyl support
   span at least **two distinct lines** (necessary, not sufficient), and enriching the object with shell
   labels or the BFS metric cannot help, since fixing a set of vectors fixes their span.
2. **Observable identification**: no theorem identifies the cumulative Gram–Schmidt span $\Sigma_c$ with
   $\dim_{\mathbb{R}} \mathfrak{su}(V_\rho)$.
   Since $n_3$ is defined as the shell where $\Sigma_c$ reaches three, the threshold is definitional at
   the level of O21.

The status-typed chain is:

BI parity (covariant, proved) →(open)→ pair fibre $c \leftrightarrow q{-}c$ →(open; **not from the O12
filtration**, Result 5)→ $V_\rho \cong \mathbb{C}^2$ →(**proved**)→
$\dim_{\mathbb{R}} \mathfrak{su}(V_\rho) = 3$ →(open)→ $\Sigma_c(n_3) = 3$.

## Spectral Observations

- $Q_8$ realises the three directions isotropically through its **full** neutral sector
  ($M \propto I$), while its two-axis generating subset shows a generating set alone need not.
- The ADE binary graphs ($2I$ in particular) show exactly three non-trivial Cayley eigenvalue classes.
  This is retained as a **consistency observation**, compatible with — but not probative of — the
  threefold structure.

## Status

- Adjoint dimension of a supplied spinor carrier: **proved** (Theorem 3.1)
- Exact stabiliser of the O12 filtration, $U \rtimes M_n \subseteq B(L)$: **proved, unconditional** (Theorem 4.2)
- No binary polyhedral group selected by that filtration: **proved, unconditional** (Corollary 4.3)
- $\mathrm{Stab}(W_{<n}) \cong C_{2q}$: **measured on the sampled blocks**, not a theorem (Remark 4.4)
- Filtration intersection $\{\pm 1\}$ for all generic blocks: **not established** (Remark 4.4)
- $\Sigma_c(n_3) = 3$: **supplied selection rule / open**, awaiting the two bridges
- Location of $n_3$ as a function of pipeline parameters: open (O22 §6.2 programme)

## Reproduction

The numerical statements of Result 5 are reproduced by the code in `code/`, which is self-contained: it reads
and writes only inside this repository and needs no external data.

```bash
cd code
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
python3 o23_filtration_stabiliser.py --q 53 --blocks 3 --n-max 8
```

See [`code/README.md`](code/README.md) for the complete list of commands, the thresholds, the seeds, and the
stored checkpoint evidence. The script carries two independent tracks — an integer-only combinatorial one with
no floating threshold, and a floating-point one — and aborts if they disagree. It contains no `HEFF_DIM`, no
fixed target rank, and no branch specific to any binary polyhedral group.

## Repository Structure

```text
o23/
├── code/     # Reproduction code for Result 5, with stored checkpoint evidence
├── out/      # Compiled O23 PDF (generated by compile.sh, git-ignored)
├── tex/      # LaTeX sources
└── README.md
```

## Citation

If you reference this work, please cite via the Zenodo concept DOI:

J. Beau,
*Three Neutral Directions from a Spinor Carrier: Conditional Status of the Threefold Admissibility Threshold*,
Zenodo, 2026. [doi:10.5281/zenodo.19375136](https://doi.org/10.5281/zenodo.19375136)

## Acknowledgements

Portions of the derivations, conceptual synthesis, and editorial refinement benefited from iterative
interactions with large language models used as analytical assistants.
All theoretical results, computations, and interpretations remain the sole responsibility of the author.
