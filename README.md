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

### 5. What the O12 filtration does and does not supply (proved / computational)

Every $k=3$ fingerprint vector of the O12 shell-span construction is a **pure Fourier mode**, of frequency
$f = \sum_i c_i b_i \bmod q$. Each **level** $W_{<n}$ is therefore a coordinate subspace in the Fourier
basis, its projector is Fourier-diagonal, and its Weyl support is a **single line** $L$. Consequently, for
every odd prime $q$, every generic block and every proper nonzero level,

$\mathrm{Stab}(W_{<n}) = U \rtimes M_n \subseteq B(L)$,   $M_n = \{s \in \mathbb{F}_q^\times : sF_n = F_n\}$

where $B(L)$ is the Borel subgroup stabilising $L$ and $U$ its unipotent radical. This is a statement about
one **level**; the full-filtration stabiliser $\bigcap_n \mathrm{Stab}(W_{<n})$ is a different, smaller
object.

**Excluded (proved, unconditional):** the **exceptional** binary polyhedral groups
$2T \cong \mathrm{SL}(2,3)$, $2O$, $2I \cong \mathrm{SL}(2,5)$, since $B(L)$ is metacyclic and none of them
is. This holds even at the primes where they exist in the ambient group, e.g. $2I$ for $q \equiv \pm1 \pmod 5$.

**Not excluded, and realised:** dicyclic (binary dihedral) groups *are* metacyclic and do carry faithful
two-dimensional spinor representations. At the explicit generic levels $(q,c,n) = (53,(47,21,32),1)$ and
$(101,(41,95,6),1)$ the stabiliser is $\mathrm{Dic}_q$ of order $4q$, and the level decomposes
multiplicity-freely under it as

$W_{<1} \cong \mathbf{1}' \oplus V_1 \oplus V_2 \oplus V_3 \oplus V_4$,   $\dim_{\mathbb{C}} V_i = 2$,   $1 + 4\times 2 = 9$

with $\langle\chi,\chi\rangle = 5$ and every multiplicity equal to one. This is a **computational result on
two named witnesses**, not a universal claim about dicyclic levels.

Exactly **two** of the four are admissible carriers. $V_j$ is faithful, hence $\mathrm{SU}(2)$-valued as
Result 1 requires, iff $j$ is odd, because the central element $a^q$ acts by $(-1)^j$; for even $j$ the
representation factors through the dihedral quotient and is not a spinor carrier. The admissible indices are
$j \in \{11,47\}$ at $q=53$ and $j \in \{41,89\}$ at $q=101$.

**Consequence.** The filtration supplies a spinorial group and two admissible carriers. Each component is
mathematically canonical — the decomposition is multiplicity-free — so what is missing is not canonicity but
a **structurally justified** reason to prefer one admissible carrier over the other. Bridge 1 is therefore
**displaced, not closed**, with a **twofold** residual ambiguity. No claim is made that such a selector
cannot exist.

One internal candidate is recorded together with its refutation: $c_\Sigma$ is a constituent in every
dicyclic level examined, and is an admissible carrier in the two deposited witnesses, but at $q=53$ with
block $(10,35,18)$ it is even, hence not faithful. The rule "$j = c_\Sigma$" therefore does not select an
admissible carrier.

## The Two Open Bridges

The conversion of the adjoint-dimension theorem into a derivation of $\Sigma_c(n_3) = 3$ requires two
bridges, both **open**:

1. **Carrier selection**: no result derives $V_\rho \cong \mathbb{C}^2$ from Born–Infeld parity.
   O18 establishes parity as a covariance of the response family only (even-order responses preserved,
   odd-order reversed); the fibre identification is a typed open classification problem (O18 Problem 2.8).
   One candidate source is now analysed exactly (Result 5): it excludes the exceptional binary polyhedral
   groups, but it *does* select a dicyclic stabiliser carrying multiplicity-free two-dimensional content at
   explicit levels. The bridge is therefore sharpened rather than closed — what is missing is a justified
   selector between the two admissible carriers. Enriching the object with shell labels or the BFS metric
   cannot enlarge the stabiliser, since fixing a set of vectors fixes their span.
2. **Observable identification**: no theorem identifies the cumulative Gram–Schmidt span $\Sigma_c$ with
   $\dim_{\mathbb{R}} \mathfrak{su}(V_\rho)$.
   Since $n_3$ is defined as the shell where $\Sigma_c$ reaches three, the threshold is definitional at
   the level of O21.

The status-typed chain is:

BI parity (covariant, proved) →(open)→ pair fibre $c \leftrightarrow q{-}c$ →(open; refined by Result 5 into an
existence part supplied and a selection part missing)→ $V_\rho \cong \mathbb{C}^2$ →(**proved**)→ $\dim_{\mathbb{R}} \mathfrak{su}(V_\rho) = 3$ →(open)→
$\Sigma_c(n_3) = 3$.

## Spectral Observations

- $Q_8$ realises the three directions isotropically through its **full** neutral sector
  ($M \propto I$), while its two-axis generating subset shows a generating set alone need not.
- The ADE binary graphs ($2I$ in particular) show exactly three non-trivial Cayley eigenvalue classes.
  This is retained as a **consistency observation**, compatible with — but not probative of — the
  threefold structure.

## Status

- Adjoint dimension of a supplied spinor carrier: **proved** (Theorem 3.1)
- Exact stabiliser of each O12 **level**, $U \rtimes M_n \subseteq B(L)$: **proved, unconditional**
- No **exceptional** binary polyhedral group ($2T$, $2O$, $2I$) selected: **proved, unconditional**
- Dicyclic stabiliser with multiplicity-free spinorial content: **computational, two named witnesses**
- Exactly two of the four constituents faithful (admissible carriers): **proved from the central character**
- A structurally justified selector between the two admissible carriers: **open** (the candidate
  "$j = c_\Sigma$" is refuted by an explicit counterexample)
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
