#!/usr/bin/env python3
"""Exact stabiliser of the deposited O12 filtration inside the Weil image.

Reproduction code for O23 v2.1, Section 4.4: Lemma 4.1 (pure Fourier modes), Theorem 4.2
(Stab(W_<n) = U >< M_n, contained in the Borel B(L)), Corollary 4.3 (no binary polyhedral
group is selected) and Remark 4.4 (measured multiplier groups).  Self-contained: every path
it reads or writes lies inside this repository.

Question. The O12 construction produces a filtration W_{<n} of subspaces of C^q.  The Weil
(metaplectic) representation rho_c of SL(2, Z/q) acts on the same C^q.  Which subgroup of
SL(2, Z/q), if any, stabilises W_{<n}?  The point of asking in this order is that a subgroup
found here is SELECTED by the deposited filtration, whereas a subgroup chosen by hand (2I, 2O,
2T, ...) would be SUPPLIED -- the HEFF_DIM = 3 circularity moved one level up.

No target dimension, no HEFF_DIM, no 2I-specific branch appears anywhere below.  The
classification of whatever is found is a posteriori (lock 4).

Method.  Membership is the invariant projector criterion (lock 3)

    g in Stab(W_{<n})   <=>   rho_c(g) P_n rho_c(g)^dag = P_n .

Executed through its exact Weyl-coefficient form.  Expanding P_n in the Weyl operator basis,
P_n = (1/q) sum_v p_n(v) W_c(v), and using rho_c(g) W_c(v) rho_c(g)^{-1} = W_c(g v), the
criterion is equivalent to

    p_n o g^{-1} = p_n   on F_q^2 ,

which costs O(q^2) per group element instead of O(q^3), and stays basis-free.  The intertwining
identity is VERIFIED on the generators before use; the script aborts if it fails, because the
reduction is not licensed otherwise.

Group elements are never enumerated blindly.  Any g in Stab permutes the level sets of p_n, so
g maps a chosen v0 in the rarest level class C into C.  The elements sending v0 to a given w form
a coset of the order-q unipotent stabiliser of v0, giving |C| * q candidates -- orbit-stabiliser,
as required, with the candidate count and ETA announced before the enumeration starts.

Usage:
    python3 front_filtration_stabiliser.py --q 53
    python3 front_filtration_stabiliser.py --q 101 --blocks 8
    python3 front_filtration_stabiliser.py --q 211 --blocks 4 --resume

Author-side calibration: every number this prints is a measurement of the DEPOSITED filtration,
not of the programme's intended carrier.  Nothing here bridges to Sigma_c(n_3) = 3; that bridge
is explicitly open (preregistration section 4).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

# Preregistered thresholds (section 5 of the preregistration; do not tune after seeing output).
TOL_RANK = 1e-9          # singular-value cut, relative to sigma_1
MIN_GAP = 1e2            # spectral-gap factor below which a rank is flagged unstable
TOL_SUPP = 1e-8          # Weyl-coefficient support threshold, relative to |p(0)|
TOL_MEMBER = 1e-8        # stabiliser membership, relative to |p(0)|
TOL_VERIFY = 1e-8        # absolute tolerance for the intertwining precheck
PRESCREEN = 64           # points used in the cheap first-stage membership filter
NORMALISER_BUDGET = 400  # max |Stab| for which the normaliser is enumerated exactly

CHECKPOINT_DIR = Path(__file__).resolve().parent / "checkpoints"


# Field helpers

def inv_mod(a, q):
    return pow(int(a) % q, q - 2, q)


def e_q(t, q):
    return np.exp(2j * np.pi * (np.asarray(t) % q) / q)


def primitive_root(q):
    for a in range(2, q):
        seen, v = set(), 1
        for _ in range(q - 1):
            v = (v * a) % q
            seen.add(v)
        if len(seen) == q - 1:
            return a
    raise ValueError(f"no primitive root found for q={q}")


# Deposited O12 construction (faithful reimplementation of spectral_O12.py)

def heisenberg_mul(u, v, q):
    a, b, gam = u
    ap, bp, gamp = v
    return ((a + ap) % q, (b + bp) % q, (gam + gamp + a * bp) % q)


def heisenberg_inv(u, q):
    a, b, gam = u
    return ((-a) % q, (-b) % q, (a * b - gam) % q)


def build_generators(q):
    """Standard symmetric generating set S = {+/-X, +/-Y}."""
    X, Y = (1, 0, 0), (0, 1, 0)
    return [X, heisenberg_inv(X, q), Y, heisenberg_inv(Y, q)]


def bfs_shells(gens, q, max_nodes):
    """BFS from the identity, capped at max_nodes visited elements."""
    identity = (0, 0, 0)
    visited = {identity}
    current = [identity]
    shells = [current]
    total = 1
    while total < max_nodes:
        nxt = []
        for u in current:
            for g in gens:
                v = heisenberg_mul(u, g, q)
                if v not in visited:
                    visited.add(v)
                    nxt.append(v)
                    total += 1
                    if total >= max_nodes:
                        break
            if total >= max_nodes:
                break
        if not nxt:
            break
        shells.append(nxt)
        current = nxt
    return shells


def heisenberg_mul_batch(u_arr, g, q):
    ap, bp, gamp = int(g[0]), int(g[1]), int(g[2])
    out = u_arr.copy()
    out[:, 2] = (u_arr[:, 2] + gamp + u_arr[:, 0] * bp) % q
    out[:, 0] = (u_arr[:, 0] + ap) % q
    out[:, 1] = (u_arr[:, 1] + bp) % q
    return out


def weil_batch_lut(a_arr, b_arr, gamma_arr, psi_table, q):
    """rho_c(a,b,gamma)|uniform> for a batch, exactly as in spectral_O12.py."""
    x_out = np.arange(q, dtype=np.int64)
    x_in = (x_out[None, :] - a_arr[:, None]) % q
    arg = (gamma_arr[:, None] + b_arr[:, None] * x_in) % q
    return psi_table[arg] / np.sqrt(q)


def fingerprint_vectors_batch(shell_arr, c_block, gens_arr, q):
    """All k=3 fingerprint vectors of a shell: (M * 4^3, q) complex."""
    psis = [np.exp(2j * np.pi * int(c) * np.arange(q, dtype=np.int64) / q) for c in c_block]
    out = []
    for s1 in gens_arr:
        ep1 = heisenberg_mul_batch(shell_arr, s1, q)
        v1 = weil_batch_lut(ep1[:, 0], ep1[:, 1], ep1[:, 2], psis[0], q)
        for s2 in gens_arr:
            ep2 = heisenberg_mul_batch(ep1, s2, q)
            v2 = weil_batch_lut(ep2[:, 0], ep2[:, 1], ep2[:, 2], psis[1], q)
            for s3 in gens_arr:
                ep3 = heisenberg_mul_batch(ep2, s3, q)
                v3 = weil_batch_lut(ep3[:, 0], ep3[:, 1], ep3[:, 2], psis[2], q)
                out.append(v1 * v2 * v3)
    return np.concatenate(out, axis=0)


def sample_generic_blocks(q, m_block, rng):
    """Deposited genericity: all c_i != 0 and sum c_i != 0 (mod q)."""
    blocks = []
    while len(blocks) < m_block:
        c = rng.integers(1, q, size=3)
        if (int(c[0]) + int(c[1]) + int(c[2])) % q != 0:
            blocks.append(c)
    return np.array(blocks, dtype=np.int64)


# Projector onto W_{<n}, with the declared rank criterion

def projector_from_vectors(vecs, q):
    """Orthogonal projector onto span(vecs) plus the rank and its spectral gap.

    vecs: (N, q) complex, treated as spanning vectors.  No basis is privileged and no
    dimension is fixed: the rank is an output (lock 3).
    """
    A = np.asarray(vecs).T                       # (q, N): spanning vectors as columns
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    if s.size == 0 or s[0] == 0:
        return np.zeros((q, q), dtype=complex), 0, np.inf
    keep = s / s[0] >= TOL_RANK
    rank = int(keep.sum())
    if rank == 0:
        return np.zeros((q, q), dtype=complex), 0, np.inf
    if rank < s.size and s[rank] > 0:
        gap = float(s[rank - 1] / s[rank])
    else:
        gap = float("inf")
    Ur = U[:, :rank]
    return Ur @ Ur.conj().T, rank, gap


# Weyl operators and the coefficient function p(v)

def weyl_op(q, c, a, b):
    """Symmetrised Weyl operator W_c(a,b): f(x) -> e_q(c b (x - a/2)) f(x - a)."""
    inv2 = inv_mod(2, q)
    mat = np.zeros((q, q), dtype=complex)
    x = np.arange(q)
    phase = (b * ((x - a * inv2) % q)) % q
    mat[x, (x - a) % q] = e_q(c * phase, q)
    return mat


def weyl_coefficients(P, q, c):
    """p(a,b) = tr(P W_c(a,b)^dag) for every v = (a,b), in O(q^2 log q).

    tr(P W^dag) = e_q(c b a/2) * sum_x P[x, x-a] e_q(-c b x), i.e. a phase times the DFT of the
    a-th generalised diagonal of P, sampled at frequency c*b.
    """
    inv2 = inv_mod(2, q)
    x = np.arange(q)
    diags = np.empty((q, q), dtype=complex)
    for a in range(q):
        diags[a] = P[x, (x - a) % q]
    dft = np.fft.fft(diags, axis=1)              # dft[a, k] = sum_x diags[a,x] e^{-2pi i k x/q}
    b = np.arange(q)
    k = (c * b) % q
    p = dft[:, k]                                # (a, b)
    p *= e_q(c * (b[None, :] * ((np.arange(q)[:, None] * inv2) % q)) % q, q)
    return p


# Weil representation of SL(2, Z/q), built from Bruhat words

class WeilRep:
    """rho_c on C^q, with the symplectic action of each generator determined numerically.

    Generators (conventions of front_equivariance_killswitch.py::metaplectic_generators):
        L(u): diag(e_q(c * (u/2) * x^2))   quadratic phase
        Dg(s): f(x) -> f(s^{-1} x)          scaling
        S:    finite Fourier at character c
    Their symplectic matrices are not hardcoded: they are read off from conjugation on the Weyl
    basis and then checked against SL(2,q).
    """

    def __init__(self, q, c):
        self.q, self.c = q, c
        self.inv2 = inv_mod(2, q)
        x = np.arange(q)
        self.S = e_q(c * np.outer(x, x), q) / np.sqrt(q)
        self._sympl_cache = {}

    def L(self, u):
        """Lower unipotent [[1,0],[u,1]]."""
        q, c = self.q, self.c
        x = np.arange(q)
        t = (u * self.inv2) % q
        return np.diag(e_q(c * t * x * x, q))

    def Dg(self, s):
        """Split torus [[s,0],[0,s^{-1}]]."""
        q = self.q
        sinv = inv_mod(s, q)
        x = np.arange(q)
        mat = np.zeros((q, q), dtype=complex)
        mat[x, (sinv * x) % q] = 1.0
        return mat

    def U(self, t):
        """Upper unipotent [[1,t],[0,1]] = S L(-t) S^{-1}."""
        return self.S @ self.L((-t) % self.q) @ self.S.conj().T

    def rho(self, g):
        """rho_c(g) for g = [[al,be],[ga,de]] in SL(2,q), via Bruhat.

        ga == 0:  g = Dg(al) U(be/al)
        ga != 0:  g = U(al/ga) S Dg(ga) U(de/ga)
        Defined up to the generator phases; the projective ambiguity is reported, never hidden.
        """
        q = self.q
        al, be, ga, de = (int(v) % q for v in g)
        if ga == 0:
            return self.Dg(al) @ self.U((be * inv_mod(al, q)) % q)
        gainv = inv_mod(ga, q)
        return (self.U((al * gainv) % q) @ self.S @ self.Dg(ga)
                @ self.U((de * gainv) % q))

    def symplectic_of(self, M):
        """Read off the 2x2 action v -> g v from conjugation M W(v) M^dag = phase * W(g v).

        Returns the matrix as a tuple, or None if the conjugation does not permute the Weyl
        basis (which would invalidate the whole reduction).
        """
        q, c = self.q, self.c
        img = {}
        for v in [(1, 0), (0, 1)]:
            W = weyl_op(q, c, v[0], v[1])
            Wc = M @ W @ M.conj().T
            # locate the Weyl operator proportional to Wc via the coefficient map
            p = weyl_coefficients(Wc / q, q, c)   # tr(Wc W^dag)/q peaks on the image point
            idx = np.unravel_index(np.argmax(np.abs(p)), p.shape)
            mag = np.abs(p[idx])
            resid = np.abs(p).copy()
            resid[idx] = 0.0
            if mag < 0.5 or resid.max() > 1e-8 * max(mag, 1.0):
                return None
            img[v] = (int(idx[0]), int(idx[1]))
        (a11, a21) = img[(1, 0)]
        (a12, a22) = img[(0, 1)]
        return (a11 % q, a12 % q, a21 % q, a22 % q)


def sl2_mul(g, h, q):
    a, b, c_, d = g
    e, f, g_, h_ = h
    return ((a * e + b * g_) % q, (a * f + b * h_) % q,
            (c_ * e + d * g_) % q, (c_ * f + d * h_) % q)


def sl2_inv(g, q):
    a, b, c_, d = g
    return (d % q, (-b) % q, (-c_) % q, a % q)


def verify_intertwining(rep, sample_v, rng):
    """Mandatory precheck: rho(g) W(v) rho(g)^{-1} = W(g v), and Bruhat words reproduce
    the generators.  Returns (ok, message)."""
    q, c = rep.q, rep.c
    s = primitive_root(q)
    named = {
        "L(1)":  (rep.L(1),  (1, 0, 1, 1)),
        "Dg(s)": (rep.Dg(s), (s % q, 0, 0, inv_mod(s, q))),
        "S":     (rep.S,     (0, q - 1, 1, 0)),
    }
    for name, (M, expected) in named.items():
        got = rep.symplectic_of(M)
        if got is None:
            return False, f"{name}: conjugation does not permute the Weyl basis"
        a, b, cc, d = got
        if (a * d - b * cc) % q != 1:
            return False, f"{name}: symplectic matrix {got} has determinant != 1"
        if got != tuple(int(v) % q for v in expected):
            return False, f"{name}: symplectic action {got} != expected {expected}"
    # end-to-end: Bruhat-built rho(g) must conjugate W(v) to W(g v) exactly
    for _ in range(6):
        g = random_sl2(q, rng)
        R = rep.rho(g)
        if np.abs(R @ R.conj().T - np.eye(q)).max() > 1e-8:
            return False, f"rho({g}) is not unitary"
        for v in sample_v:
            W = weyl_op(q, c, v[0], v[1])
            lhs = R @ W @ R.conj().T
            gv = ((g[0] * v[0] + g[1] * v[1]) % q, (g[2] * v[0] + g[3] * v[1]) % q)
            rhs = weyl_op(q, c, gv[0], gv[1])
            ov = np.vdot(rhs, lhs) / q
            # the reduction p o g^{-1} = p needs ov = 1 exactly, NOT merely |ov| = 1:
            # a residual phase would rescale the Weyl coefficients and break the equivalence.
            if abs(ov - 1.0) > TOL_VERIFY:
                return False, (f"rho({g}) W({v}) rho^-1 = ({ov:.6f}) W({gv}); the reduction "
                               f"requires the phase to be exactly 1")
    return True, "exact intertwining (phase = 1) verified on generators and Bruhat words"


def random_sl2(q, rng):
    while True:
        a, b, c_ = (int(rng.integers(0, q)) for _ in range(3))
        if a % q != 0:
            d = ((1 + b * c_) * inv_mod(a, q)) % q
            return (a % q, b % q, c_ % q, d)


# Stabiliser: orbit-stabiliser over the rarest level class of p

def value_key(z):
    """Hash a complex value at the declared membership tolerance."""
    return (round(float(z.real) / TOL_MEMBER), round(float(z.imag) / TOL_MEMBER))


def unipotent_stabiliser(v0, q):
    """The order-q subgroup fixing v0 (v0 != 0), as a list of SL(2,q) tuples."""
    a, b = int(v0[0]) % q, int(v0[1]) % q
    # complete v0 to a basis B = [v0 v1] with det 1
    if a != 0:
        v1 = (0, inv_mod(a, q))
    else:
        v1 = ((-inv_mod(b, q)) % q, 0)
    B = (a, v1[0], b, v1[1])
    assert (B[0] * B[3] - B[1] * B[2]) % q == 1
    Binv = sl2_inv(B, q)
    return [sl2_mul(sl2_mul(B, (1, t, 0, 1), q), Binv, q) for t in range(q)]


def map_v0_to_w(v0, w, q):
    """One g in SL(2,q) with g v0 = w (both nonzero)."""
    def basis(v):
        a, b = int(v[0]) % q, int(v[1]) % q
        if a != 0:
            v1 = (0, inv_mod(a, q))
        else:
            v1 = ((-inv_mod(b, q)) % q, 0)
        return (a, v1[0], b, v1[1])
    return sl2_mul(basis(w), sl2_inv(basis(v0), q), q)


_WORKER = {}


def _worker_init(p_flat, q, prescreen_idx, tol_abs):
    _WORKER["p"] = np.asarray(p_flat).reshape(q, q)
    _WORKER["q"] = q
    _WORKER["pre"] = prescreen_idx
    _WORKER["tol"] = tol_abs
    aa, bb = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    _WORKER["aa"] = aa
    _WORKER["bb"] = bb


def _test_chunk(chunk):
    """Two-stage membership test for a chunk of candidate group elements."""
    p, q, pre, tol = _WORKER["p"], _WORKER["q"], _WORKER["pre"], _WORKER["tol"]
    aa, bb = _WORKER["aa"], _WORKER["bb"]
    pre_a, pre_b = pre
    p_pre = p[pre_a, pre_b]
    out = []
    for g in chunk:
        gi = sl2_inv(tuple(g), q)
        a2 = (gi[0] * pre_a + gi[1] * pre_b) % q
        b2 = (gi[2] * pre_a + gi[3] * pre_b) % q
        if np.abs(p[a2, b2] - p_pre).max() > tol:
            continue
        A2 = (gi[0] * aa + gi[1] * bb) % q
        B2 = (gi[2] * aa + gi[3] * bb) % q
        if np.abs(p[A2, B2] - p).max() <= tol:
            out.append(tuple(int(v) for v in g))
    return out


def stabiliser_of_p(p, q, rng, workers, label="", announce=True):
    """Exact setwise stabiliser {g in SL(2,q) : p o g^{-1} = p}, by orbit-stabiliser."""
    p0 = abs(p[0, 0])
    tol_abs = TOL_MEMBER * max(p0, 1.0)
    mask = np.abs(p) > TOL_SUPP * max(p0, 1.0)
    mask[0, 0] = False
    supp = np.argwhere(mask)
    if supp.size == 0:
        return None, {"support": 0, "note": "p supported at v=0 only: stabiliser is all of SL(2,q)"}

    classes = {}
    for a, b in supp:
        classes.setdefault(value_key(p[a, b]), []).append((int(a), int(b)))
    rarest = min(classes.values(), key=len)
    v0 = rarest[0]

    Uv0 = unipotent_stabiliser(v0, q)
    n_cand = len(rarest) * q
    if announce:
        print(f"    [{label}] |supp p| = {len(supp)}, level classes = {len(classes)}, "
              f"rarest class = {len(rarest)}")
        print(f"    [{label}] candidates = |C| * q = {len(rarest)} * {q} = {n_cand}")

    cands = np.empty((n_cand, 4), dtype=np.int64)
    i = 0
    for w in rarest:
        g0 = map_v0_to_w(v0, w, q)
        for u in Uv0:
            cands[i] = sl2_mul(g0, u, q)
            i += 1

    pre_a = rng.integers(0, q, size=PRESCREEN)
    pre_b = rng.integers(0, q, size=PRESCREEN)
    prescreen = (pre_a, pre_b)

    t0 = time.time()
    probe = min(2000, n_cand)
    _worker_init(p.ravel(), q, prescreen, tol_abs)
    found = _test_chunk(cands[:probe])
    dt = time.time() - t0
    eta = dt / max(probe, 1) * n_cand / max(workers, 1)
    if announce:
        print(f"    [{label}] ETA for the full candidate sweep: {eta:6.1f} s "
              f"on {workers} core(s)  ({n_cand} candidates)")

    if n_cand <= probe:
        stab = found
    else:
        rest = cands[probe:]
        chunks = np.array_split(rest, max(workers * 8, 1))
        stab = list(found)
        done, total = probe, n_cand
        if workers > 1:
            with ProcessPoolExecutor(
                    max_workers=workers, initializer=_worker_init,
                    initargs=(p.ravel(), q, prescreen, tol_abs)) as ex:
                for res in ex.map(_test_chunk, chunks):
                    stab.extend(res)
                    done += len(chunks[0]) if chunks else 0
                    if announce:
                        el = time.time() - t0
                        frac = min(done / total, 1.0)
                        rem = el / max(frac, 1e-9) - el
                        print(f"    [{label}] {100*frac:5.1f}%  elapsed {el:6.1f}s  "
                              f"ETA {rem:6.1f}s", end="\r", flush=True)
        else:
            for ch in chunks:
                stab.extend(_test_chunk(ch))
        if announce:
            print(" " * 78, end="\r")
    stab = sorted(set(tuple(int(v) for v in g) for g in stab))
    return stab, {"support": int(len(supp)), "classes": len(classes),
                  "rarest": len(rarest), "candidates": int(n_cand)}


def q_of(P):
    return P.shape[0]


def pointwise_stabilisers(stab, rep, P, rank):
    """Lock 1, third item: the pointwise and projectively-pointwise stabilisers, MEASURED.

    pointwise:              rho(g)|_W = id          <=>  rho(g) P = P
    projectively pointwise: rho(g)|_W = lambda id   <=>  rho(g) P = lambda P
    The linear version depends on the generator phases of the Bruhat-built rho; the projective
    version does not.  Both are reported, and lambda with them.
    """
    if rank == q_of(P):
        # W = C^q: the pointwise stabilisers are DEFINED, and equal the kernels of rho_c.
        # SL(2,q) is quasi-simple for q > 3, so both are read off lock 2 rather than enumerated.
        ker = weil_projective_kernel(rep)
        n_lin = 1                       # rho_c(-1) = parity != 1, so only the identity acts as id
        n_proj = len(ker["projective_kernel"])
        return {"pointwise": n_lin, "projective_pointwise": n_proj,
                "note": "W = C^q: pointwise stabilisers are the kernels of rho_c (lock 2)"}
    if stab is None or rank == 0:
        return {"pointwise": None, "projective_pointwise": None,
                "note": "rank 0: no subspace to act on"}
    lin, proj, lambdas = [], [], []
    for g in stab:
        R = rep.rho(g)
        RP = R @ P
        num = float(np.linalg.norm(RP - P))
        if num <= TOL_MEMBER * max(rank, 1):
            lin.append(g)
        lam = np.trace(P.conj().T @ RP) / max(rank, 1)
        if float(np.linalg.norm(RP - lam * P)) <= TOL_MEMBER * max(rank, 1):
            proj.append(g)
            lambdas.append(complex(lam))
    return {"pointwise": len(lin), "projective_pointwise": len(proj),
            "projective_phases": sorted({(round(z.real, 6), round(z.imag, 6)) for z in lambdas})}


def weil_projective_kernel(rep):
    """Lock 2: the projective kernel of rho_c, which is NOT the kernel of SL(2,q) -> PSL(2,q).

    SL(2,q) is quasi-simple for q > 3, so any normal subgroup is central; it therefore suffices
    to test the centre element -1.  rho_c(-1) is the parity operator, not a scalar, so the two
    kernels genuinely differ and must be reported separately.
    """
    q = rep.q
    minus_i = (q - 1, 0, 0, q - 1)
    R = rep.rho(minus_i)
    lam = np.trace(R) / q
    scalar_defect = float(np.linalg.norm(R - lam * np.eye(q)))
    x = np.arange(q)
    parity = np.zeros((q, q), dtype=complex)
    parity[x, (-x) % q] = 1.0
    ov = np.vdot(parity, R) / q
    return {
        "centre_of_SL2": ["+1", "-1"],
        "rho(-1)_is_scalar": bool(scalar_defect <= 1e-8),
        "rho(-1)_scalar_defect": scalar_defect,
        "rho(-1)_is_parity_up_to_phase": bool(abs(abs(ov) - 1.0) <= 1e-8),
        "projective_kernel": (["+1", "-1"] if scalar_defect <= 1e-8 else ["+1"]),
    }


def exact_frequency_set(shell_arr, c_block, gens_arr, q):
    """The EXACT set of Fourier frequencies contributed by a shell, by integer arithmetic.

    Each k=3 fingerprint vector is the pointwise product of three Schroedinger vectors whose
    x-dependence is the single phase c_i * b_i; the product therefore carries the single
    frequency sum_i c_i b_i mod q.  No floating point and no threshold enter here.
    """
    c1, c2, c3 = (int(v) for v in c_block)
    out = set()
    for s1 in gens_arr:
        ep1 = heisenberg_mul_batch(shell_arr, s1, q)
        for s2 in gens_arr:
            ep2 = heisenberg_mul_batch(ep1, s2, q)
            for s3 in gens_arr:
                ep3 = heisenberg_mul_batch(ep2, s3, q)
                f = (c1 * ep1[:, 1] + c2 * ep2[:, 1] + c3 * ep3[:, 1]) % q
                out.update(int(v) for v in np.unique(f))
    return out


def multiplier_group(F, q):
    """M = {s in F_q^* : s F = F}, exactly."""
    return sorted(s for s in range(1, q) if {(s * f) % q for f in F} == F)


def multiplier_sweep(q, n_blocks, n_max, shells, gens_arr, rng):
    """How often is M_n = {+-1}?  Integer-only survey over many generic blocks.

    M_n = {+-1} was observed on the documented sample, but genericity (c_i != 0, sum c_i != 0)
    does NOT force it: extra multipliers do occur, typically at shallow depth, and collapse
    later.  This sweep measures the frequency instead of asserting uniformity.  The Borel
    containment is untouched either way -- U >< M_n lies in B(L) for every M_n.
    """
    blocks = sample_generic_blocks(q, n_blocks, rng)
    per_level = [{"levels": 0, "extra": 0, "orders": {}} for _ in range(n_max)]
    filt_orders = {}
    n_blocks_with_extra = 0
    for c_block in blocks:
        F, saw_extra = set(), False
        M_filt = None
        for n in range(1, n_max + 1):
            F |= exact_frequency_set(np.array(shells[n - 1], dtype=np.int64),
                                     c_block, gens_arr, q)
            if len(F) == q:
                break
            M = multiplier_group(F, q)
            rec = per_level[n - 1]
            rec["levels"] += 1
            rec["orders"][len(M)] = rec["orders"].get(len(M), 0) + 1
            if len(M) > 2:
                rec["extra"] += 1
                saw_extra = True
            M_filt = set(M) if M_filt is None else (M_filt & set(M))
        if saw_extra:
            n_blocks_with_extra += 1
        if M_filt is not None:
            filt_orders[len(M_filt)] = filt_orders.get(len(M_filt), 0) + 1
    return {"blocks": int(n_blocks), "blocks_with_extra_multipliers": n_blocks_with_extra,
            "per_level": [{"n": i + 1, **r} for i, r in enumerate(per_level) if r["levels"]],
            "filtration_multiplier_orders": filt_orders}


def predicted_stabiliser(F, q):
    """U >< M, the exact prediction once W = span{e_f : f in F} is a Fourier-coordinate subspace.

    L = <(1,0)> is the Weyl-support line, B(L) = upper triangular matrices.  Its unipotent part
    U = {[[1,t],[0,1]]} acts diagonally in the Fourier basis and so preserves EVERY coordinate
    subspace; its torus element diag(s, s^{-1}) sends e_f to e_{sf} and so preserves W exactly
    when s in M.  Returns None when W is the whole space (every g stabilises it).
    """
    if len(F) == q:
        return None, list(range(1, q))
    M = multiplier_group(F, q)
    U = [(1, t, 0, 1) for t in range(q)]
    T = [(s, 0, 0, inv_mod(s, q)) for s in M]
    return sorted({sl2_mul(u, d, q) for u in U for d in T}), M


def cross_check_projector(stab, rep, P, limit=8):
    """Independent confirmation on the returned stabiliser: ||(1-P) rho(g) P||_2 = 0."""
    q = rep.q
    I = np.eye(q)
    worst = 0.0
    for g in stab[:limit]:
        R = rep.rho(g)
        worst = max(worst, float(np.linalg.norm((I - P) @ R @ P, 2)))
    return worst


# A posteriori classification (lock 4: no 2I-specific branch)

def element_order(g, q):
    e = (1, 0, 0, 1)
    for k in range(1, 2 * q * q):
        e = sl2_mul(e, g, q)
        if e == (1, 0, 0, 1):
            return k
    return -1


def classify_subgroup(stab, q):
    """Order, element-order profile, centre, abelian/solvable flags, Dickson-type reading."""
    if stab is None:
        return {"order": q * (q * q - 1), "type": "all of SL(2,q)"}
    n = len(stab)
    codes = set(stab)
    orders = {}
    for g in stab:
        orders[element_order(g, q)] = orders.get(element_order(g, q), 0) + 1
    abelian = all(sl2_mul(g, h, q) == sl2_mul(h, g, q)
                  for g in stab for h in stab) if n <= 120 else None
    centre = [g for g in stab
              if all(sl2_mul(g, h, q) == sl2_mul(h, g, q) for h in stab)] if n <= 2000 else None
    minus_i = ((q - 1) % q, 0, 0, (q - 1) % q)
    # fixed lines: a subgroup inside a Borel preserves a line of F_q^2
    lines = []
    for v in [(1, t) for t in range(q)] + [(0, 1)]:
        if all(((g[0] * v[0] + g[1] * v[1]) % q) * v[1]
               == ((g[2] * v[0] + g[3] * v[1]) % q) * v[0] % q for g in stab):
            lines.append(v)
    info = {
        "order": n,
        "element_order_profile": dict(sorted(orders.items())),
        "abelian": abelian,
        "centre_order": (len(centre) if centre is not None else None),
        "contains_minus_identity": minus_i in codes,
        "invariant_lines": lines,
        "inside_borel": len(lines) > 0,
    }
    # Dickson-type reading, stated as a possibility list, never as a search target
    cand = []
    if n == 1:
        cand.append("trivial")
    if n == 2 and minus_i in codes:
        cand.append("centre {+-1}")
    if info["abelian"]:
        cand.append("abelian (cyclic or elementary abelian)")
    if len(lines) >= 1:
        cand.append("contained in a Borel (metacyclic): no binary polyhedral group embeds")
    if n in (24, 48, 120):
        cand.append({24: "order 24 (2T ~ SL(2,3) possible)",
                     48: "order 48 (2O possible)",
                     120: "order 120 (2I ~ SL(2,5) possible)"}[n])
    info["dickson_reading"] = cand or ["unclassified: inspect the profile"]
    return info


def generating_set(stab, q):
    """A small generating set of the subgroup stab (normalising it suffices to normalise stab)."""
    target = set(stab)
    gens, gen = [], {(1, 0, 0, 1)}
    for g in stab:
        if g in gen:
            continue
        gens.append(g)
        frontier = set(gen)
        while frontier:                      # close the current generators under multiplication
            new = set()
            for a in frontier:
                for b in gens:
                    for prod in (sl2_mul(a, b, q), sl2_mul(b, a, q)):
                        if prod not in gen:
                            gen.add(prod)
                            new.add(prod)
            frontier = new
        if len(gen) == len(target):
            break
    assert gen == target, "generating-set closure did not reproduce the stabiliser"
    return gens


def normaliser(stab, q, budget=NORMALISER_BUDGET):
    """Exact normaliser in SL(2,q), by vectorised enumeration (integer arithmetic only)."""
    if stab is None:
        return {"order": q * (q * q - 1), "note": "normaliser of the whole group"}
    gens = generating_set(stab, q)
    if len(gens) > budget:
        return {"order": None,
                "note": f"generating set of size {len(gens)} exceeds the declared budget "
                        f"{budget}; normaliser not enumerated"}
    a = np.arange(q, dtype=np.int64)
    A, C = np.meshgrid(a, a, indexing="ij")
    A, C = A.ravel(), C.ravel()
    keep = ~((A == 0) & (C == 0))
    A, C = A[keep], C[keep]
    elems = []
    # second column (B,D) with A*D - B*C = 1: one particular solution plus t*(A,C)
    Ainv = np.zeros(q, dtype=np.int64)
    for v in range(1, q):
        Ainv[v] = inv_mod(v, q)
    B0 = np.zeros_like(A)
    D0 = np.zeros_like(A)
    nz = A != 0
    D0[nz] = Ainv[A[nz]]
    B0[~nz] = (q - Ainv[C[~nz]]) % q
    for t in range(q):
        elems.append(np.stack([A, (B0 + t * A) % q, C, (D0 + t * C) % q], axis=1))
    G = np.concatenate(elems, axis=0)
    assert ((G[:, 0] * G[:, 3] - G[:, 1] * G[:, 2]) % q == 1).all()
    code = lambda M: (M[:, 0] + q * M[:, 1] + q * q * M[:, 2] + q ** 3 * M[:, 3])
    stab_codes = np.sort(np.array([g[0] + q * g[1] + q * q * g[2] + q ** 3 * g[3]
                                   for g in stab], dtype=np.int64))
    ok = np.ones(G.shape[0], dtype=bool)
    Gi = np.stack([G[:, 3], (q - G[:, 1]) % q, (q - G[:, 2]) % q, G[:, 0]], axis=1)
    for s in gens:
        S = np.array(s, dtype=np.int64)
        m0 = (G[:, 0] * S[0] + G[:, 1] * S[2]) % q
        m1 = (G[:, 0] * S[1] + G[:, 1] * S[3]) % q
        m2 = (G[:, 2] * S[0] + G[:, 3] * S[2]) % q
        m3 = (G[:, 2] * S[1] + G[:, 3] * S[3]) % q
        c0 = (m0 * Gi[:, 0] + m1 * Gi[:, 2]) % q
        c1 = (m0 * Gi[:, 1] + m1 * Gi[:, 3]) % q
        c2 = (m2 * Gi[:, 0] + m3 * Gi[:, 2]) % q
        c3 = (m2 * Gi[:, 1] + m3 * Gi[:, 3]) % q
        codes = c0 + q * c1 + q * q * c2 + q ** 3 * c3
        ok &= np.isin(codes, stab_codes)
        if not ok.any():
            break
    N = G[ok]
    return {"order": int(N.shape[0]),
            "index_over_stab": int(N.shape[0]) // max(len(stab), 1)}


# Driver

def run_q(q, args):
    rng = np.random.default_rng(args.seed)
    CHECKPOINT_DIR.mkdir(exist_ok=True)

    # Checkpoint identity carries EVERY result-affecting parameter, so the central-character
    # controls cannot overwrite the primary run and --resume cannot silently mix conventions.
    manifest = {"q": q, "seed": args.seed, "central": args.central, "blocks": args.blocks,
                "n_max": args.n_max, "bfs_frac": args.bfs_frac, "max_nodes": args.max_nodes,
                "chunk": args.chunk, "sweep": args.sweep,
                "sweep_seed": (args.seed + 1 if args.sweep else None),
                "thresholds": {"rank": TOL_RANK, "support": TOL_SUPP,
                               "member": TOL_MEMBER, "verify": TOL_VERIFY}}
    tag = (f"q{q}_seed{args.seed}_c-{args.central}_b{args.blocks}_n{args.n_max}"
           f"_f{args.bfs_frac}_m{args.max_nodes}_s{args.sweep}")
    ckpt = CHECKPOINT_DIR / f"stab_{tag}.json"
    results = {}
    if args.resume and ckpt.exists():
        stored = json.loads(ckpt.read_text())
        if stored.get("manifest") != manifest:
            print(f"[q={q}] REFUSING to resume: {ckpt.name} was written with a different "
                  f"parameter set.\n  stored: {stored.get('manifest')}\n  current: {manifest}")
            sys.exit(4)
        results = stored.get("blocks", {})
        print(f"[q={q}] resumed {len(results)} completed block(s) from {ckpt.name}")

    workers = args.workers or max(1, (os.cpu_count() or 2) - 1)
    gens = build_generators(q)
    gens_arr = np.array(gens, dtype=np.int64)
    max_nodes = min(int(args.bfs_frac * q ** 3), args.max_nodes)
    print(f"[q={q}] BFS budget {max_nodes} nodes ({args.bfs_frac} * q^3 capped at {args.max_nodes})")
    shells = bfs_shells(gens, q, max_nodes)
    print(f"[q={q}] {len(shells)} BFS shells, sizes {[len(s) for s in shells[:8]]}...")
    n_max = min(args.n_max, len(shells))

    blocks = sample_generic_blocks(q, args.blocks, rng)
    for bi, c_block in enumerate(blocks):
        key = f"block{bi}"
        if key in results:
            print(f"[q={q}] {key} already done, skipping")
            continue
        c_amb = int(c_block.sum()) % q if args.central == "sum" else (
            int(c_block[0]) % q if args.central == "c1" else 1)
        print(f"[q={q}] {key}: c_block = {c_block.tolist()}, ambient character c = {c_amb}")

        rep = WeilRep(q, c_amb)
        ok, msg = verify_intertwining(rep, [(1, 0), (0, 1), (1, 1), (2, 3)], rng)
        print(f"    precheck: {msg}")
        if not ok:
            print(f"[q={q}] ABORT: the Weyl-coefficient reduction is not licensed.")
            sys.exit(2)

        ker = weil_projective_kernel(rep)
        print(f"    kernels (lock 2): centre of SL(2,q) = {{+1,-1}}; "
              f"rho(-1) scalar? {ker['rho(-1)_is_scalar']} "
              f"(defect {ker['rho(-1)_scalar_defect']:.2e}, parity up to phase: "
              f"{ker['rho(-1)_is_parity_up_to_phase']}); "
              f"projective kernel of rho_c = {ker['projective_kernel']}")

        acc = np.zeros((0, q), dtype=complex)
        freqs = set()
        per_level, stab_filt = [], None
        for n in range(1, n_max + 1):
            shell_arr = np.array(shells[n - 1], dtype=np.int64)
            # chunk the shell: a deep shell would otherwise allocate |shell| * 64 * q complexes
            for lo in range(0, shell_arr.shape[0], args.chunk):
                vecs = fingerprint_vectors_batch(shell_arr[lo:lo + args.chunk],
                                                 c_block, gens_arr, q)
                acc = np.concatenate([acc, vecs], axis=0)
                if acc.shape[0] > 4 * q:              # keep the span, bound the memory
                    U, s, _ = np.linalg.svd(acc.T, full_matrices=False)
                    keep = int((s / s[0] >= TOL_RANK).sum()) if s.size and s[0] > 0 else 0
                    acc = U[:, :keep].T.copy()
            P, rank, gap = projector_from_vectors(acc, q)
            p = weyl_coefficients(P, q, c_amb)
            stab, meta = stabiliser_of_p(p, q, rng, workers, label=f"q={q} n<{n}")
            worst = cross_check_projector(stab or [], rep, P)

            # exact combinatorial track: F_n by integer arithmetic, then Stab = U >< M_n
            freqs |= exact_frequency_set(shell_arr, c_block, gens_arr, q)
            pred, M_n = predicted_stabiliser(freqs, q)
            dim_ok = (len(freqs) == rank)
            stab_ok = (pred is None and stab is None) or (
                pred is not None and stab is not None and set(pred) == set(stab))
            pw = pointwise_stabilisers(stab, rep, P, rank)

            flag = "" if gap >= MIN_GAP else "  [UNSTABLE RANK: gap below threshold]"
            print(f"    n<{n}: dim W = {rank}/{q} (exact |F_n| = {len(freqs)}, "
                  f"match {dim_ok}), gap = {gap:.3e}{flag}")
            print(f"          |Stab_set| = {'ALL' if stab is None else len(stab)}, "
                  f"exact U><M_n predicts {'ALL' if pred is None else len(pred)} "
                  f"(M_n = {M_n if len(M_n) <= 6 else str(len(M_n)) + ' elements'}), "
                  f"match {stab_ok}")
            print(f"          pointwise = {pw['pointwise']}, projectively pointwise = "
                  f"{pw['projective_pointwise']}, projector cross-check = {worst:.2e}")
            if not (dim_ok and stab_ok):
                print(f"    [q={q}] ABORT: exact combinatorial track disagrees with the "
                      f"numerical measurement.")
                sys.exit(3)
            per_level.append({
                "n": n, "dim": rank, "exact_dim": len(freqs), "dim_match": dim_ok,
                "gap": (None if math.isinf(gap) else gap),
                "stab_order": (None if stab is None else len(stab)),
                "predicted_order": (None if pred is None else len(pred)),
                "multiplier_group": M_n, "stab_match": stab_ok,
                "pointwise": pw, "meta": meta, "cross_check": worst,
            })
            stab_set = None if stab is None else set(stab)
            stab_filt = stab_set if stab_filt is None else (
                stab_filt if stab_set is None else stab_filt & stab_set)

        stab_list = None if stab_filt is None else sorted(stab_filt)
        cls = classify_subgroup(stab_list, q)
        nrm = normaliser(stab_list, q)
        print(f"    FILTRATION stabiliser: order {cls['order']}, "
              f"reading {cls['dickson_reading']}")
        print(f"    normaliser: {nrm}")
        results[key] = {
            "c_block": c_block.tolist(), "c_ambient": c_amb,
            "kernels": ker, "levels": per_level,
            "filtration_stabiliser_order": cls["order"],
            "classification": cls, "normaliser": nrm,
        }
        ckpt.write_text(json.dumps({"manifest": manifest, "blocks": results},
                                   indent=2, default=str))
        print(f"    checkpoint written to {ckpt}")

    if args.sweep and "_sweep" not in results:
        # The sweep RNG is seeded at seed+1, recorded in the manifest: reproducing the reported
        # rates with `--seed 12345` replays this exact stream, so a genuinely independent check
        # requires a different --seed, not a different call into multiplier_sweep.
        print(f"[q={q}] multiplier sweep over {args.sweep} generic blocks "
              f"(integer arithmetic only, sweep seed {args.seed + 1})")
        sw = multiplier_sweep(q, args.sweep, n_max, shells, gens_arr,
                              np.random.default_rng(args.seed + 1))
        print(f"    blocks showing M_n != {{+-1}} at some proper level: "
              f"{sw['blocks_with_extra_multipliers']}/{sw['blocks']}")
        for r in sw["per_level"]:
            print(f"      n={r['n']}: {r['extra']}/{r['levels']} levels with extra "
                  f"multipliers; |M_n| histogram {r['orders']}")
        print(f"    |M| of the filtration intersection, histogram: "
              f"{sw['filtration_multiplier_orders']}")
        results["_sweep"] = sw
        ckpt.write_text(json.dumps({"manifest": manifest, "blocks": results},
                                   indent=2, default=str))
        print(f"    sweep evidence stored in {ckpt.name}")
    elif args.sweep:
        print(f"[q={q}] sweep already present in {ckpt.name}, not recomputed")

    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, required=True, help="odd prime (53, then 101, then 211)")
    ap.add_argument("--blocks", type=int, default=6, help="generic O12 blocks sampled")
    ap.add_argument("--n-max", type=int, default=10, dest="n_max", help="BFS depths examined")
    ap.add_argument("--bfs-frac", type=float, default=0.30, dest="bfs_frac")
    ap.add_argument("--max-nodes", type=int, default=400_000, dest="max_nodes")
    ap.add_argument("--chunk", type=int, default=1500, help="shell elements per batch")
    ap.add_argument("--central", choices=["sum", "c1", "one"], default="sum",
                    help="ambient central character: sum c_i (primary), c_1 or 1 (controls)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--workers", type=int, default=0, help="0 = cpu_count - 1")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--sweep", type=int, default=0,
                    help="survey M_n over this many extra generic blocks (integer-only)")
    args = ap.parse_args()

    print(f"O23 v2.1 Section 4.4 reproduction  (seed {args.seed})")
    print(f"thresholds: rank {TOL_RANK}, support {TOL_SUPP}, membership {TOL_MEMBER}")
    t0 = time.time()
    run_q(args.q, args)
    print(f"[q={args.q}] total {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
