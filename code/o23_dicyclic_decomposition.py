#!/usr/bin/env python3
"""Decomposition of an O12 filtration level under its selected dicyclic stabiliser.

Reproduction code for O23 v2.1, subsection "A proved obstruction on Bridge 1":
Proposition "Two exact dicyclic witnesses with multiplicity-free spinorial content".

Context.  o23_filtration_stabiliser.py proves that every proper nonzero O12 level has exact
stabiliser U >< M_n inside a Borel B(L).  Containment in a metacyclic Borel excludes the
EXCEPTIONAL binary polyhedral groups 2T, 2O, 2I, but NOT the dicyclic groups Dic_q, which are
metacyclic and do carry faithful two-dimensional spinor representations.  This script exhibits
levels whose stabiliser is dicyclic and decomposes the level as a representation of it.

What is computed, for each witness (q, block c, level n):
  1. the exact frequency set F_n, the multiplier group M_n and Stab = U >< M_n (integer only);
  2. verification that Stab is dicyclic of order 4q: an element a of order 2q, an element x of
     order 4 with x^2 = -I and x^{-1} a x = a^{-1};
  3. the invariance defect ||(1-P) rho(g) P|| over the whole stabiliser;
  4. the homomorphism defect of the Bruhat-word lift over all ordered pairs of the stabiliser;
  5. the character of W_{<n}, its norm <chi,chi>, and the multiplicities of the one- and
     two-dimensional irreducible representations of Dic_q;
  6. the isotypic projectors, checked to be idempotent, to sum to the identity on W_{<n}, and to
     have rank d_i * m_i.

Epistemic status of the lift.  A genuine linear Weil representation of SL(2,q) exists for q an
odd prime (Weil 1964; Gerardin 1977).  That THIS Bruhat-word implementation realises that
genuine lift, rather than a projective twist of it, is verified numerically here -- the
homomorphism defect vanishes to ~1e-14 over ALL |G|^2 ordered pairs of the stabiliser, which is
exhaustive on the group that matters -- and is NOT proved; the character computation is
conditional on it.  A tempting robustness argument does not close this gap and must not be
quoted as if it did: two maps that are BOTH homomorphisms lifting the same projective
representation differ by a linear character, so twisting shows only that the component LABELS
are convention-dependent once the map is known to be a homomorphism.  It says nothing at all if
the constructed map fails to be one.

Usage:
    python3 o23_dicyclic_decomposition.py                 # both deposited witnesses
    python3 o23_dicyclic_decomposition.py --q 53 --block 47,21,32 --level 1
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from o23_filtration_stabiliser import (
    bfs_shells, build_generators, element_order, exact_frequency_set, inv_mod,
    multiplier_group, predicted_stabiliser, sl2_inv, sl2_mul, WeilRep,
)

# Declared thresholds (reported with every run; the scientific claims are integer-valued)
TOL_INVARIANCE = 1e-8     # ||(1-P) rho(g) P||
TOL_HOM = 1e-8            # ||rho(g)rho(h) - rho(gh)||
TOL_MULT = 1e-6           # distance of a multiplicity from the nearest integer
TOL_PROJ = 1e-8           # idempotence / completeness of the isotypic projectors

# The two witnesses deposited in O23 v2.1 (dicyclic-witnesses proposition), each with the
# EXACT published shape.  The kill-switch below compares against these, so a future change that
# altered the result would fail the run instead of silently rewriting the evidence.
# faithful_j: V_j is faithful (an admissible SU(2)-valued carrier) iff j is odd, since the
# central element a^q acts by (-1)^j.
WITNESSES = [(53, (47, 21, 32), 1), (101, (41, 95, 6), 1), (53, (10, 35, 18), 1)]
PUBLISHED = {
    (53, (47, 21, 32), 1): {"dim_W": 9, "stabiliser_order": 212, "chi_norm": 5.0,
                            "one_dim": 1, "two_dim": 4, "two_dim_j": [6, 11, 42, 47],
                            "faithful_j": [11, 47]},
    (101, (41, 95, 6), 1): {"dim_W": 9, "stabiliser_order": 404, "chi_norm": 5.0,
                            "one_dim": 1, "two_dim": 4, "two_dim_j": [12, 41, 60, 89],
                            "faithful_j": [41, 89]},
    # Cited in the paper as the counterexample to the ORIENTED rule j = c_Sigma: here c_Sigma
    # is a constituent but is even, hence not faithful.  The SIGN-COMPLETED rule
    # j in {c_Sigma, q - c_Sigma} then take the odd one still selects correctly (43).
    (53, (10, 35, 18), 1): {"dim_W": 9, "stabiliser_order": 212, "chi_norm": 5.0,
                            "one_dim": 1, "two_dim": 4, "two_dim_j": [10, 17, 36, 43],
                            "faithful_j": [17, 43]},
}

CKPT_DIR = Path(__file__).resolve().parent / "checkpoints"
OUT = CKPT_DIR / "dicyclic_decomposition.json"


def stabiliser_at_level(q, block, level):
    """Exact F_n, M_n and Stab = U >< M_n at the requested level, by integer arithmetic."""
    gens = build_generators(q)
    gens_arr = np.array(gens, dtype=np.int64)
    shells = bfs_shells(gens, q, min(int(0.30 * q ** 3), 400_000))
    freqs = set()
    for n in range(1, level + 1):
        freqs |= exact_frequency_set(np.array(shells[n - 1], dtype=np.int64),
                                     np.array(block), gens_arr, q)
    stab, M = predicted_stabiliser(freqs, q)
    return freqs, M, stab


def identify_dicyclic(stab, M, q):
    """Check the standard dicyclic presentation: |G| = 4q, a of order 2q, x of order 4 with
    x^2 = -I and x^{-1} a x = a^{-1}.  Returns (is_dicyclic, a, x, report)."""
    minus_i = ((q - 1) % q, 0, 0, (q - 1) % q)
    a = sl2_mul(minus_i, (1, 1, 0, 1), q)          # (-I) * unipotent: order 2q
    roots = [s for s in M if (s * s) % q == (q - 1) % q]
    report = {"order": len(stab), "expected_4q": 4 * q, "multiplier_order": len(M),
              "sqrt_minus_one_in_M": roots, "order_a": element_order(a, q)}
    if len(stab) != 4 * q or not roots:
        report["is_dicyclic"] = False
        return False, a, None, report
    s = roots[0]
    x = (s, 0, 0, inv_mod(s, q))
    ok = (element_order(x, q) == 4
          and sl2_mul(x, x, q) == minus_i
          and sl2_mul(sl2_mul(sl2_inv(x, q), a, q), x, q) == sl2_inv(a, q)
          and element_order(a, q) == 2 * q
          and set(stab) == {sl2_mul(pow_sl2(a, i, q), pow_sl2(x, j, q), q)
                            for i in range(2 * q) for j in range(2)})
    report.update({"order_x": element_order(x, q),
                   "x2_is_minus_I": sl2_mul(x, x, q) == minus_i,
                   "x_inverts_a": sl2_mul(sl2_mul(sl2_inv(x, q), a, q), x, q) == sl2_inv(a, q),
                   "generated_by_a_x": ok, "is_dicyclic": bool(ok)})
    return bool(ok), a, x, report


def pow_sl2(g, k, q):
    out = (1, 0, 0, 1)
    for _ in range(k % (2 * 2 * q)):
        out = sl2_mul(out, g, q)
    return out


def fourier_basis(freqs, q):
    """Orthonormal basis of W = span{e_f : f in F}, as columns."""
    x = np.arange(q)
    cols = [np.exp(2j * np.pi * f * x / q) / np.sqrt(q) for f in sorted(freqs)]
    return np.stack(cols, axis=1)


def dicyclic_characters(q, a, x, stab):
    """Character table of Dic_q on the given element list.

    Dic_q = <a, x | a^{2q} = 1, x^2 = a^q, x^{-1} a x = a^{-1}> has 4 linear characters
    (G/[G,G] = C_4, since [G,G] = <a^2> has order q) and q-1 two-dimensional irreducibles
    rho_j, j = 1..q-1, with chi_j(a^m) = 2 cos(pi j m / q) and chi_j = 0 off <a>.
    """
    apow, e = {}, (1, 0, 0, 1)
    for m in range(2 * q):
        apow[e] = m
        e = sl2_mul(e, a, q)
    # linear characters: g = a^m or a^m x, and the class of g in G/[G,G] = C_4 is
    # m mod 2 (from a) plus 1 (from x) when g lies outside <a>
    nu = {}
    for g in stab:
        if g in apow:
            nu[g] = (2 * apow[g]) % 4 if False else (apow[g] % 2) * 2
        else:
            h = sl2_mul(g, sl2_inv(x, q), q)
            nu[g] = ((apow[h] % 2) * 2 + 1) % 4
    lin = [{g: np.exp(2j * np.pi * k * nu[g] / 4) for g in stab} for k in range(4)]
    two = []
    for j in range(1, q):
        two.append({g: (2 * np.cos(np.pi * j * apow[g] / q) if g in apow else 0.0)
                    for g in stab})
    return lin, two, apow


def decompose(q, block, level, verbose=True):
    freqs, M, stab = stabiliser_at_level(q, block, level)
    c_sigma = int(sum(block)) % q
    rep = WeilRep(q, c_sigma)
    is_dic, a, x, dic_report = identify_dicyclic(stab, M, q)

    B = fourier_basis(freqs, q)
    I_q = np.eye(q)
    P = B @ B.conj().T
    rhos = {g: rep.rho(g) for g in stab}
    inv_defect = max(float(np.linalg.norm((I_q - P) @ rhos[g] @ P)) for g in stab)

    hom_defect = 0.0
    for g, h in itertools.product(stab, repeat=2):
        d = float(np.linalg.norm(rhos[g] @ rhos[h] - rhos[sl2_mul(g, h, q)]))
        hom_defect = max(hom_defect, d)

    W = {g: B.conj().T @ rhos[g] @ B for g in stab}          # |F| x |F| blocks
    chi = {g: complex(np.trace(W[g])) for g in stab}
    G = len(stab)
    chi_norm = float(sum(abs(v) ** 2 for v in chi.values()) / G)

    lin, two, _ = dicyclic_characters(q, a, x, stab)

    def multiplicity(table):
        """<chi, chi_i> as a COMPLEX number, with its imaginary part and its distance to the
        nearest integer reported.  A multiplicity is only meaningful if both are negligible;
        rounding first would hide a failure of the character identification."""
        z = complex(sum(chi[g] * np.conj(table[g]) for g in stab) / G)
        return {"value": z, "imag": abs(z.imag),
                "integer_defect": abs(z.real - round(z.real))}

    mult_lin_c = [multiplicity(t) for t in lin]
    mult_two_c = [multiplicity(t) for t in two]
    worst_imag = max([m["imag"] for m in mult_lin_c + mult_two_c])
    worst_int = max([m["integer_defect"] for m in mult_lin_c + mult_two_c])
    mult_lin = [m["value"].real for m in mult_lin_c]
    mult_two = [m["value"].real for m in mult_two_c]

    # isotypic projectors on W: P_i = (d_i/|G|) sum_g conj(chi_i(g)) rho_W(g)
    dim = B.shape[1]
    total = np.zeros((dim, dim), dtype=complex)
    proj_report = []
    for name, tab, d_i, m_i in (
            [(f"lin{k}", lin[k], 1, mult_lin[k]) for k in range(4)]
            + [(f"two{j+1}", two[j], 2, mult_two[j]) for j in range(q - 1)]):
        if abs(m_i) < TOL_MULT:
            continue
        Pi = (d_i / G) * sum(np.conj(tab[g]) * W[g] for g in stab)
        total += Pi
        proj_report.append({
            "irrep": name, "dim": d_i, "multiplicity": round(m_i, 9),
            "idempotence_defect": float(np.linalg.norm(Pi @ Pi - Pi)),
            "rank": int(np.linalg.matrix_rank(Pi, tol=1e-8)),
            "expected_rank": int(round(d_i * m_i)),
        })
    completeness = float(np.linalg.norm(total - np.eye(dim)))

    dim_from_mult = sum(p["dim"] * p["multiplicity"] for p in proj_report)
    two_dim_j = sorted(int(p["irrep"][3:]) for p in proj_report if p["dim"] == 2)
    faithful_j = [j for j in two_dim_j if j % 2 == 1]   # a^q acts by (-1)^j
    result = {
        "q": q, "block": list(block), "level": level, "c_sigma": c_sigma,
        "dim_W": dim, "multiplier_group": M, "stabiliser_order": G,
        "dicyclic": dic_report,
        "invariance_defect": inv_defect,
        "homomorphism_defect": hom_defect,
        "chi_norm": chi_norm,
        "constituents": proj_report,
        "isotypic_completeness_defect": completeness,
        "dimension_check": {"sum_d_times_m": round(dim_from_mult, 9), "dim_W": dim,
                            "ok": abs(dim_from_mult - dim) < TOL_MULT},
        "multiplicity_free": all(abs(p["multiplicity"] - 1) < TOL_MULT for p in proj_report),
        "worst_multiplicity_imaginary_part": worst_imag,
        "worst_multiplicity_integer_defect": worst_int,
        "two_dim_j": two_dim_j,
        "faithful_j": faithful_j,
        "admissible_carrier_count": len(faithful_j),
        "c_sigma_is_a_constituent": c_sigma in two_dim_j,
        "c_sigma_is_an_admissible_carrier": c_sigma in faithful_j,
        # oriented rule "j = c_Sigma" vs sign-completed rule
        # "j in {c_Sigma, q - c_Sigma}, take the odd one"
        "sign_pair": sorted({c_sigma, (q - c_sigma) % q}),
        "sign_pair_are_constituents": set([c_sigma, (q - c_sigma) % q]) <= set(two_dim_j),
        "sign_completed_selection": next(
            (j for j in (c_sigma, (q - c_sigma) % q) if j % 2 == 1), None),
        "sign_completed_rule_selects_an_admissible_carrier": bool(
            next((j for j in (c_sigma, (q - c_sigma) % q) if j % 2 == 1), None) in faithful_j),
        "thresholds": {"invariance": TOL_INVARIANCE, "homomorphism": TOL_HOM,
                       "multiplicity": TOL_MULT, "projector": TOL_PROJ},
    }
    passed = (is_dic
              and inv_defect < TOL_INVARIANCE
              and hom_defect < TOL_HOM
              and completeness < TOL_PROJ
              and result["dimension_check"]["ok"]
              and result["multiplicity_free"]                       # was computed, never tested
              and worst_imag < TOL_MULT                             # multiplicities are real
              and worst_int < TOL_MULT                              # and are integers
              and all(p["idempotence_defect"] < TOL_PROJ for p in proj_report)
              and all(p["rank"] == p["expected_rank"] for p in proj_report))

    # Kill-switch against the exact published shape, for the deposited witnesses only.
    key = (q, tuple(int(v) for v in block), level)
    exp = PUBLISHED.get(key)
    shape_report = None
    if exp is not None:
        n1 = len([p for p in proj_report if p["dim"] == 1])
        n2 = len([p for p in proj_report if p["dim"] == 2])
        shape_report = {
            "dim_W": (dim == exp["dim_W"]),
            "stabiliser_order": (G == exp["stabiliser_order"]),
            "chi_norm": (abs(chi_norm - exp["chi_norm"]) < TOL_MULT),
            "one_dim_count": (n1 == exp["one_dim"]),
            "two_dim_count": (n2 == exp["two_dim"]),
            "two_dim_j": (two_dim_j == exp["two_dim_j"]),
            "faithful_j": (faithful_j == exp["faithful_j"]),
        }
        passed = passed and all(shape_report.values())
    result["published_shape_check"] = shape_report
    result["all_checks_passed"] = bool(passed)

    if verbose:
        n2 = [p for p in proj_report if p["dim"] == 2]
        n1 = [p for p in proj_report if p["dim"] == 1]
        print(f"q={q} block={tuple(block)} level={level}  c_Sigma={c_sigma}")
        print(f"  dim W_<{level} = {dim};  M_n = {M};  |Stab| = {G} "
              f"(= 4q = {4*q}: {G == 4*q});  dicyclic: {is_dic}")
        print(f"  a has order {dic_report['order_a']} = 2q; x has order "
              f"{dic_report.get('order_x')}; x^2 = -I: {dic_report.get('x2_is_minus_I')}; "
              f"x^-1 a x = a^-1: {dic_report.get('x_inverts_a')}")
        print(f"  invariance defect {inv_defect:.2e}; homomorphism defect over all "
              f"{G*G} pairs {hom_defect:.2e}")
        print(f"  <chi,chi> = {chi_norm:.6f};  constituents: {len(n1)} of dim 1, "
              f"{len(n2)} of dim 2;  multiplicity-free: {result['multiplicity_free']}")
        print(f"  two-dimensional indices j = {two_dim_j}; faithful (a^q acts by -1, "
              f"admissible carriers) j = {faithful_j}: {len(faithful_j)} of {len(two_dim_j)}")
        print(f"  c_Sigma = {c_sigma}: constituent {result['c_sigma_is_a_constituent']}, "
              f"admissible carrier {result['c_sigma_is_an_admissible_carrier']}")
        print(f"  sign pair {{c_Sigma, q-c_Sigma}} = {result['sign_pair']}: both constituents "
              f"{result['sign_pair_are_constituents']}; sign-completed rule selects j = "
              f"{result['sign_completed_selection']}, admissible: "
              f"{result['sign_completed_rule_selects_an_admissible_carrier']}")
        print(f"  multiplicities: worst |Im| {worst_imag:.2e}, worst integer defect "
              f"{worst_int:.2e}")
        if result["published_shape_check"] is not None:
            print(f"  published-shape kill-switch: "
                  f"{all(result['published_shape_check'].values())} "
                  f"{result['published_shape_check']}")
        print(f"  dimension check: {' + '.join(str(p['dim']) for p in proj_report)} "
              f"= {int(round(dim_from_mult))} = dim W ({result['dimension_check']['ok']})")
        print(f"  isotypic projectors: completeness defect {completeness:.2e}; "
              f"all ranks as expected: "
              f"{all(p['rank'] == p['expected_rank'] for p in proj_report)}")
        print(f"  two-dimensional constituents: "
              f"{[p['irrep'] for p in n2]}")
        print(f"  ALL CHECKS PASSED: {passed}")
    return result


def selector_sweep(per_prime, primes=(53, 101), seed=99, max_draws=20000):
    """Audit the sign-completed selector rule on dicyclic level-1 blocks.

    Oriented rule:        j = c_Sigma.
    Sign-completed rule:  j in {c_Sigma, q - c_Sigma}, then take the unique odd member.
    Exactly one member is odd because q is odd.  Reports both, so the refutation of the
    oriented rule and the survival of the sign-completed one are both traceable.
    """
    from o23_filtration_stabiliser import bfs_shells as _bfs
    rows, ok, ok_oriented, failed_checks = [], 0, 0, 0
    per_prime_found, exhausted = {}, {}
    distinct = set()
    for q in primes:
        rng = np.random.default_rng(seed)
        gens = build_generators(q)
        gens_arr = np.array(gens, dtype=np.int64)
        shells = _bfs(gens, q, min(int(0.30 * q ** 3), 400_000))
        seen = set()
        for _ in range(max_draws):
            c = tuple(int(v) for v in rng.integers(1, q, size=3))
            if sum(c) % q == 0 or c in seen:
                continue
            F = exact_frequency_set(np.array(shells[0], dtype=np.int64),
                                    np.array(c), gens_arr, q)
            if len(F) >= q:
                continue
            M = multiplier_group(F, q)
            if len(M) != 4 or not any((s * s) % q == q - 1 for s in M):
                continue
            seen.add(c)
            distinct.add((q, frozenset(int(f) for f in F)))
            r = decompose(q, c, 1, verbose=False)
            rows.append({k: r[k] for k in
                         ("q", "block", "c_sigma", "two_dim_j", "faithful_j", "sign_pair",
                          "sign_pair_are_constituents", "sign_completed_selection",
                          "sign_completed_rule_selects_an_admissible_carrier",
                          "c_sigma_is_an_admissible_carrier", "multiplicity_free",
                          "all_checks_passed")})
            # A row only counts as a success if the decomposition it rests on passed ITS OWN
            # checks.  Counting the rule alone made the guard vacuous: with a threshold forced
            # to zero every row failed while the sweep still reported full success.
            if not r["all_checks_passed"]:
                failed_checks += 1
            ok += bool(r["sign_completed_rule_selects_an_admissible_carrier"]
                       and r["sign_pair_are_constituents"]
                       and r["all_checks_passed"])
            ok_oriented += bool(r["c_sigma_is_an_admissible_carrier"]
                                and r["all_checks_passed"])
            if len(seen) >= per_prime:
                break
        per_prime_found[q] = len(seen)
        exhausted[q] = len(seen) < per_prime
    return {
        "manifest": {"per_prime_requested": per_prime, "primes": list(primes), "seed": seed,
                     "max_draws": max_draws, "level": 1,
                     "thresholds": {"invariance": TOL_INVARIANCE, "homomorphism": TOL_HOM,
                                    "multiplicity": TOL_MULT, "projector": TOL_PROJ}},
        "primes": list(primes), "levels": len(rows),
        "levels_per_prime": per_prime_found,
        "search_exhausted_before_quota": exhausted,
        "distinct_exact_subspaces": len(distinct),
        "rows_failing_own_checks": failed_checks,
        "successes": ok, "oriented_successes": ok_oriented, "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int)
    ap.add_argument("--block", type=str, help="comma-separated c1,c2,c3")
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--selector-sweep", type=int, default=0, dest="sweep",
                    help="audit the sign-completed selector rule on this many further "
                         "dicyclic level-1 blocks per prime, and deposit the evidence")
    ap.add_argument("--sweep-seed", type=int, default=99, dest="sweep_seed",
                    help="seed for the selector sweep (recorded in its manifest); the other "
                         "script's runs use their own --seed, default 12345")
    args = ap.parse_args()

    if args.q and args.block:
        cases = [(args.q, tuple(int(v) for v in args.block.split(",")), args.level)]
        # An ad-hoc case must never overwrite the deposited combined evidence file.
        out = CKPT_DIR / (f"dicyclic_q{args.q}_c{args.block.replace(',', '-')}"
                          f"_n{args.level}.json")
    else:
        cases = WITNESSES
        out = OUT

    if args.sweep:
        out = CKPT_DIR / f"selector_rule_sweep_{args.sweep}_seed{args.sweep_seed}.json"
        sweep = selector_sweep(args.sweep, seed=args.sweep_seed)
        CKPT_DIR.mkdir(exist_ok=True)
        out.write_text(json.dumps(sweep, indent=2, default=str))
        print(f"\nsign-completed rule: {sweep['successes']}/{sweep['levels']} rows "
              f"({sweep['primes']}); oriented rule j = c_Sigma: "
              f"{sweep['oriented_successes']}/{sweep['levels']}")
        print(f"rows per prime: {sweep['levels_per_prime']}; distinct exact W_<1 subspaces: "
              f"{sweep['distinct_exact_subspaces']}; rows failing their own checks: "
              f"{sweep['rows_failing_own_checks']}")
        for q_, ex in sweep["search_exhausted_before_quota"].items():
            if ex:
                print(f"  WARNING: search exhausted at q={q_} before reaching the quota "
                      f"({sweep['levels_per_prime'][q_]} of {args.sweep} found)")
        print(f"evidence written to {out.name}")
        raise SystemExit(0 if (sweep["successes"] == sweep["levels"]
                               and sweep["rows_failing_own_checks"] == 0) else 6)

    results = [decompose(q, c, n) for q, c, n in cases]
    CKPT_DIR.mkdir(exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nevidence written to {out.name}")
    if not all(r["all_checks_passed"] for r in results):
        raise SystemExit(5)


if __name__ == "__main__":
    main()
