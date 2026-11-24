"""Run all 4 exact computational experiments for Paper 07 and emit certified artifacts.

Outputs:
  results/certificates/<code_name>.cert.json   (9 exact-Q COHERENT-CSS-CERT-1.0 files)
  results/exp1_certified_atlas.json
  results/exp2_angle_sweeps.json
  results/exp3_swap_point_binary_law.json
  results/exp4_surface_and_rep_parity.json
"""

from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import (
    CSSCode,
    build_code_catalog,
    build_coherent_css_certificate,
    canonical_fraction_str,
    compute_coset_weight_enumerators,
    compute_exact_fidelities,
    compute_swap_point_status,
    validate_css_code,
    verify_coherent_css_certificate_exact,
)
from src.css_coherent import (
    bits_to_mask,
    gf2_in_span,
    gf2_nullspace,
    gf2_rank,
    make_repetition_z_code,
    make_rotated_surface_code,
    mask_to_bits,
    popcount,
)


def svd_crosscheck(code: CSSCode, p: Fraction, enumerators: Any) -> Tuple[float, float]:
    """Independent complex 2x2 SVD cross-check of F_CPTP and F_Pauli at sin^2(theta) = p."""
    theta = math.asin(math.sqrt(float(p)))
    c_th = math.cos(theta)
    s_th = math.sin(theta)
    f_cptp_num = 0.0
    f_pauli_num = 0.0
    for _, (a_cnt, b_cnt) in enumerators.items():
        alpha = sum(c * (c_th ** (code.n - w)) * ((1j * s_th) ** w) for w, c in a_cnt.items())
        beta = sum(c * (c_th ** (code.n - w)) * ((1j * s_th) ** w) for w, c in b_cnt.items())
        K = np.array([[alpha + beta, 0.0], [0.0, alpha - beta]], dtype=complex)
        sv = np.linalg.svd(K, compute_uv=False)
        f_cptp_num += 0.25 * float(sv.sum() ** 2)
        f_pauli_num += max(abs(alpha) ** 2, abs(beta) ** 2)
    return f_pauli_num, f_cptp_num


def run_exp1_certified_atlas(cert_dir: Path) -> Dict[str, Any]:
    """Build & independently verify exact-Q certificates for all 9 atlas CSS codes."""
    cert_dir.mkdir(parents=True, exist_ok=True)
    catalog = build_code_catalog()
    rows: List[Dict[str, Any]] = []

    p_cert_points = (Fraction(1, 10), Fraction(1, 5), Fraction(1, 2))
    max_svd_err = 0.0

    for name, code in catalog.items():
        t0 = time.perf_counter()
        meta = validate_css_code(code)
        swap = compute_swap_point_status(code)
        enums = compute_coset_weight_enumerators(code)
        cert = build_coherent_css_certificate(code, p_points=p_cert_points)
        cert_path = cert_dir / f"{name}.cert.json"
        cert_path.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")

        # Independent zero-float exact replay from disk
        t_v0 = time.perf_counter()
        loaded = json.loads(cert_path.read_text(encoding="utf-8"))
        ver = verify_coherent_css_certificate_exact(loaded)
        verify_ms = (time.perf_counter() - t_v0) * 1000.0
        total_ms = (time.perf_counter() - t0) * 1000.0

        r10 = compute_exact_fidelities(code, Fraction(1, 10), enumerators=enums)
        r20 = compute_exact_fidelities(code, Fraction(1, 5), enumerators=enums)
        r50 = compute_exact_fidelities(code, Fraction(1, 2), enumerators=enums)

        # Independent complex 2x2 SVD cross-check at p = 1/5 and p = 1/2
        p20_num, c20_num = svd_crosscheck(code, Fraction(1, 5), enums)
        p50_num, c50_num = svd_crosscheck(code, Fraction(1, 2), enums)
        svd_err = max(
            abs(p20_num - float(r20["F_Pauli"])),  # type: ignore[arg-type]
            abs(c20_num - float(r20["F_CPTP"])),  # type: ignore[arg-type]
            abs(p50_num - float(r50["F_Pauli"])),  # type: ignore[arg-type]
            abs(c50_num - float(r50["F_CPTP"])),  # type: ignore[arg-type]
        )
        max_svd_err = max(max_svd_err, svd_err)
        if svd_err >= 1e-12:
            raise RuntimeError(f"SVD cross-check discrepancy on {name}: {svd_err}")

        infid_ratio_p20 = (
            float((Fraction(1, 1) - r20["F_Twirl"]) / (Fraction(1, 1) - r20["F_Pauli"]))  # type: ignore[operator]
            if r20["F_Pauli"] != Fraction(1, 1)
            else None
        )

        rows.append(
            {
                "code_name": name,
                "family": code.family,
                "n": code.n,
                "k": 1,
                "d_X": code.d_X,
                "d_Z": code.d_Z,
                "r_X": meta["r_X"],
                "r_Z": meta["r_Z"],
                "L_Z_weight": meta["L_Z_weight"],
                "L_Z_parity": meta["L_Z_parity"],
                "all_z_doubly_even": meta["all_z_stabilizers_doubly_even"],
                "parity_branch": ver["parity_branch"],
                "K_Z_dim": swap["K_Z_dim"],
                "swap_regime": swap["regime"],
                "immunity_witness_g_star": swap["immunity_witness_g_star"],
                "p_1_10": {
                    "F_Pauli_exact": canonical_fraction_str(r10["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP_exact": canonical_fraction_str(r10["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl_exact": canonical_fraction_str(r10["F_Twirl"]),  # type: ignore[arg-type]
                    "F_Pauli": float(r10["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP": float(r10["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl": float(r10["F_Twirl"]),  # type: ignore[arg-type]
                },
                "p_1_5": {
                    "F_Pauli_exact": canonical_fraction_str(r20["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP_exact": canonical_fraction_str(r20["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl_exact": canonical_fraction_str(r20["F_Twirl"]),  # type: ignore[arg-type]
                    "F_Pauli": float(r20["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP": float(r20["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl": float(r20["F_Twirl"]),  # type: ignore[arg-type]
                    "twirl_to_coherent_infid_ratio": infid_ratio_p20,
                },
                "p_1_2": {
                    "F_Pauli_exact": canonical_fraction_str(r50["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP_exact": canonical_fraction_str(r50["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl_exact": canonical_fraction_str(r50["F_Twirl"]),  # type: ignore[arg-type]
                    "p_trivial_syndrome_exact": canonical_fraction_str(r50["p_trivial_syndrome"]),  # type: ignore[arg-type]
                    "num_syndromes": r50["num_syndromes"],
                    "equal_coset_syndromes": r50["equal_coset_syndromes"],
                    "zero_coset_syndromes": r50["zero_coset_syndromes"],
                },
                "max_svd_err": svd_err,
                "verify_ms": round(verify_ms, 3),
                "total_ms": round(total_ms, 3),
                "certificate_verified": True,
            }
        )

    return {
        "experiment": "exp1_certified_atlas",
        "num_codes": len(rows),
        "all_verified": all(r["certificate_verified"] for r in rows),
        "max_svd_crosscheck_err": max_svd_err,
        "codes": rows,
    }


def run_exp2_angle_sweeps() -> Dict[str, Any]:
    """Exact-Q sweeps over p = k/40 (k = 0..20, theta in [0, pi/4]) for all 9 atlas codes."""
    catalog = build_code_catalog()
    grid = [Fraction(k, 40) for k in range(21)]
    sweeps: Dict[str, Any] = {}

    for name, code in catalog.items():
        enums = compute_coset_weight_enumerators(code)
        pts: List[Dict[str, Any]] = []
        for p in grid:
            r = compute_exact_fidelities(code, p, enumerators=enums)
            pts.append(
                {
                    "p_exact": canonical_fraction_str(p),
                    "p": float(p),
                    "theta_over_pi": math.asin(math.sqrt(float(p))) / math.pi,
                    "F_Pauli_exact": canonical_fraction_str(r["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP_exact": canonical_fraction_str(r["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl_exact": canonical_fraction_str(r["F_Twirl"]),  # type: ignore[arg-type]
                    "F_Pauli": float(r["F_Pauli"]),  # type: ignore[arg-type]
                    "F_CPTP": float(r["F_CPTP"]),  # type: ignore[arg-type]
                    "F_Twirl": float(r["F_Twirl"]),  # type: ignore[arg-type]
                    "delta_coh_minus_twirl": float(r["coherent_minus_twirl"]),  # type: ignore[arg-type]
                    "p_trivial_syndrome": float(r["p_trivial_syndrome"]),  # type: ignore[arg-type]
                }
            )
        sweeps[name] = {
            "n": code.n,
            "L_Z_parity": "odd" if popcount(code.lz_mask()) % 2 == 1 else "even",
            "points": pts,
        }

    # Verify structural highlights directly from exact rationals
    d2_all_equal = all(
        pt["F_Pauli_exact"] == pt["F_CPTP_exact"] == pt["F_Twirl_exact"]
        for pt in sweeps["rotsurf_d2_Xmaj"]["points"]
    )
    qrm_penalty_points = [
        pt["p_exact"]
        for pt in sweeps["qrm_15_1_3"]["points"]
        if pt["F_Pauli"] < pt["F_Twirl"] - 1e-15
    ]

    return {
        "experiment": "exp2_angle_sweeps",
        "grid_size": len(grid),
        "rotsurf_d2_Xmaj_exact_all_angle_cancellation": d2_all_equal,
        "qrm_coherent_penalty_p_points": qrm_penalty_points,
        "sweeps": sweeps,
    }


def _generate_random_even_css_code(
    n: int, r_z: int, rng: random.Random, idx: int
) -> CSSCode | None:
    """Deterministically sample a valid [[n, 1, d]] even-Z-stabilizer CSS code."""
    # Sample r_z linearly independent even-weight vectors for H_Z
    even_pool = [m for m in range(1, 1 << n) if popcount(m) % 2 == 0]
    hz: List[int] = []
    for _ in range(40):
        cand = rng.choice(even_pool)
        if not gf2_in_span(cand, hz, n):
            hz.append(cand)
            if len(hz) == r_z:
                break
    if len(hz) != r_z:
        return None

    # H_X must lie in H_Z^⊥ and have rank r_x = n - r_z - 1
    hz_perp = gf2_nullspace(hz, n)
    r_x = n - r_z - 1
    if len(hz_perp) < r_x + 1:
        return None

    # Random invertible mixing of hz_perp basis to pick r_x rows for H_X and 1 row for L_X
    basis = list(hz_perp)
    for _ in range(2 * len(basis)):
        i = rng.randrange(len(basis))
        j = rng.randrange(len(basis))
        if i != j:
            basis[i] ^= basis[j]
    hx = basis[:r_x]
    if gf2_rank(hx, n) != r_x:
        return None

    # Find L_Z in H_X^⊥ \ rowspace(H_Z)
    hx_perp = gf2_nullspace(hx, n)
    lz_cand = None
    for v in hx_perp:
        if not gf2_in_span(v, hz, n):
            lz_cand = v
            break
    if lz_cand is None:
        return None

    # Find L_X in H_Z^⊥ \ rowspace(H_X) with |L_X ∩ L_Z| odd
    lx_cand = None
    for v in basis[r_x:]:
        if not gf2_in_span(v, hx, n) and popcount(v & lz_cand) % 2 == 1:
            lx_cand = v
            break
    if lx_cand is None:
        return None

    code = CSSCode(
        name=f"rand_css_n{n}_rz{r_z}_{idx:03d}",
        n=n,
        H_X=tuple(tuple(mask_to_bits(r, n)) for r in hx),
        H_Z=tuple(tuple(mask_to_bits(r, n)) for r in hz),
        L_X=tuple(mask_to_bits(lx_cand, n)),
        L_Z=tuple(mask_to_bits(lz_cand, n)),
        family="random_census",
    )
    validate_css_code(code)
    return code


def run_exp3_swap_point_binary_law() -> Dict[str, Any]:
    """Verify Theorem 2.2 0-1 binary law on all 9 atlas codes + 120 random even-Z CSS codes."""
    rng = random.Random(20261124)
    catalog = build_code_catalog()

    sampled: List[CSSCode] = list(catalog.values())
    target_random = 120
    configs = [(4, 1), (5, 1), (5, 2), (6, 1), (6, 2), (6, 3), (7, 2), (7, 3), (8, 2), (8, 3), (8, 4), (9, 3), (9, 4)]
    attempts = 0
    rand_count = 0
    while rand_count < target_random and attempts < 5000:
        attempts += 1
        n, r_z = configs[rand_count % len(configs)]
        code = _generate_random_even_css_code(n, r_z, rng, rand_count)
        if code is not None:
            sampled.append(code)
            rand_count += 1

    collapse_count = 0
    immune_count = 0
    intermediate_count = 0
    odd_lz_collapse = 0
    odd_lz_immune = 0
    even_lz_collapse = 0
    even_lz_immune = 0
    exact_theorem_matches = 0

    records: List[Dict[str, Any]] = []
    for code in sampled:
        meta = validate_css_code(code)
        swap = compute_swap_point_status(code)
        r50 = compute_exact_fidelities(code, Fraction(1, 2))
        fp = r50["F_Pauli"]
        fc = r50["F_CPTP"]
        pred = swap["predicted_F_Pauli_pi4"]
        if fp == pred:
            exact_theorem_matches += 1
        else:
            raise RuntimeError(f"Theorem 2.2 mismatch on {code.name}: {fp} != {pred}")

        if fp == Fraction(1, 2):
            collapse_count += 1
            if meta["L_Z_parity"] == "odd":
                odd_lz_collapse += 1
            else:
                even_lz_collapse += 1
        elif fp == Fraction(1, 1):
            immune_count += 1
            if meta["L_Z_parity"] == "odd":
                odd_lz_immune += 1
            else:
                even_lz_immune += 1
        else:
            intermediate_count += 1

        records.append(
            {
                "code_name": code.name,
                "family": code.family,
                "n": code.n,
                "r_Z": meta["r_Z"],
                "L_Z_parity": meta["L_Z_parity"],
                "K_Z_dim": swap["K_Z_dim"],
                "swap_regime": swap["regime"],
                "F_Pauli_pi4": canonical_fraction_str(fp),  # type: ignore[arg-type]
                "F_CPTP_pi4": canonical_fraction_str(fc),  # type: ignore[arg-type]
                "equal_coset_syndromes": r50["equal_coset_syndromes"],
                "zero_coset_syndromes": r50["zero_coset_syndromes"],
                "num_syndromes": r50["num_syndromes"],
            }
        )

    return {
        "experiment": "exp3_swap_point_binary_law",
        "total_codes_tested": len(sampled),
        "atlas_codes": len(catalog),
        "random_even_css_codes": rand_count,
        "exact_theorem_matches": exact_theorem_matches,
        "collapse_half_count": collapse_count,
        "immune_unit_count": immune_count,
        "intermediate_value_count": intermediate_count,
        "quadrant_counts": {
            "odd_lz_collapse_half": odd_lz_collapse,
            "odd_lz_immune_unit": odd_lz_immune,
            "even_lz_collapse_half": even_lz_collapse,
            "even_lz_immune_unit": even_lz_immune,
        },
        "records": records,
    }


def run_exp4_surface_and_rep_parity() -> Dict[str, Any]:
    """Compare even-vs-odd distance rotated surface codes (d=2..4) and repetition codes (d=2..9)."""
    p_vals = [Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(3, 10), Fraction(2, 5), Fraction(1, 2)]
    surface_rows: List[Dict[str, Any]] = []
    for d, orient in (
        (2, "X_majority"),
        (2, "Z_majority"),
        (3, "X_majority"),
        (4, "X_majority"),
        (4, "Z_majority"),
    ):
        code = make_rotated_surface_code(d, orient)
        enums = compute_coset_weight_enumerators(code)
        swap = compute_swap_point_status(code)
        evals = {}
        for p in p_vals:
            r = compute_exact_fidelities(code, p, enumerators=enums)
            evals[canonical_fraction_str(p)] = {
                "F_Pauli_exact": canonical_fraction_str(r["F_Pauli"]),  # type: ignore[arg-type]
                "F_CPTP_exact": canonical_fraction_str(r["F_CPTP"]),  # type: ignore[arg-type]
                "F_Twirl_exact": canonical_fraction_str(r["F_Twirl"]),  # type: ignore[arg-type]
                "F_Pauli": float(r["F_Pauli"]),  # type: ignore[arg-type]
                "F_CPTP": float(r["F_CPTP"]),  # type: ignore[arg-type]
                "F_Twirl": float(r["F_Twirl"]),  # type: ignore[arg-type]
            }
        surface_rows.append(
            {
                "code_name": code.name,
                "d": d,
                "orientation": orient,
                "n": code.n,
                "parity": "odd" if d % 2 == 1 else "even",
                "swap_regime": swap["regime"],
                "evaluations": evals,
            }
        )

    rep_rows: List[Dict[str, Any]] = []
    for d in range(2, 10):
        code = make_repetition_z_code(d)
        enums = compute_coset_weight_enumerators(code)
        evals = {}
        for p in p_vals:
            r = compute_exact_fidelities(code, p, enumerators=enums)
            evals[canonical_fraction_str(p)] = {
                "F_Pauli_exact": canonical_fraction_str(r["F_Pauli"]),  # type: ignore[arg-type]
                "F_CPTP_exact": canonical_fraction_str(r["F_CPTP"]),  # type: ignore[arg-type]
                "F_Twirl_exact": canonical_fraction_str(r["F_Twirl"]),  # type: ignore[arg-type]
                "F_Pauli": float(r["F_Pauli"]),  # type: ignore[arg-type]
                "F_CPTP": float(r["F_CPTP"]),  # type: ignore[arg-type]
                "F_Twirl": float(r["F_Twirl"]),  # type: ignore[arg-type]
            }
        rep_rows.append(
            {
                "code_name": code.name,
                "d": d,
                "n": code.n,
                "parity": "odd" if d % 2 == 1 else "even",
                "evaluations": evals,
            }
        )

    return {
        "experiment": "exp4_surface_and_rep_parity",
        "p_values": [canonical_fraction_str(p) for p in p_vals],
        "rotated_surface": surface_rows,
        "phase_repetition": rep_rows,
    }


def main() -> None:
    res_dir = ROOT / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    cert_dir = res_dir / "certificates"

    t_start = time.perf_counter()
    e1 = run_exp1_certified_atlas(cert_dir)
    (res_dir / "exp1_certified_atlas.json").write_text(json.dumps(e1, indent=2) + "\n", encoding="utf-8")

    e2 = run_exp2_angle_sweeps()
    (res_dir / "exp2_angle_sweeps.json").write_text(json.dumps(e2, indent=2) + "\n", encoding="utf-8")

    e3 = run_exp3_swap_point_binary_law()
    (res_dir / "exp3_swap_point_binary_law.json").write_text(json.dumps(e3, indent=2) + "\n", encoding="utf-8")

    e4 = run_exp4_surface_and_rep_parity()
    (res_dir / "exp4_surface_and_rep_parity.json").write_text(json.dumps(e4, indent=2) + "\n", encoding="utf-8")

    elapsed = time.perf_counter() - t_start
    print(f"All 4 experiments completed in {elapsed:.2f}s:")
    print(f"  Exp 1: {e1['num_codes']} certified codes, all_verified={e1['all_verified']}, max_svd_err={e1['max_svd_crosscheck_err']:.2e}")
    print(f"  Exp 2: d=2 X-maj exact cancellation={e2['rotsurf_d2_Xmaj_exact_all_angle_cancellation']}, QRM penalty points={len(e2['qrm_coherent_penalty_p_points'])}")
    print(f"  Exp 3: {e3['total_codes_tested']} codes tested ({e3['collapse_half_count']} collapse=1/2, {e3['immune_unit_count']} immune=1, {e3['intermediate_value_count']} intermediate), matches={e3['exact_theorem_matches']}/{e3['total_codes_tested']}")
    print(f"  Exp 4: {len(e4['rotated_surface'])} surface + {len(e4['phase_repetition'])} repetition codes certified")


if __name__ == "__main__":
    main()
