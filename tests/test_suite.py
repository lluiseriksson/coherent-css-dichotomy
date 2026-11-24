"""Comprehensive pytest suite for Paper 07, including exact CSS distance checks,
complex 2x2 SVD physics cross-checks, and 37 live adversarial forgery,
canonical-parser, schema, and strict type-injection attacks."""

from __future__ import annotations

import copy
import math
from decimal import Decimal
from fractions import Fraction
from typing import Any, Dict

import numpy as np
import pytest

from src import (
    build_code_catalog,
    build_coherent_css_certificate,
    compute_coset_weight_enumerators,
    compute_exact_fidelities,
    compute_swap_point_status,
    parse_canonical_fraction,
    validate_css_code,
    verify_coherent_css_certificate_exact,
)
from src.css_coherent import compute_exact_css_distances


def test_catalog_css_invariants_distances_and_svd_crosscheck() -> None:
    """All 9 catalog codes satisfy CSS axioms, exact (d_X, d_Z), and match complex 2x2 SVD."""
    catalog = build_code_catalog()
    assert len(catalog) == 9
    expected_distances = {
        "steane_7_1_3": (3, 3),
        "shor_9_1_3": (3, 3),
        "qrm_15_1_3": (3, 7),
        "rotsurf_d2_Xmaj": (2, 2),
        "rotsurf_d2_Zmaj": (2, 2),
        "rotsurf_d3": (3, 3),
        "rotsurf_d4_Xmaj": (4, 4),
        "rotsurf_d4_Zmaj": (4, 4),
        "swap_immune_8_1_2_odd": (2, 3),
    }
    for name, code in catalog.items():
        meta = validate_css_code(code)
        assert meta["k"] == 1
        assert meta["all_z_stabilizers_even"] is True
        assert compute_exact_css_distances(code) == expected_distances[name]
        enums = compute_coset_weight_enumerators(code)
        for p in (Fraction(1, 20), Fraction(1, 5), Fraction(3, 10), Fraction(1, 2)):
            res = compute_exact_fidelities(code, p, enumerators=enums)
            theta = math.asin(math.sqrt(float(p)))
            c_th, s_th = math.cos(theta), math.sin(theta)
            f_cptp_num = 0.0
            f_pauli_num = 0.0
            for _, (a_cnt, b_cnt) in enums.items():
                alpha = sum(c * (c_th ** (code.n - w)) * ((1j * s_th) ** w) for w, c in a_cnt.items())
                beta = sum(c * (c_th ** (code.n - w)) * ((1j * s_th) ** w) for w, c in b_cnt.items())
                K = np.array([[alpha + beta, 0.0], [0.0, alpha - beta]], dtype=complex)
                sv = np.linalg.svd(K, compute_uv=False)
                f_cptp_num += 0.25 * float(sv.sum() ** 2)
                f_pauli_num += max(abs(alpha) ** 2, abs(beta) ** 2)
            assert abs(f_cptp_num - float(res["F_CPTP"])) < 1e-12, (name, p)  # type: ignore[arg-type]
            assert abs(f_pauli_num - float(res["F_Pauli"])) < 1e-12, (name, p)  # type: ignore[arg-type]


def test_exact_benchmark_phenomena() -> None:
    """Verify Steane p_0(pi/4)=1, RotSurf d=2 cancellation, Shor p=1/4 unit fidelity, QRM penalty, and swap immunity."""
    cat = build_code_catalog()

    # 1. Steane [[7,1,3]]: p_0(pi/4) = 1, F_Pauli(pi/4) = 1/2, F_CPTP = 1
    r_steane = compute_exact_fidelities(cat["steane_7_1_3"], Fraction(1, 2))
    assert r_steane["p_trivial_syndrome"] == Fraction(1, 1)
    assert r_steane["F_Pauli"] == Fraction(1, 2)
    assert r_steane["F_CPTP"] == Fraction(1, 1)

    # 2. RotSurf d=2 X-majority: exact all-angle cancellation F_Pauli == F_CPTP == F_Twirl
    d2 = cat["rotsurf_d2_Xmaj"]
    d2_enums = compute_coset_weight_enumerators(d2)
    for k in range(0, 11):
        r = compute_exact_fidelities(d2, Fraction(k, 20), enumerators=d2_enums)
        assert r["F_Pauli"] == r["F_CPTP"] == r["F_Twirl"]

    # 2b. Even-d 90-degree rotated surface code orientation phase transition at pi/4
    for d_str in ("d2", "d4"):
        r_x = compute_exact_fidelities(cat[f"rotsurf_{d_str}_Xmaj"], Fraction(1, 2))
        r_z = compute_exact_fidelities(cat[f"rotsurf_{d_str}_Zmaj"], Fraction(1, 2))
        assert r_x["F_Pauli"] == Fraction(1, 2) and r_x["F_CPTP"] == Fraction(1, 2)
        assert r_z["F_Pauli"] == Fraction(1, 1) and r_z["F_CPTP"] == Fraction(1, 1)

    # 3. Shor [[9,1,3]]: >100x suppression at p=1/5 AND exact unit Pauli fidelity at p=1/4 (theta=pi/6)!
    r_shor_20 = compute_exact_fidelities(cat["shor_9_1_3"], Fraction(1, 5))
    assert r_shor_20["F_CPTP"] == Fraction(1, 1)
    assert (Fraction(1, 1) - r_shor_20["F_Twirl"]) > 100 * (Fraction(1, 1) - r_shor_20["F_Pauli"])  # type: ignore[operator]
    r_shor_25 = compute_exact_fidelities(cat["shor_9_1_3"], Fraction(1, 4))
    assert r_shor_25["F_Pauli"] == Fraction(1, 1)

    # 4. QRM [[15,1,3]]: coherent penalty F_Pauli < F_Twirl at p = 1/5 despite doubly-even S_Z
    r_qrm = compute_exact_fidelities(cat["qrm_15_1_3"], Fraction(1, 5))
    assert r_qrm["F_Pauli"] < r_qrm["F_Twirl"]  # type: ignore[operator]
    assert r_qrm["F_CPTP"] == Fraction(1, 1)

    # 5. Swap-immune witnesses: exact F_Pauli(pi/4) = 1 in BOTH even-|L_Z| and odd-|L_Z| branches
    for name in ("rotsurf_d2_Zmaj", "rotsurf_d4_Zmaj", "swap_immune_8_1_2_odd"):
        sw = compute_swap_point_status(cat[name])
        assert sw["regime"] == "swap_immune_unit"
        r50 = compute_exact_fidelities(cat[name], Fraction(1, 2))
        assert r50["F_Pauli"] == Fraction(1, 1)
        assert r50["F_CPTP"] == Fraction(1, 1)
        assert r50["zero_coset_syndromes"] == r50["num_syndromes"]


def test_all_9_certificates_roundtrip_verified() -> None:
    """Every catalog code builds a valid certificate accepted by verify_coherent_css_certificate_exact."""
    for name, code in build_code_catalog().items():
        cert = build_coherent_css_certificate(code, p_points=(Fraction(1, 5), Fraction(1, 2)))
        out = verify_coherent_css_certificate_exact(cert)
        assert out["verified"] is True
        assert out["code_name"] == name


def _base_certs() -> Dict[str, Dict[str, Any]]:
    cat = build_code_catalog()
    return {
        "steane": build_coherent_css_certificate(cat["steane_7_1_3"]),
        "d4": build_coherent_css_certificate(cat["rotsurf_d4_Xmaj"]),
        "immune": build_coherent_css_certificate(cat["rotsurf_d4_Zmaj"]),
    }


def test_37_adversarial_forgery_and_type_attacks_rejected() -> None:
    """37 distinct adversarial attacks must ALL be rejected with ValueError."""
    base = _base_certs()

    attacks: list[tuple[str, Dict[str, Any]]] = []

    # 1. Wrong schema_version
    c = copy.deepcopy(base["steane"]); c["schema_version"] = "COHERENT-CSS-CERT-0.9"
    attacks.append(("bad_schema_version", c))

    # 2. Missing required key
    c = copy.deepcopy(base["steane"]); del c["swap_point"]
    attacks.append(("missing_required_key", c))

    # 3. Float injection in n
    c = copy.deepcopy(base["steane"]); c["n"] = 7.0
    attacks.append(("float_in_n", c))

    # 4. Float injection inside nested evaluation field
    c = copy.deepcopy(base["steane"]); c["evaluations"][0]["F_Pauli"] = 0.5
    attacks.append(("float_in_eval", c))

    # 5. Bool injection (True instead of 1) in H_X
    c = copy.deepcopy(base["steane"]); c["H_X"][0][0] = True
    attacks.append(("bool_in_H_X", c))

    # 6. Bool injection in k
    c = copy.deepcopy(base["steane"]); c["k"] = True
    attacks.append(("bool_in_k", c))

    # 7. Decimal injection
    c = copy.deepcopy(base["steane"]); c["evaluations"][0]["p"] = Decimal("0.5")
    attacks.append(("decimal_injection", c))

    # 8. Tuple injection
    c = copy.deepcopy(base["steane"]); c["L_Z"] = tuple(c["L_Z"])
    attacks.append(("tuple_injection", c))

    # 9. Non-string dict key injection
    c = copy.deepcopy(base["steane"]); c[123] = "rogue"  # type: ignore[index]
    attacks.append(("non_str_dict_key", c))

    # 10. Unreduced rational string "2/4"
    c = copy.deepcopy(base["steane"]); c["evaluations"][1]["p"] = "2/4"
    attacks.append(("unreduced_fraction", c))

    # 11. Trailing "/1" in canonical rational ("1/1")
    c = copy.deepcopy(base["steane"]); c["evaluations"][0]["F_CPTP"] = "1/1"
    attacks.append(("trailing_slash_one", c))

    # 12. Leading '+' in rational string
    c = copy.deepcopy(base["steane"]); c["evaluations"][1]["p"] = "+1/2"
    attacks.append(("leading_plus_rational", c))

    # 13. Leading whitespace in rational string
    c = copy.deepcopy(base["steane"]); c["evaluations"][1]["p"] = " 1/2"
    attacks.append(("whitespace_rational", c))

    # 14. Trailing newline injection in fraction ("1/2\n")
    c = copy.deepcopy(base["steane"]); c["evaluations"][1]["p"] = "1/2\n"
    attacks.append(("trailing_newline_fraction", c))

    # 15. Trailing newline injection in integer ("1\n")
    c = copy.deepcopy(base["steane"]); c["evaluations"][0]["F_CPTP"] = "1\n"
    attacks.append(("trailing_newline_int", c))

    # 16. Non-commuting H_X / H_Z bit flip
    c = copy.deepcopy(base["steane"]); c["H_X"][0][6] = 1
    attacks.append(("non_commuting_HX_HZ", c))

    # 17. Odd-weight Z-stabilizer row (flip 1 bit of H_Z)
    c = copy.deepcopy(base["steane"]); c["H_Z"][0][0] = 0
    attacks.append(("odd_weight_HZ", c))

    # 18. L_Z inside rowspace(H_Z)
    c = copy.deepcopy(base["steane"]); c["L_Z"] = list(c["H_Z"][0])
    attacks.append(("LZ_in_HZ_span", c))

    # 19. L_X inside rowspace(H_X)
    c = copy.deepcopy(base["steane"]); c["L_X"] = list(c["H_X"][0])
    attacks.append(("LX_in_HX_span", c))

    # 20. Forged parity_branch on RotSurf d=4 X-maj
    c = copy.deepcopy(base["d4"]); c["parity_branch"] = "odd_lz_unitary_cptp"
    attacks.append(("forged_parity_branch", c))

    # 21. Forged F_CPTP = "1" on even-|L_Z| RotSurf d=4 X-maj
    c = copy.deepcopy(base["d4"]); c["evaluations"][0]["F_CPTP"] = "1"
    attacks.append(("forged_F_CPTP_even_lz", c))

    # 22. Perturbed F_Pauli by 1/1000000
    c = copy.deepcopy(base["steane"])
    orig = parse_canonical_fraction(c["evaluations"][0]["F_Pauli"])
    pert = orig + Fraction(1, 1000000)
    c["evaluations"][0]["F_Pauli"] = f"{pert.numerator}/{pert.denominator}"
    attacks.append(("perturbed_F_Pauli", c))

    # 23. Forged swap_point regime on Steane
    c = copy.deepcopy(base["steane"])
    c["swap_point"]["regime"] = "swap_immune_unit"
    c["swap_point"]["predicted_F_Pauli_pi4"] = "1"
    attacks.append(("forged_swap_regime", c))

    # 24. Missing immunity witness on rotsurf_d4_Zmaj
    c = copy.deepcopy(base["immune"])
    c["swap_point"]["immunity_witness_g_star"] = None
    attacks.append(("missing_immunity_witness", c))

    # 25. Forged witness not in S_Z on rotsurf_d4_Zmaj
    c = copy.deepcopy(base["immune"])
    c["swap_point"]["immunity_witness_g_star"] = [1] + [0] * 15
    attacks.append(("witness_not_in_SZ", c))

    # 26. Spurious non-null witness on collapse-regime Steane
    c = copy.deepcopy(base["steane"])
    c["swap_point"]["immunity_witness_g_star"] = list(c["H_Z"][0])
    attacks.append(("spurious_witness_on_collapse", c))

    # 27. Missing mandatory p = "1/2" evaluation
    c = copy.deepcopy(base["steane"])
    c["evaluations"] = [ev for ev in c["evaluations"] if ev["p"] != "1/2"]
    attacks.append(("missing_pi4_eval", c))

    # 28. Duplicate weight entry in trivial_syndrome_enumerator
    c = copy.deepcopy(base["steane"])
    c["trivial_syndrome_enumerator"]["coset_0"].append([0, 1])
    attacks.append(("duplicate_weight_enumerator", c))

    # 29. Reversed (non-ascending) weight ordering in trivial_syndrome_enumerator
    c = copy.deepcopy(base["steane"])
    c["trivial_syndrome_enumerator"]["coset_0"] = list(reversed(c["trivial_syndrome_enumerator"]["coset_0"]))
    attacks.append(("reversed_weight_enumerator", c))

    # 30. Non-int string in trivial_syndrome_enumerator (must raise ValueError, NOT TypeError)
    c = copy.deepcopy(base["steane"])
    c["trivial_syndrome_enumerator"]["coset_0"][0] = ["0", 1]
    attacks.append(("string_weight_enumerator", c))

    # 31. Forged swap_point.K_Z_basis string
    c = copy.deepcopy(base["steane"])
    c["swap_point"]["K_Z_basis"] = "TOTALLY_FAKE_BASIS"
    attacks.append(("fake_KZ_basis_string", c))

    # 32. Forged swap_point.K_Z_basis linearly dependent / wrong subspace rows
    c = copy.deepcopy(base["steane"])
    c["swap_point"]["K_Z_basis"][0] = list(c["swap_point"]["K_Z_basis"][1])
    attacks.append(("dependent_KZ_basis", c))

    # 33. Missing swap_point.K_Z_basis key
    c = copy.deepcopy(base["steane"])
    del c["swap_point"]["K_Z_basis"]
    attacks.append(("missing_KZ_basis", c))

    # 34. Extra rogue key at top level
    c = copy.deepcopy(base["steane"]); c["rogue_top_key"] = "forged_theorem"
    attacks.append(("rogue_top_key", c))

    # 35. Extra rogue key in swap_point and evaluations[0]
    c = copy.deepcopy(base["steane"]); c["swap_point"]["rogue_swap"] = 1
    attacks.append(("rogue_swap_key", c))

    # 36. Extra rogue key in evaluations[0]
    c = copy.deepcopy(base["steane"]); c["evaluations"][0]["rogue_eval"] = "0"
    attacks.append(("rogue_eval_key", c))

    # 37. Non-string code_name / family injection
    c = copy.deepcopy(base["steane"]); c["code_name"] = None; c["family"] = {"evil": 1}
    attacks.append(("non_str_code_name_family", c))

    assert len(attacks) == 37
    for attack_name, payload in attacks:
        with pytest.raises(ValueError, match=r".+"):
            verify_coherent_css_certificate_exact(payload)
