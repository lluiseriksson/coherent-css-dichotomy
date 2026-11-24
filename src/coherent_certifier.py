"""Zero-floating-point-trust exact-Q certificate builder and verifier for coherent CSS decoding.

Schema: COHERENT-CSS-CERT-1.0
Security hardening (learned from hostile ARR-ASSESS-1.0 audits):
  - Strict recursive type whitelisting: ONLY plain JSON types (dict with str keys, list,
    str, strict int [excluding bool!], None) are allowed anywhere in the certificate tree.
    Floats, bools, Decimals, complex numbers, sets, tuples, or custom objects are rejected.
  - Strict canonical rational string enforcement: every rational field MUST match
    `"0"`, `"-?[1-9][0-9]*"`, or `"-?[1-9][0-9]*/[2-9][0-9]*"` / reduced coprime `"p/q"`
    with q >= 2, gcd(|p|, q) == 1, no whitespace, no `+` sign, no leading zeros, no `/1`.
  - Complete independent replay from the raw CSS matrices (H_X, H_Z, L_X, L_Z):
      1. Validates all F_2 CSS axioms (commutativity, full row ranks, k=1, logical pair,
         all Z-stabilizers even weight).
      2. Recomputes K_Z = S_Z ∩ S_Z^⊥ and verifies the Theorem 2.2 swap-point regime
         (`swap_collapse_half` vs `swap_immune_unit`) and any claimed witness g*.
      3. Recomputes the exact coset weight enumerators across all 2^{r_X} syndromes
         (or verifies the embedded digest & exact partition sum == 2^n).
      4. Independently replays F_e^{Pauli}(p), F_e^{CPTP}(p), F_e^{twirl}(p), and
         p_0(p) at every certified rational point p ∈ Q ∩ [0, 1] (mandating p = 1/2
         to verify the Theorem 2.2 binary law on every certificate!).
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from typing import Any, Dict, List, Sequence, Tuple

from .css_coherent import (
    CSSCode,
    bits_to_mask,
    compute_coset_weight_enumerators,
    compute_exact_fidelities,
    compute_self_orthogonal_subcode,
    compute_swap_point_status,
    gf2_in_span,
    gf2_rank,
    popcount,
    validate_css_code,
)

SCHEMA_VERSION = "COHERENT-CSS-CERT-1.0"
_INT_RE = re.compile(r"\A(?:0|-?[1-9][0-9]*)\Z")
_FRAC_RE = re.compile(r"\A(?P<num>-?[1-9][0-9]*)/(?P<den>[1-9][0-9]*)\Z")

_TOP_KEYS = frozenset({
    "schema_version",
    "code_name",
    "family",
    "n",
    "k",
    "r_X",
    "r_Z",
    "H_X",
    "H_Z",
    "L_X",
    "L_Z",
    "L_Z_weight",
    "parity_branch",
    "swap_point",
    "trivial_syndrome_enumerator",
    "total_error_count",
    "evaluations",
})
_SWAP_KEYS = frozenset({
    "regime",
    "predicted_F_Pauli_pi4",
    "K_Z_dim",
    "K_Z_basis",
    "immunity_witness_g_star",
})
_TRIV_KEYS = frozenset({"coset_0", "coset_1"})
_EVAL_KEYS = frozenset({
    "p",
    "F_Pauli",
    "F_CPTP",
    "F_Twirl",
    "p_trivial_syndrome",
    "num_syndromes",
    "equal_coset_syndromes",
    "zero_coset_syndromes",
})


def _assert_strict_json_types(obj: Any, path: str = "$") -> None:
    """Recursively reject floats, bools, Decimals, sets, tuples, or non-str dict keys."""
    if obj is None:
        return
    # Note: isinstance(True, int) is True in Python, so reject bool BEFORE int!
    if type(obj) is bool:
        raise ValueError(f"Forbidden bool at {path}: {obj!r}")
    if type(obj) is int:
        return
    if type(obj) is str:
        return
    if type(obj) is list:
        for idx, item in enumerate(obj):
            _assert_strict_json_types(item, f"{path}[{idx}]")
        return
    if type(obj) is dict:
        for k, v in obj.items():
            if type(k) is not str:
                raise ValueError(f"Non-string dict key at {path}: {k!r}")
            _assert_strict_json_types(v, f"{path}.{k}")
        return
    raise ValueError(f"Forbidden type {type(obj).__name__} at {path}: {obj!r}")


def canonical_fraction_str(x: Fraction) -> str:
    """Format an exact Fraction in strict reduced canonical form ('p' or 'p/q')."""
    if not isinstance(x, Fraction):
        raise TypeError(f"Expected Fraction, got {type(x).__name__}")
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def parse_canonical_fraction(s: Any, field_name: str = "rational") -> Fraction:
    """Parse a strictly canonical rational string ('0', 'p', or reduced 'p/q' with q >= 2)."""
    if type(s) is not str:
        raise ValueError(f"{field_name} must be a canonical string, got {type(s).__name__}")
    if _INT_RE.match(s):
        return Fraction(int(s), 1)
    m = _FRAC_RE.match(s)
    if not m:
        raise ValueError(f"Non-canonical rational string for {field_name}: {s!r}")
    num = int(m.group("num"))
    den = int(m.group("den"))
    if den < 2:
        raise ValueError(f"Denominator must be >= 2 for fraction {field_name}: {s!r}")
    if math.gcd(abs(num), den) != 1:
        raise ValueError(f"Unreduced rational string for {field_name}: {s!r}")
    return Fraction(num, den)


def build_coherent_css_certificate(
    code: CSSCode,
    p_points: Sequence[Fraction] = (Fraction(1, 5), Fraction(1, 2)),
) -> Dict[str, Any]:
    """Construct a self-contained, zero-float COHERENT-CSS-CERT-1.0 certificate."""
    meta = validate_css_code(code)
    if not meta["all_z_stabilizers_even"]:
        raise ValueError(f"Code {code.name} violates even-Z-stabilizer hypothesis")

    swap_info = compute_swap_point_status(code)
    enumerators = compute_coset_weight_enumerators(code)

    # Ensure p = 1/2 is always included to certify Theorem 2.2
    pts = list(p_points)
    if Fraction(1, 2) not in pts:
        pts.append(Fraction(1, 2))

    evaluations: List[Dict[str, Any]] = []
    for p in pts:
        res = compute_exact_fidelities(code, p, enumerators=enumerators)
        evaluations.append(
            {
                "p": canonical_fraction_str(p),
                "F_Pauli": canonical_fraction_str(res["F_Pauli"]),  # type: ignore[arg-type]
                "F_CPTP": canonical_fraction_str(res["F_CPTP"]),  # type: ignore[arg-type]
                "F_Twirl": canonical_fraction_str(res["F_Twirl"]),  # type: ignore[arg-type]
                "p_trivial_syndrome": canonical_fraction_str(res["p_trivial_syndrome"]),  # type: ignore[arg-type]
                "num_syndromes": int(res["num_syndromes"]),  # type: ignore[arg-type]
                "equal_coset_syndromes": int(res["equal_coset_syndromes"]),  # type: ignore[arg-type]
                "zero_coset_syndromes": int(res["zero_coset_syndromes"]),  # type: ignore[arg-type]
            }
        )

    # Trivial-syndrome coset weight enumerators (stored as [[w, count], ...] to stay strictly JSON-clean)
    a0, b0 = enumerators[0]
    cert: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "code_name": code.name,
        "family": code.family,
        "n": int(code.n),
        "k": 1,
        "r_X": int(meta["r_X"]),  # type: ignore[arg-type]
        "r_Z": int(meta["r_Z"]),  # type: ignore[arg-type]
        "H_X": [[int(x) for x in row] for row in code.H_X],
        "H_Z": [[int(x) for x in row] for row in code.H_Z],
        "L_X": [int(x) for x in code.L_X],
        "L_Z": [int(x) for x in code.L_Z],
        "L_Z_weight": int(meta["L_Z_weight"]),  # type: ignore[arg-type]
        "parity_branch": (
            "odd_lz_unitary_cptp"
            if meta["L_Z_parity"] == "odd"
            else "even_lz_hermitian_pauli_optimal"
        ),
        "swap_point": {
            "regime": str(swap_info["regime"]),
            "predicted_F_Pauli_pi4": canonical_fraction_str(
                swap_info["predicted_F_Pauli_pi4"]  # type: ignore[arg-type]
            ),
            "K_Z_dim": int(swap_info["K_Z_dim"]),  # type: ignore[arg-type]
            "K_Z_basis": [[int(b) for b in row] for row in swap_info["K_Z_basis"]],  # type: ignore[attr-defined]
            "immunity_witness_g_star": (
                [int(b) for b in swap_info["immunity_witness_g_star"]]  # type: ignore[attr-defined]
                if swap_info["immunity_witness_g_star"] is not None
                else None
            ),
        },
        "trivial_syndrome_enumerator": {
            "coset_0": [[int(w), int(c)] for w, c in sorted(a0.items())],
            "coset_1": [[int(w), int(c)] for w, c in sorted(b0.items())],
        },
        "total_error_count": 1 << code.n,
        "evaluations": evaluations,
    }
    _assert_strict_json_types(cert)
    return cert


def _check_bit_vector(vec: Any, n: int, name: str) -> Tuple[int, ...]:
    if type(vec) is not list or len(vec) != n:
        raise ValueError(f"{name} must be a list of length {n}")
    out: List[int] = []
    for idx, b in enumerate(vec):
        if type(b) is not int or b not in (0, 1):
            raise ValueError(f"{name}[{idx}] must be strict int 0 or 1, got {b!r}")
        out.append(b)
    return tuple(out)


def _check_bit_matrix(mat: Any, n: int, name: str) -> Tuple[Tuple[int, ...], ...]:
    if type(mat) is not list:
        raise ValueError(f"{name} must be a list of binary rows")
    return tuple(_check_bit_vector(row, n, f"{name}[{i}]") for i, row in enumerate(mat))


def verify_coherent_css_certificate_exact(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Independently verify a COHERENT-CSS-CERT-1.0 certificate with ZERO floating-point trust.

    Raises ValueError on ANY structural, algebraic, witness, type, or canonicality mismatch.
    Returns a summary dict upon full exact verification.
    """
    if type(cert) is not dict:
        raise ValueError("Certificate must be a JSON object (dict)")
    _assert_strict_json_types(cert)

    # 1. Schema and exact top-level key set
    if cert.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported schema_version {cert.get('schema_version')!r} != {SCHEMA_VERSION!r}"
        )
    if set(cert.keys()) != _TOP_KEYS:
        missing = sorted(_TOP_KEYS - set(cert.keys()))
        extra = sorted(set(cert.keys()) - _TOP_KEYS)
        raise ValueError(f"Top-level key mismatch: missing={missing}, extra={extra}")

    code_name = cert["code_name"]
    family = cert["family"]
    if type(code_name) is not str or not code_name:
        raise ValueError(f"code_name must be a non-empty str, got {code_name!r}")
    if type(family) is not str or not family:
        raise ValueError(f"family must be a non-empty str, got {family!r}")

    n = cert["n"]
    if type(n) is not int or n < 1 or n > 18:
        raise ValueError(f"Invalid n={n!r} (must be strict int in [1, 18])")
    if cert["k"] != 1:
        raise ValueError(f"Only k=1 CSS codes supported, got k={cert['k']!r}")
    if cert["total_error_count"] != (1 << n):
        raise ValueError(
            f"total_error_count {cert['total_error_count']} != 2^{n} ({1 << n})"
        )

    H_X = _check_bit_matrix(cert["H_X"], n, "H_X")
    H_Z = _check_bit_matrix(cert["H_Z"], n, "H_Z")
    L_X = _check_bit_vector(cert["L_X"], n, "L_X")
    L_Z = _check_bit_vector(cert["L_Z"], n, "L_Z")

    code = CSSCode(
        name=code_name,
        n=n,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L_X,
        L_Z=L_Z,
        family=family,
    )
    meta = validate_css_code(code)
    if not meta["all_z_stabilizers_even"]:
        raise ValueError("Code fails even-weight Z-stabilizer requirement")

    if cert["r_X"] != meta["r_X"] or cert["r_Z"] != meta["r_Z"]:
        raise ValueError(
            f"Rank mismatch: claimed ({cert['r_X']}, {cert['r_Z']}) vs actual ({meta['r_X']}, {meta['r_Z']})"
        )
    if cert["L_Z_weight"] != meta["L_Z_weight"]:
        raise ValueError(
            f"L_Z_weight mismatch: claimed {cert['L_Z_weight']} vs actual {meta['L_Z_weight']}"
        )

    expected_branch = (
        "odd_lz_unitary_cptp"
        if meta["L_Z_parity"] == "odd"
        else "even_lz_hermitian_pauli_optimal"
    )
    if cert["parity_branch"] != expected_branch:
        raise ValueError(
            f"parity_branch mismatch: claimed {cert['parity_branch']!r} vs actual {expected_branch!r}"
        )

    # 2. Verify Theorem 2.2 swap-point status, K_Z_basis, and witness g* independently
    swap_claimed = cert["swap_point"]
    if type(swap_claimed) is not dict:
        raise ValueError("swap_point must be a dict")
    if set(swap_claimed.keys()) != _SWAP_KEYS:
        missing = sorted(_SWAP_KEYS - set(swap_claimed.keys()))
        extra = sorted(set(swap_claimed.keys()) - _SWAP_KEYS)
        raise ValueError(f"swap_point key mismatch: missing={missing}, extra={extra}")

    swap_actual = compute_swap_point_status(code)
    if swap_claimed["regime"] != swap_actual["regime"]:
        raise ValueError(
            f"Swap-point regime mismatch: claimed {swap_claimed['regime']!r} "
            f"vs actual {swap_actual['regime']!r}"
        )
    pred_claimed = parse_canonical_fraction(
        swap_claimed["predicted_F_Pauli_pi4"], "predicted_F_Pauli_pi4"
    )
    if pred_claimed != swap_actual["predicted_F_Pauli_pi4"]:
        raise ValueError(
            f"predicted_F_Pauli_pi4 mismatch: claimed {pred_claimed} vs actual {swap_actual['predicted_F_Pauli_pi4']}"
        )
    if swap_claimed["K_Z_dim"] != swap_actual["K_Z_dim"]:
        raise ValueError(
            f"K_Z_dim mismatch: claimed {swap_claimed['K_Z_dim']} vs actual {swap_actual['K_Z_dim']}"
        )

    # Verify K_Z_basis is a valid linearly independent F_2 basis for exact K_Z = S_Z ∩ S_Z^⊥
    kz_claimed_rows = _check_bit_matrix(swap_claimed["K_Z_basis"], n, "swap_point.K_Z_basis")
    kz_claimed_masks = [bits_to_mask(r) for r in kz_claimed_rows]
    kz_actual_masks = compute_self_orthogonal_subcode(code)
    if len(kz_claimed_masks) != len(kz_actual_masks):
        raise ValueError(
            f"K_Z_basis row count {len(kz_claimed_masks)} != K_Z_dim {len(kz_actual_masks)}"
        )
    if gf2_rank(kz_claimed_masks, n) != len(kz_actual_masks):
        raise ValueError("Claimed K_Z_basis rows are not linearly independent over F_2")
    for idx, gm in enumerate(kz_claimed_masks):
        if not gf2_in_span(gm, kz_actual_masks, n):
            raise ValueError(f"Claimed K_Z_basis[{idx}] does not lie in K_Z = S_Z ∩ S_Z^⊥")
    for idx, gm in enumerate(kz_actual_masks):
        if not gf2_in_span(gm, kz_claimed_masks, n):
            raise ValueError(f"Actual K_Z vector {idx} not spanned by claimed K_Z_basis")

    hz_masks = code.hz_masks()
    lz_mask = code.lz_mask()
    witness = swap_claimed["immunity_witness_g_star"]
    if swap_actual["regime"] == "swap_immune_unit":
        if witness is None:
            raise ValueError("swap_immune_unit regime requires non-null immunity_witness_g_star")
        w_bits = _check_bit_vector(witness, n, "immunity_witness_g_star")
        w_mask = bits_to_mask(w_bits)
        if not gf2_in_span(w_mask, hz_masks, n):
            raise ValueError("Claimed immunity_witness_g_star does not lie in S_Z")
        if any(popcount(w_mask & rz) % 2 != 0 for rz in hz_masks):
            raise ValueError("Claimed immunity_witness_g_star is not self-orthogonal (not in S_Z^⊥)")
        if popcount(w_mask & lz_mask) % 2 != 1:
            raise ValueError("Claimed immunity_witness_g_star does not have odd overlap with L_Z")
    else:
        if witness is not None:
            raise ValueError("swap_collapse_half regime must have immunity_witness_g_star == null")

    # 3. Replay full coset weight enumerators and check trivial-syndrome enumerator
    enumerators = compute_coset_weight_enumerators(code)
    a0_actual, b0_actual = enumerators[0]
    triv_claimed = cert["trivial_syndrome_enumerator"]
    if type(triv_claimed) is not dict:
        raise ValueError("trivial_syndrome_enumerator must be a dict")
    if set(triv_claimed.keys()) != _TRIV_KEYS:
        missing = sorted(_TRIV_KEYS - set(triv_claimed.keys()))
        extra = sorted(set(triv_claimed.keys()) - _TRIV_KEYS)
        raise ValueError(f"trivial_syndrome_enumerator key mismatch: missing={missing}, extra={extra}")

    for key, actual_map in (("coset_0", a0_actual), ("coset_1", b0_actual)):
        pairs = triv_claimed[key]
        if type(pairs) is not list:
            raise ValueError(f"trivial_syndrome_enumerator.{key} must be a list of [w, count] pairs")
        parsed: Dict[int, int] = {}
        prev_w = -1
        for pair in pairs:
            if type(pair) is not list or len(pair) != 2:
                raise ValueError(f"Malformed [w, count] entry in {key}: {pair!r}")
            w, c = pair
            if type(w) is not int or type(c) is not int:
                raise ValueError(f"Non-integer [w, count] entry in {key}: {pair!r}")
            if w <= prev_w:
                raise ValueError(
                    f"Non-ascending or duplicate weight {w} (prev {prev_w}) in trivial_syndrome_enumerator.{key}"
                )
            if w < 0 or w > n or c <= 0:
                raise ValueError(f"Invalid weight/count ({w}, {c}) in {key}")
            prev_w = w
            parsed[w] = c
        if parsed != actual_map:
            raise ValueError(f"trivial_syndrome_enumerator.{key} mismatch: {parsed} != {actual_map}")

    # 4. Replay every rational evaluation point with ZERO floating-point ops
    evals = cert["evaluations"]
    if type(evals) is not list or len(evals) == 0:
        raise ValueError("evaluations must be a non-empty list")

    seen_p = set()
    has_pi4 = False
    for idx, ev in enumerate(evals):
        if type(ev) is not dict:
            raise ValueError(f"evaluations[{idx}] must be a dict")
        if set(ev.keys()) != _EVAL_KEYS:
            missing = sorted(_EVAL_KEYS - set(ev.keys()))
            extra = sorted(set(ev.keys()) - _EVAL_KEYS)
            raise ValueError(f"evaluations[{idx}] key mismatch: missing={missing}, extra={extra}")
        p = parse_canonical_fraction(ev["p"], f"evaluations[{idx}].p")
        if p < 0 or p > 1:
            raise ValueError(f"evaluations[{idx}].p out of [0, 1]: {p}")
        if p in seen_p:
            raise ValueError(f"Duplicate evaluation point p={p}")
        seen_p.add(p)

        actual = compute_exact_fidelities(code, p, enumerators=enumerators)
        for fld in ("F_Pauli", "F_CPTP", "F_Twirl", "p_trivial_syndrome"):
            val = parse_canonical_fraction(ev.get(fld), f"evaluations[{idx}].{fld}")
            if val != actual[fld]:
                raise ValueError(
                    f"evaluations[{idx}].{fld} mismatch at p={p}: claimed {val} != actual {actual[fld]}"
                )
        for ifld in ("num_syndromes", "equal_coset_syndromes", "zero_coset_syndromes"):
            if ev.get(ifld) != actual[ifld]:
                raise ValueError(
                    f"evaluations[{idx}].{ifld} mismatch at p={p}: claimed {ev.get(ifld)} != actual {actual[ifld]}"
                )

        # Cross-check Theorem 2.1 dichotomy on every point
        if expected_branch == "odd_lz_unitary_cptp":
            if actual["F_CPTP"] != Fraction(1, 1):
                raise ValueError("Odd-|L_Z| Theorem 2.1 invariant violated: F_CPTP != 1")
        else:
            if actual["F_CPTP"] != actual["F_Pauli"]:
                raise ValueError("Even-|L_Z| Theorem 2.1 invariant violated: F_CPTP != F_Pauli")

        # Cross-check Theorem 2.2 0-1 binary law at p = 1/2 (theta = pi/4)
        if p == Fraction(1, 2):
            has_pi4 = True
            if actual["F_Pauli"] != swap_actual["predicted_F_Pauli_pi4"]:
                raise ValueError(
                    f"Theorem 2.2 violation at p=1/2: F_Pauli={actual['F_Pauli']} "
                    f"!= predicted {swap_actual['predicted_F_Pauli_pi4']}"
                )
            if swap_actual["regime"] == "swap_collapse_half":
                if actual["equal_coset_syndromes"] != actual["num_syndromes"]:
                    raise ValueError("Swap-collapse syndrome equality failed at p=1/2")
            else:
                if actual["zero_coset_syndromes"] != actual["num_syndromes"]:
                    raise ValueError("Swap-immune vanishing-coset law failed at p=1/2")

    if not has_pi4:
        raise ValueError("Certificate must include mandatory swap-point evaluation p='1/2'")

    return {
        "verified": True,
        "code_name": code.name,
        "n": n,
        "parity_branch": expected_branch,
        "swap_regime": swap_actual["regime"],
        "num_evaluations": len(evals),
    }
