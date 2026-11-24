"""Exact F_2 linear algebra and zero-float-trust exact-Q coherent decoding engine for CSS codes.

Implements:
  - Exact F_2 CSS validation (commutativity, ranks, k=1 logical pair, even-weight S_Z check).
  - Exact coset weight enumerators A_s(w), B_s(w) across all 2^{r_X} syndromes.
  - Exact self-orthogonal Z-stabilizer subcode K_Z = S_Z ∩ S_Z^⊥ and the exact
    Theorem 2.2 0-1 swap-point criticality criterion:
        F_e^{Pauli}(pi/4) = 1/2  <=>  L_Z ⊥ (S_Z ∩ S_Z^⊥)  <=>  (L_X + S_X) ∩ S_Z = ∅
        F_e^{Pauli}(pi/4) = 1    <=>  ∃ g* ∈ S_Z ∩ S_Z^⊥ with |L_Z ∩ g*| ≡ 1 (mod 2).
  - Zero-floating-point exact-Q evaluation of F_e^{Pauli}(p), F_e^{CPTP}(p), F_e^{twirl}(p)
    at any rational p = sin^2(theta) ∈ Q ∩ [0, 1].
  - Standard benchmark CSS code catalog: Steane [[7,1,3]], Shor [[9,1,3]], QRM [[15,1,3]],
    Rotated Surface Codes d=2, 3, 4, plus swap-immune CSS codes (e.g. [[6,1,2]] and [[11,1,3]]-like
    constructions with explicit witness g* ∈ (L_X + S_X) ∩ S_Z) and repetition-code family.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple


def bits_to_mask(bits: Sequence[int]) -> int:
    """Convert MSB-first bit list (length n) to integer bitmask."""
    mask = 0
    for b in bits:
        if b not in (0, 1):
            raise ValueError(f"Non-binary entry {b!r} in bit vector")
        mask = (mask << 1) | int(b)
    return mask


def mask_to_bits(mask: int, n: int) -> List[int]:
    """Convert integer bitmask to MSB-first bit list of length n."""
    return [(mask >> (n - 1 - i)) & 1 for i in range(n)]


def popcount(x: int) -> int:
    """Exact Hamming weight of a nonnegative integer bitmask."""
    return int(x).bit_count()


def gf2_rref(rows: Sequence[int], n: int) -> Tuple[List[int], List[int]]:
    """Compute reduced row echelon form over F_2 for bitmask rows of width n.

    Returns (basis_rows, pivot_bit_positions) where pivot_bit_positions are
    descending bit indices in [0, n-1].
    """
    mat = [int(r) for r in rows if int(r) != 0]
    basis: List[int] = []
    pivots: List[int] = []
    row_idx = 0
    for bit in range(n - 1, -1, -1):
        pivot_row = None
        for r in range(row_idx, len(mat)):
            if (mat[r] >> bit) & 1:
                pivot_row = r
                break
        if pivot_row is None:
            continue
        mat[row_idx], mat[pivot_row] = mat[pivot_row], mat[row_idx]
        p = mat[row_idx]
        for r in range(len(mat)):
            if r != row_idx and ((mat[r] >> bit) & 1):
                mat[r] ^= p
        pivots.append(bit)
        row_idx += 1
    basis = [r for r in mat[:row_idx] if r != 0]
    return basis, pivots


def gf2_rank(rows: Sequence[int], n: int) -> int:
    """Exact rank over F_2."""
    basis, _ = gf2_rref(rows, n)
    return len(basis)


def gf2_in_span(vec: int, rows: Sequence[int], n: int) -> bool:
    """Return True iff vec lies in rowspace_F2(rows)."""
    basis, pivots = gf2_rref(rows, n)
    v = int(vec)
    for row, bit in zip(basis, pivots):
        if (v >> bit) & 1:
            v ^= row
    return v == 0


def gf2_nullspace(rows: Sequence[int], n: int) -> List[int]:
    """Compute a basis for {x ∈ F_2^n : row · x ≡ 0 (mod 2) for all row in rows}."""
    basis, pivots = gf2_rref(rows, n)
    pivot_set = set(pivots)
    free_bits = [bit for bit in range(n - 1, -1, -1) if bit not in pivot_set]
    null_basis: List[int] = []
    for fbit in free_bits:
        vec = 1 << fbit
        for row, pbit in zip(basis, pivots):
            if (row >> fbit) & 1:
                vec |= 1 << pbit
        null_basis.append(vec)
    return null_basis


def gf2_subspace_intersection(rows_a: Sequence[int], rows_b: Sequence[int], n: int) -> List[int]:
    """Compute an F_2 basis for span(rows_a) ∩ span(rows_b) = (span(rows_a)^⊥ + span(rows_b)^⊥)^⊥."""
    null_a = gf2_nullspace(rows_a, n)
    null_b = gf2_nullspace(rows_b, n)
    return gf2_nullspace(null_a + null_b, n)


def gf2_enumerate_span(basis: Sequence[int]) -> List[int]:
    """Enumerate all 2^{len(basis)} vectors in the F_2 linear span of basis."""
    elems = [0]
    for b in basis:
        elems.extend(x ^ b for x in list(elems))
    return elems


@dataclass(frozen=True)
class CSSCode:
    """Exact specification of an [[n, 1, d]] CSS quantum code."""

    name: str
    n: int
    H_X: Tuple[Tuple[int, ...], ...]
    H_Z: Tuple[Tuple[int, ...], ...]
    L_X: Tuple[int, ...]
    L_Z: Tuple[int, ...]
    d_X: Optional[int] = None
    d_Z: Optional[int] = None
    family: str = "general"

    def hx_masks(self) -> List[int]:
        return [bits_to_mask(r) for r in self.H_X]

    def hz_masks(self) -> List[int]:
        return [bits_to_mask(r) for r in self.H_Z]

    def lx_mask(self) -> int:
        return bits_to_mask(self.L_X)

    def lz_mask(self) -> int:
        return bits_to_mask(self.L_Z)


def compute_exact_css_distances(code: CSSCode) -> Tuple[int, int]:
    """Compute exact (d_X, d_Z) by exhaustive enumeration over L_X ⊕ S_X and L_Z ⊕ S_Z."""
    hx = code.hx_masks()
    hz = code.hz_masks()
    lx = code.lx_mask()
    lz = code.lz_mask()
    if len(hx) > 16 or len(hz) > 16:
        raise ValueError("Exact distance enumeration requires r_X, r_Z <= 16")
    sx_elems = gf2_enumerate_span(hx)
    sz_elems = gf2_enumerate_span(hz)
    d_x = min(popcount(lx ^ g) for g in sx_elems)
    d_z = min(popcount(lz ^ g) for g in sz_elems)
    return d_x, d_z


def validate_css_code(code: CSSCode) -> Dict[str, object]:
    """Exhaustively validate [[n, 1, d]] CSS code invariants over F_2.

    Raises ValueError on any violation; returns structural metadata dict.
    """
    n = code.n
    if n < 1:
        raise ValueError("n must be >= 1")
    for label, mat in (("H_X", code.H_X), ("H_Z", code.H_Z)):
        for idx, row in enumerate(mat):
            if len(row) != n:
                raise ValueError(f"{label}[{idx}] length {len(row)} != n={n}")
            if any(b not in (0, 1) for b in row):
                raise ValueError(f"{label}[{idx}] contains non-binary entries")
    for label, vec in (("L_X", code.L_X), ("L_Z", code.L_Z)):
        if len(vec) != n:
            raise ValueError(f"{label} length {len(vec)} != n={n}")
        if any(b not in (0, 1) for b in vec):
            raise ValueError(f"{label} contains non-binary entries")

    hx = code.hx_masks()
    hz = code.hz_masks()
    lx = code.lx_mask()
    lz = code.lz_mask()

    # Commutativity H_X H_Z^T = 0 mod 2
    for i, rx in enumerate(hx):
        for j, rz in enumerate(hz):
            if popcount(rx & rz) % 2 != 0:
                raise ValueError(f"H_X[{i}] and H_Z[{j}] anticommute")

    # Full row rank and k = n - r_X - r_Z = 1
    r_x = gf2_rank(hx, n)
    r_z = gf2_rank(hz, n)
    if r_x != len(hx):
        raise ValueError(f"H_X rows are not linearly independent (rank {r_x} != {len(hx)})")
    if r_z != len(hz):
        raise ValueError(f"H_Z rows are not linearly independent (rank {r_z} != {len(hz)})")
    k = n - r_x - r_z
    if k != 1:
        raise ValueError(f"Expected k=1 logical qubit, got k={k} (n={n}, r_X={r_x}, r_Z={r_z})")

    # Logical operators
    for i, rx in enumerate(hx):
        if popcount(rx & lz) % 2 != 0:
            raise ValueError(f"L_Z anticommutes with H_X[{i}]")
    for j, rz in enumerate(hz):
        if popcount(rz & lx) % 2 != 0:
            raise ValueError(f"L_X anticommutes with H_Z[{j}]")
    if gf2_in_span(lz, hz, n):
        raise ValueError("L_Z lies in rowspace(H_Z)")
    if gf2_in_span(lx, hx, n):
        raise ValueError("L_X lies in rowspace(H_X)")
    if popcount(lx & lz) % 2 != 1:
        raise ValueError("L_X and L_Z do not anticommute")

    all_z_even = all(popcount(rz) % 2 == 0 for rz in hz)
    all_z_doubly_even = all(popcount(rz) % 4 == 0 for rz in hz) and all(
        popcount(hz[i] & hz[j]) % 2 == 0 for i in range(len(hz)) for j in range(i + 1, len(hz))
    )
    lz_wt = popcount(lz)
    lx_wt = popcount(lx)

    # If d_X / d_Z are declared and exact enumeration is tractable, assert exact match!
    if (code.d_X is not None or code.d_Z is not None) and r_x <= 14 and r_z <= 14:
        actual_dx, actual_dz = compute_exact_css_distances(code)
        if code.d_X is not None and code.d_X != actual_dx:
            raise ValueError(f"Declared d_X={code.d_X} != exact d_X={actual_dx} for {code.name}")
        if code.d_Z is not None and code.d_Z != actual_dz:
            raise ValueError(f"Declared d_Z={code.d_Z} != exact d_Z={actual_dz} for {code.name}")

    return {
        "n": n,
        "k": k,
        "r_X": r_x,
        "r_Z": r_z,
        "all_z_stabilizers_even": all_z_even,
        "all_z_stabilizers_doubly_even": all_z_doubly_even,
        "L_Z_weight": lz_wt,
        "L_Z_parity": "odd" if lz_wt % 2 == 1 else "even",
        "L_X_weight": lx_wt,
        "L_X_parity": "odd" if lx_wt % 2 == 1 else "even",
    }


def compute_self_orthogonal_subcode(code: CSSCode) -> List[int]:
    """Compute an F_2 basis for K_Z = S_Z ∩ S_Z^⊥."""
    hz = code.hz_masks()
    hz_perp = gf2_nullspace(hz, code.n)
    return gf2_subspace_intersection(hz, hz_perp, code.n)


def compute_swap_point_status(code: CSSCode) -> Dict[str, object]:
    """Compute the exact Theorem 2.2 0-1 swap-point classification at theta = pi/4.

    By Theorem 2.2 (exact F_2 harmonic analysis):
      F_e^{Pauli}(pi/4) = 1/2  <=>  L_Z ⊥ K_Z   (where K_Z = S_Z ∩ S_Z^⊥)
      F_e^{Pauli}(pi/4) = 1    <=>  ∃ g* ∈ K_Z with |L_Z ∩ g*| ≡ 1 (mod 2)
                                    (equivalently g* ∈ (L_X + S_X) ∩ S_Z).
    """
    kz_basis = compute_self_orthogonal_subcode(code)
    lz = code.lz_mask()
    witness_mask: Optional[int] = None
    for g in kz_basis:
        if popcount(g & lz) % 2 == 1:
            witness_mask = g
            break

    if witness_mask is None:
        return {
            "regime": "swap_collapse_half",
            "predicted_F_Pauli_pi4": Fraction(1, 2),
            "K_Z_dim": len(kz_basis),
            "K_Z_basis": [mask_to_bits(g, code.n) for g in kz_basis],
            "immunity_witness_g_star": None,
        }
    return {
        "regime": "swap_immune_unit",
        "predicted_F_Pauli_pi4": Fraction(1, 1),
        "K_Z_dim": len(kz_basis),
        "K_Z_basis": [mask_to_bits(g, code.n) for g in kz_basis],
        "immunity_witness_g_star": mask_to_bits(witness_mask, code.n),
    }


def compute_coset_weight_enumerators(
    code: CSSCode,
) -> Dict[int, Tuple[Dict[int, int], Dict[int, int]]]:
    """Compute exact coset weight enumerators (A_s, B_s) for all 2^{r_X} X-syndromes.

    For each syndrome s ∈ {0, ..., 2^{r_X}-1}, selects the canonical representative
    e_s = min{e ∈ {0, ..., 2^n-1} : H_X e = s} and partitions the syndrome fiber into:
      C_s^{(0)} = e_s ⊕ S_Z,    A_s(w) = |{e ∈ C_s^{(0)} : wt(e) = w}|,
      C_s^{(1)} = e_s ⊕ L_Z ⊕ S_Z,  B_s(w) = |{e ∈ C_s^{(1)} : wt(e) = w}|.
    """
    validate_css_code(code)
    n = code.n
    hx = code.hx_masks()
    hz = code.hz_masks()
    lz = code.lz_mask()
    r_x = len(hx)
    num_syn = 1 << r_x
    sz_elems = gf2_enumerate_span(hz)

    # Find canonical minimum-integer representative e_s for every syndrome s
    reps: Dict[int, int] = {}
    # Direct linear scan finds all 2^{r_X} representatives cleanly; fast for n <= 16.
    # For even faster exact construction at any n, solve H_X e = s via right-inverse
    # and then take min over the coset — or scan when 2^n <= 262144.
    if n <= 18:
        remaining = num_syn
        for e in range(1 << n):
            s = 0
            for rx in hx:
                s = (s << 1) | (popcount(rx & e) & 1)
            if s not in reps:
                reps[s] = e
                remaining -= 1
                if remaining == 0:
                    break
    else:
        raise ValueError(f"Exact full coset enumeration requires n <= 18 (got n={n})")

    if len(reps) != num_syn:
        raise RuntimeError(f"Found {len(reps)} syndromes, expected {num_syn}")

    enumerators: Dict[int, Tuple[Dict[int, int], Dict[int, int]]] = {}
    total_counted = 0
    for s in range(num_syn):
        e_s = reps[s]
        a_counts: Dict[int, int] = {}
        b_counts: Dict[int, int] = {}
        e_s_lz = e_s ^ lz
        for g in sz_elems:
            w0 = popcount(e_s ^ g)
            a_counts[w0] = a_counts.get(w0, 0) + 1
            w1 = popcount(e_s_lz ^ g)
            b_counts[w1] = b_counts.get(w1, 0) + 1
        total_counted += sum(a_counts.values()) + sum(b_counts.values())
        enumerators[s] = (dict(sorted(a_counts.items())), dict(sorted(b_counts.items())))

    if total_counted != (1 << n):
        raise RuntimeError(f"Coset partition sum {total_counted} != 2^{n}")

    return enumerators


def coset_abs_sq_exact(
    counts: Dict[int, int],
    n: int,
    p: Fraction,
) -> Fraction:
    """Compute exact |alpha_s(p)|^2 ∈ Q for an even-parity coset weight enumerator.

    Requires all nonzero-count weights in `counts` to share the same parity (guaranteed
    whenever all Z-stabilizers have even weight by Theorem 2.1 / Lean Theorem 2).
    """
    if not counts:
        return Fraction(0, 1)
    weights = list(counts.keys())
    parity0 = weights[0] & 1
    for w in weights:
        if (w & 1) != parity0:
            raise ValueError(
                f"Coset has mixed weight parities ({weights[0]} vs {w}); "
                "even-Z-stabilizer hypothesis violated"
            )
    q = Fraction(1, 1) - p
    # Precompute powers of p and q up to n
    p_pow = [Fraction(1, 1)] * (n + 1)
    q_pow = [Fraction(1, 1)] * (n + 1)
    for k in range(1, n + 1):
        p_pow[k] = p_pow[k - 1] * p
        q_pow[k] = q_pow[k - 1] * q

    total = Fraction(0, 1)
    items = list(counts.items())
    for i, (w1, c1) in enumerate(items):
        # Diagonal term (w1 == w2): sign = (-1)^0 = +1
        total += (c1 * c1) * q_pow[n - w1] * p_pow[w1]
        for j in range(i + 1, len(items)):
            w2, c2 = items[j]
            half_sum = (w1 + w2) // 2
            half_diff = (w1 - w2) // 2
            sign = 1 if (half_diff % 2 == 0) else -1
            # Factor 2 for ordered pairs (w1, w2) and (w2, w1) since (-1)^{(w2-w1)/2} = (-1)^{(w1-w2)/2}
            total += (2 * sign * c1 * c2) * q_pow[n - half_sum] * p_pow[half_sum]
    return total


def coset_twirl_prob_exact(counts: Dict[int, int], n: int, p: Fraction) -> Fraction:
    """Compute exact Pauli-twirled coset probability ∑_w A_s(w) (1-p)^{n-w} p^w ∈ Q."""
    q = Fraction(1, 1) - p
    total = Fraction(0, 1)
    for w, c in counts.items():
        total += c * (q ** (n - w)) * (p ** w)
    return total


def compute_exact_fidelities(
    code: CSSCode,
    p: Fraction,
    enumerators: Optional[Dict[int, Tuple[Dict[int, int], Dict[int, int]]]] = None,
) -> Dict[str, object]:
    """Compute exact-Q F_e^{Pauli}(p), F_e^{CPTP}(p), F_e^{twirl}(p), and per-syndrome stats.

    Works with ZERO floating-point operations for every even-Z-stabilizer [[n, 1, d]] CSS code
    and every rational p = sin^2(theta) ∈ Q ∩ [0, 1].
    """
    if not isinstance(p, Fraction) or p < 0 or p > 1:
        raise ValueError(f"p must be a Fraction in [0, 1], got {p!r}")
    meta = validate_css_code(code)
    if not meta["all_z_stabilizers_even"]:
        raise ValueError(f"Code {code.name} has odd-weight Z-stabilizers")

    if enumerators is None:
        enumerators = compute_coset_weight_enumerators(code)

    n = code.n
    lz_odd = meta["L_Z_parity"] == "odd"

    f_pauli = Fraction(0, 1)
    f_twirl = Fraction(0, 1)
    total_coherent_prob = Fraction(0, 1)
    total_twirl_prob = Fraction(0, 1)
    equal_coset_syndromes = 0
    zero_coset_syndromes = 0
    p_trivial = Fraction(0, 1)

    for s, (a_counts, b_counts) in enumerators.items():
        a_sq = coset_abs_sq_exact(a_counts, n, p)
        b_sq = coset_abs_sq_exact(b_counts, n, p)
        if a_sq < 0 or b_sq < 0:
            raise RuntimeError(f"Negative squared amplitude on syndrome {s}: {a_sq}, {b_sq}")
        total_coherent_prob += a_sq + b_sq
        f_pauli += a_sq if a_sq >= b_sq else b_sq
        if a_sq == b_sq:
            equal_coset_syndromes += 1
        if a_sq == 0 or b_sq == 0:
            zero_coset_syndromes += 1
        if s == 0:
            p_trivial = a_sq + b_sq

        ta = coset_twirl_prob_exact(a_counts, n, p)
        tb = coset_twirl_prob_exact(b_counts, n, p)
        total_twirl_prob += ta + tb
        f_twirl += ta if ta >= tb else tb

    if total_coherent_prob != Fraction(1, 1):
        raise RuntimeError(f"Coherent probability conservation failed: {total_coherent_prob} != 1")
    if total_twirl_prob != Fraction(1, 1):
        raise RuntimeError(f"Twirled probability conservation failed: {total_twirl_prob} != 1")

    # By Theorem 2.1 (machine-checked in Lean 4 Theorems 2-6):
    #   Odd |L_Z|  => K_s is pure logical Z-rotation => F_e^{CPTP}(p) = ∑_s (|α_s|^2 + |β_s|^2) = 1
    #   Even |L_Z| => K_s is Hermitian up to phase   => F_e^{CPTP}(p) = ∑_s max(|α_s|^2, |β_s|^2) = F_e^{Pauli}(p)
    f_cptp = Fraction(1, 1) if lz_odd else f_pauli

    return {
        "p": p,
        "F_Pauli": f_pauli,
        "F_CPTP": f_cptp,
        "F_Twirl": f_twirl,
        "p_trivial_syndrome": p_trivial,
        "num_syndromes": len(enumerators),
        "equal_coset_syndromes": equal_coset_syndromes,
        "zero_coset_syndromes": zero_coset_syndromes,
        "coherent_minus_twirl": f_pauli - f_twirl,
        "cptp_minus_pauli": f_cptp - f_pauli,
    }


# ============================================================================
# EXACT CSS CODE CONSTRUCTORS (BENCHMARK ATLAS)
# ============================================================================


def make_steane_7_1_3() -> CSSCode:
    """Steane [[7, 1, 3]] doubly-even self-dual CSS code."""
    H = (
        (1, 1, 1, 1, 0, 0, 0),
        (1, 1, 0, 0, 1, 1, 0),
        (1, 0, 1, 0, 1, 0, 1),
    )
    L = (1, 1, 1, 1, 1, 1, 1)
    return CSSCode(
        name="steane_7_1_3",
        n=7,
        H_X=H,
        H_Z=H,
        L_X=L,
        L_Z=L,
        d_X=3,
        d_Z=3,
        family="color_self_dual",
    )


def make_shor_9_1_3() -> CSSCode:
    """Shor [[9, 1, 3]] Bacon-Shor / concatenated phase-bit-flip CSS code (Z-checks wt 2, X-checks wt 6)."""
    H_Z = (
        (1, 1, 0, 0, 0, 0, 0, 0, 0),
        (0, 1, 1, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 1, 1, 0, 0, 0, 0),
        (0, 0, 0, 0, 1, 1, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 1, 1, 0),
        (0, 0, 0, 0, 0, 0, 0, 1, 1),
    )
    H_X = (
        (1, 1, 1, 1, 1, 1, 0, 0, 0),
        (0, 0, 0, 1, 1, 1, 1, 1, 1),
    )
    L_Z = (1, 0, 0, 1, 0, 0, 1, 0, 0)
    L_X = (1, 1, 1, 0, 0, 0, 0, 0, 0)
    return CSSCode(
        name="shor_9_1_3",
        n=9,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L_X,
        L_Z=L_Z,
        d_X=3,
        d_Z=3,
        family="concatenated",
    )


def make_qrm_15_1_3() -> CSSCode:
    """Quantum Reed-Muller [[15, 1, 3]] CSS code (transversal T; doubly-even S_Z of weight 8)."""
    # Columns 1..15 in 4-bit binary (1-indexed: x1, x2, x3, x4)
    cols = [((j >> 3) & 1, (j >> 2) & 1, (j >> 1) & 1, j & 1) for j in range(1, 16)]
    # Simplex code RM(1,4)* has 4 generators (weight 8) -> H_Z
    H_Z = tuple(tuple(cols[i][b] for i in range(15)) for b in range(4))
    # Punctured RM(2,4)* has 4 linear + 6 quadratic generators (weight 8 and 4) -> H_X
    quad_rows = []
    for b1 in range(4):
        for b2 in range(b1 + 1, 4):
            quad_rows.append(tuple(cols[i][b1] & cols[i][b2] for i in range(15)))
    H_X = H_Z + tuple(quad_rows)
    L = tuple(1 for _ in range(15))
    return CSSCode(
        name="qrm_15_1_3",
        n=15,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L,
        L_Z=L,
        d_X=3,
        d_Z=7,
        family="reed_muller",
    )


def make_rotated_surface_code(d: int, orientation: str = "X_majority") -> CSSCode:
    """Construct the exact [[d^2, 1, d]] rotated surface code for any d >= 2.

    For odd d, r_X = r_Z = (d^2 - 1)/2 and the two 90-degree orientations are isomorphic.
    For even d, r_X + r_Z = d^2 - 1 is ODD, breaking X <-> Z duality into two non-equivalent
    orientations under pure coherent Z-rotations:
      - orientation="X_majority" (r_X = d^2/2, r_Z = d^2/2 - 1):
          L_Z ⊥ K_Z => swap-point collapse F_e^{Pauli}(pi/4) = F_e^{CPTP}(pi/4) = 1/2
          (and at d=2, exact all-angle coherent-twirled cancellation, Lean Theorem 16).
      - orientation="Z_majority" (r_Z = d^2/2, r_X = d^2/2 - 1):
          The left boundary column g* = ∏_{j=0}^{d/2-1} Z_{(2j,0)} Z_{(2j+1,0)} = Z^{supp(X_L)}
          lies in K_Z = S_Z ∩ S_Z^⊥ with |L_Z ∩ g*| = 1 (odd!) => COMPLETE SWAP-POINT IMMUNITY
          F_e^{Pauli}(pi/4) = F_e^{CPTP}(pi/4) = 1!
    """
    if d < 2:
        raise ValueError(f"Rotated surface code requires d >= 2, got {d}")
    if orientation not in ("X_majority", "Z_majority"):
        raise ValueError(f"orientation must be 'X_majority' or 'Z_majority', got {orientation!r}")
    n = d * d

    def vec_from_coords(coords: Sequence[Tuple[int, int]]) -> Tuple[int, ...]:
        arr = [0] * n
        for r, c in coords:
            arr[r * d + c] = 1
        return tuple(arr)

    # First construct the Z-majority (bulk (r+c) even => X, odd => Z; left/right Z boundaries)
    # Note: for even d, this has r_Z = d^2/2, r_X = d^2/2 - 1.
    hx_list: List[Tuple[int, ...]] = []
    hz_list: List[Tuple[int, ...]] = []

    for r in range(d - 1):
        for c in range(d - 1):
            coords = [(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)]
            if (r + c) % 2 == 0:
                hx_list.append(vec_from_coords(coords))
            else:
                hz_list.append(vec_from_coords(coords))

    for c in range(d - 1):
        if (-1 + c) % 2 == 0:
            hx_list.append(vec_from_coords([(0, c), (0, c + 1)]))
        if (d - 1 + c) % 2 == 0:
            hx_list.append(vec_from_coords([(d - 1, c), (d - 1, c + 1)]))

    for r in range(d - 1):
        if (r - 1) % 2 == 1:
            hz_list.append(vec_from_coords([(r, 0), (r + 1, 0)]))
        if (r + d - 1) % 2 == 1:
            hz_list.append(vec_from_coords([(r, d - 1), (r + 1, d - 1)]))

    L_Z = vec_from_coords([(0, c) for c in range(d)])
    L_X = vec_from_coords([(r, 0) for r in range(d)])

    if orientation == "X_majority":
        # Swap X <-> Z (equivalent to 90-degree patch rotation):
        # now r_X = d^2/2, r_Z = d^2/2 - 1 for even d!
        hx_list, hz_list = hz_list, hx_list
        L_X, L_Z = L_Z, L_X

    suffix = "" if d % 2 == 1 else ("_Xmaj" if orientation == "X_majority" else "_Zmaj")
    code = CSSCode(
        name=f"rotsurf_d{d}{suffix}",
        n=n,
        H_X=tuple(hx_list),
        H_Z=tuple(hz_list),
        L_X=L_X,
        L_Z=L_Z,
        d_X=d,
        d_Z=d,
        family="rotated_surface",
    )
    validate_css_code(code)
    return code


def make_swap_immune_6_1_2() -> CSSCode:
    """Explicit [[6, 1, 2]] even-Z-stabilizer CSS code achieving EXACT F_e^{Pauli}(pi/4) = 1!

    Sharpness witness for Theorem 2.2:
      g* = (1,1,1,1,0,0) ∈ S_Z ∩ S_Z^⊥ is ALSO a logical-X operator (L_X = g*),
      with |L_Z ∩ g*| = |(1,0,0,0,1,0) ∩ (1,1,1,1,0,0)| = 1 (ODD!),
      forcing α_s(pi/4) * β_s(pi/4) = 0 on 100% of syndromes and F_e^{Pauli}(pi/4) = 1!
      Notice |L_Z| = 2 (EVEN!), so by Theorem 2.1 F_e^{CPTP}(θ) = F_e^{Pauli}(θ) for all θ,
      yet at θ = pi/4 it achieves UNIT fidelity (1) instead of collapsing to 1/2!
    """
    H_Z = (
        (1, 1, 1, 1, 0, 0),
        (0, 0, 0, 0, 1, 1),
    )
    # H_X must commute with H_Z (even weight on {0,1,2,3} and {4,5}), have rank 3,
    # and commute with L_Z = (1, 0, 0, 0, 1, 0) -> odd on {0} + odd on {4} => even overlap!
    H_X = (
        (0, 1, 1, 0, 0, 0),
        (0, 0, 1, 1, 0, 0),
        (1, 1, 0, 0, 1, 1),
    )
    L_X = (1, 1, 1, 1, 0, 0)  # Equals H_Z[0] ∈ S_Z ∩ S_Z^⊥!
    L_Z = (1, 0, 0, 0, 1, 0)  # Weight 2 (even!), |L_Z ∩ L_X| = 1
    code = CSSCode(
        name="swap_immune_6_1_2",
        n=6,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L_X,
        L_Z=L_Z,
        d_X=2,
        d_Z=2,
        family="swap_immune_witness",
    )
    validate_css_code(code)
    return code


def make_swap_immune_8_1_2_odd() -> CSSCode:
    """Explicit [[8, 1, 2]] (exact d_X=2, d_Z=3) doubly-even swap-immune CSS code
    with ODD |L_Z| = 3 and witness g* = (1,1,1,1,0,0,0,0) ∈ K_Z = S_Z ∩ S_Z^⊥.

    Every column of H_Z is nonzero (guaranteeing d_X >= 2, exact d_X = 2), and every
    element of L_Z ⊕ S_Z has weight 3 or 5 (exact d_Z = 3). Demonstrates swap-point
    immunity F_e^{Pauli}(pi/4) = 1 in the odd-|L_Z| branch with genuine distance d = 2!
    """
    H_Z = (
        (1, 1, 1, 1, 0, 0, 0, 0),
        (0, 0, 1, 1, 1, 1, 0, 0),
        (0, 0, 0, 0, 1, 1, 1, 1),
    )
    L_Z = (1, 0, 0, 0, 1, 0, 1, 0)  # Weight 3 (ODD!), exact d_Z = 3!
    L_X = (1, 1, 1, 1, 0, 0, 0, 0)  # Equals H_Z[0] ∈ S_Z ∩ S_Z^⊥, exact d_X = 2!
    hx_basis = gf2_nullspace(
        [bits_to_mask(r) for r in H_Z] + [bits_to_mask(L_Z)], 8
    )
    H_X = tuple(tuple(mask_to_bits(r, 8)) for r in hx_basis)
    code = CSSCode(
        name="swap_immune_8_1_2_odd",
        n=8,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L_X,
        L_Z=L_Z,
        d_X=2,
        d_Z=3,
        family="swap_immune_witness",
    )
    validate_css_code(code)
    return code


def make_repetition_z_code(d: int) -> CSSCode:
    """Phase-flip repetition code [[d, 1, d]]_Z (H_Z = empty? Wait: Z-checks for Z-errors = X-checks!).

    Wait: for pure Z-rotations U(theta) = exp(i theta ∑ Z_j), phase-flip detection uses X-stabilizers!
    Here H_X has (d-1) weight-2 X-checks X_j X_{j+1}, H_Z is empty (r_Z = 0, trivially all-even!),
    L_Z = Z_1 ... Z_d (weight d: ODD for odd d, EVEN for even d!), L_X = X_1!
    Every coset has size 1 (no Z-stabilizers => no coherent interference within cosets => F_Pauli == F_Twirl!),
    YET Theorem 2.1 STILL gives the exact even-odd parity dichotomy:
      - Odd d:  F_CPTP(theta) = 1 > F_Pauli(theta) = F_Twirl(theta) (reaches 1/2 at pi/4)!
      - Even d: F_CPTP(theta) = F_Pauli(theta) = F_Twirl(theta) (collapses to 1/2 at pi/4)!
    """
    if d < 2:
        raise ValueError("Repetition code requires d >= 2")
    H_X = tuple(
        tuple(1 if j in (i, i + 1) else 0 for j in range(d))
        for i in range(d - 1)
    )
    H_Z: Tuple[Tuple[int, ...], ...] = ()
    L_Z = tuple(1 for _ in range(d))
    L_X = tuple(1 if j == 0 else 0 for j in range(d))
    code = CSSCode(
        name=f"rep_phase_d{d}",
        n=d,
        H_X=H_X,
        H_Z=H_Z,
        L_X=L_X,
        L_Z=L_Z,
        d_X=1,
        d_Z=d,
        family="phase_repetition",
    )
    validate_css_code(code)
    return code


def build_code_catalog() -> Dict[str, CSSCode]:
    """Return the primary certified atlas of 9 CSS codes (including both even-d surface orientations)."""
    codes = [
        make_steane_7_1_3(),
        make_shor_9_1_3(),
        make_qrm_15_1_3(),
        make_rotated_surface_code(2, "X_majority"),
        make_rotated_surface_code(2, "Z_majority"),
        make_rotated_surface_code(3, "X_majority"),
        make_rotated_surface_code(4, "X_majority"),
        make_rotated_surface_code(4, "Z_majority"),
        make_swap_immune_8_1_2_odd(),
    ]
    return {c.name: c for c in codes}
