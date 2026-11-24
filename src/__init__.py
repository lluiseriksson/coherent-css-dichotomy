"""Paper 07: Universal Even-Odd Coherent-Error Parity Dichotomy, Swap-Point Criticality,
and Exact-Q Decoding Certificates for CSS Quantum Codes."""

from .css_coherent import (
    CSSCode,
    build_code_catalog,
    compute_coset_weight_enumerators,
    compute_exact_fidelities,
    compute_self_orthogonal_subcode,
    compute_swap_point_status,
    validate_css_code,
)
from .coherent_certifier import (
    SCHEMA_VERSION,
    build_coherent_css_certificate,
    canonical_fraction_str,
    parse_canonical_fraction,
    verify_coherent_css_certificate_exact,
)

__all__ = [
    "CSSCode",
    "SCHEMA_VERSION",
    "build_code_catalog",
    "build_coherent_css_certificate",
    "canonical_fraction_str",
    "compute_coset_weight_enumerators",
    "compute_exact_fidelities",
    "compute_self_orthogonal_subcode",
    "compute_swap_point_status",
    "parse_canonical_fraction",
    "validate_css_code",
    "verify_coherent_css_certificate_exact",
]
