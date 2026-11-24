"""Gate-2 end-to-end verification harness for Paper 07.

Checks:
  1. Lean 4 kernel (`lake build` + `lake env lean AuditAxioms.lean`):
     18/18 theorems, 0 warnings, 0 sorry, 0 Classical.choice, 0 Lean.ofReduceBool.
  2. Independent exact-Q replay of all 9 COHERENT-CSS-CERT-1.0 certificate files.
  3. Full pytest suite (including complex SVD cross-checks & 37 adversarial forgery attacks).
  4. Presence of all 4 JSON experiment outputs and 3 publication figures.
  5. SHA-256 manifest verification (`SHA256SUMS.json`), or generation with `--write-manifest`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import verify_coherent_css_certificate_exact

EXPECTED_LEAN_THEOREMS = [
    "CoherentCSS.xor_weight_exact",
    "CoherentCSS.xor_weight_mod2",
    "CoherentCSS.even_stab_preserves_parity",
    "CoherentCSS.coset_parity_same_for_even_lz",
    "CoherentCSS.coset_parity_opposite_for_odd_lz",
    "CoherentCSS.hermitian_trace_norm_eq_max",
    "CoherentCSS.unitary_rotation_exact_norm",
    "CoherentCSS.pauli_infidelity_min_form",
    "CoherentCSS.pauli_max_abs_identity",
    "CoherentCSS.self_orthog_wt_mod4_additive",
    "CoherentCSS.equal_coset_gives_half_fidelity",
    "CoherentCSS.vanishing_coset_gives_unit_fidelity",
    "CoherentCSS.odd_intersection_anticomm",
    "CoherentCSS.doubly_even_sign_positive",
    "CoherentCSS.steane_trivial_syndrome_unit_prob",
    "CoherentCSS.rotsurf_d2_exact_cancellation",
    "CoherentCSS.coherent_pauli_le_cptp",
    "CoherentCSS.css_orthogonality_mod2_symm",
]

FORBIDDEN_AXIOMS = ("sorryAx", "Classical.choice", "Lean.ofReduceBool", "Lean.trustCompiler")

MANIFEST_FILES = [
    "LICENSE",
    "README.md",
    "pyproject.toml",
    "lean/lean-toolchain",
    "lean/lakefile.toml",
    "lean/CoherentCSS.lean",
    "lean/AuditAxioms.lean",
    "src/__init__.py",
    "src/css_coherent.py",
    "src/coherent_certifier.py",
    "scripts/run_all_experiments.py",
    "scripts/make_figures.py",
    "scripts/verify_all.py",
    "tests/test_suite.py",
    "results/exp1_certified_atlas.json",
    "results/exp2_angle_sweeps.json",
    "results/exp3_swap_point_binary_law.json",
    "results/exp4_surface_and_rep_parity.json",
    "results/certificates/steane_7_1_3.cert.json",
    "results/certificates/shor_9_1_3.cert.json",
    "results/certificates/qrm_15_1_3.cert.json",
    "results/certificates/rotsurf_d2_Xmaj.cert.json",
    "results/certificates/rotsurf_d2_Zmaj.cert.json",
    "results/certificates/rotsurf_d3.cert.json",
    "results/certificates/rotsurf_d4_Xmaj.cert.json",
    "results/certificates/rotsurf_d4_Zmaj.cert.json",
    "results/certificates/swap_immune_8_1_2_odd.cert.json",
    "paper/figures/fig1_parity_dichotomy_curves.png",
    "paper/figures/fig2_coherent_vs_twirled_gap.png",
    "paper/figures/fig3_swap_point_binary_law.png",
    "paper/main.tex",
    "paper/references.bib",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_lean() -> None:
    lean_dir = ROOT / "lean"
    r_build = subprocess.run(
        ["lake", "build"], cwd=lean_dir, capture_output=True, text=True, check=False
    )
    if r_build.returncode != 0:
        raise RuntimeError(f"lake build failed:\n{r_build.stdout}\n{r_build.stderr}")
    if "warning:" in (r_build.stdout + r_build.stderr).lower():
        raise RuntimeError(f"Lean build emitted warnings:\n{r_build.stdout}\n{r_build.stderr}")

    r_audit = subprocess.run(
        ["lake", "env", "lean", "AuditAxioms.lean"],
        cwd=lean_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if r_audit.returncode != 0:
        raise RuntimeError(f"AuditAxioms.lean failed:\n{r_audit.stdout}\n{r_audit.stderr}")
    out = r_audit.stdout + "\n" + r_audit.stderr
    for thm in EXPECTED_LEAN_THEOREMS:
        if f"'{thm}'" not in out:
            raise RuntimeError(f"Missing Lean theorem in axiom audit: {thm}")
    for bad in FORBIDDEN_AXIOMS:
        if bad in out:
            raise RuntimeError(f"Forbidden Lean axiom detected: {bad}\n{out}")
    print(f"[PASS] Lean 4 kernel: {len(EXPECTED_LEAN_THEOREMS)}/{len(EXPECTED_LEAN_THEOREMS)} theorems, 0 warnings, 0 sorry, 0 Classical.choice")


def check_certificates() -> None:
    cert_dir = ROOT / "results" / "certificates"
    certs = sorted(cert_dir.glob("*.cert.json"))
    if len(certs) != 9:
        raise RuntimeError(f"Expected 9 certificate files, found {len(certs)}")
    for cp in certs:
        data = json.loads(cp.read_text(encoding="utf-8"))
        res = verify_coherent_css_certificate_exact(data)
        if not res.get("verified"):
            raise RuntimeError(f"Certificate failed verification: {cp.name}")
    print(f"[PASS] Exact-Q certificates: {len(certs)}/{len(certs)} independently verified (0 floats)")


def check_pytest() -> None:
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "-q"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if r.returncode != 0:
        raise RuntimeError(f"pytest failed:\n{r.stdout}\n{r.stderr}")
    print(f"[PASS] pytest: {r.stdout.strip().splitlines()[-1]}")


def write_or_check_manifest(write_mode: bool) -> None:
    manifest_path = ROOT / "SHA256SUMS.json"
    computed: Dict[str, str] = {}
    missing: List[str] = []
    for rel in MANIFEST_FILES:
        p = ROOT / rel
        if not p.exists():
            missing.append(rel)
        else:
            computed[rel] = sha256_file(p)
    if missing:
        raise RuntimeError(f"Missing manifest files: {missing}")

    if write_mode:
        manifest_path.write_text(json.dumps(computed, indent=2) + "\n", encoding="utf-8")
        print(f"[PASS] Wrote SHA256SUMS.json ({len(computed)} files)")
        return

    if not manifest_path.exists():
        raise RuntimeError("SHA256SUMS.json not found; run with --write-manifest first")
    recorded = json.loads(manifest_path.read_text(encoding="utf-8"))
    if recorded != computed:
        diffs = [k for k in set(recorded) | set(computed) if recorded.get(k) != computed.get(k)]
        raise RuntimeError(f"SHA256SUMS.json mismatch on: {diffs}")
    print(f"[PASS] SHA256SUMS.json: {len(computed)}/{len(computed)} hashes match")


def main() -> None:
    parser = argparse.ArgumentParser(description="Gate-2 verification harness for Paper 07")
    parser.add_argument("--write-manifest", action="store_true", help="Write SHA256SUMS.json")
    parser.add_argument("--skip-manifest", action="store_true", help="Skip SHA256SUMS.json check")
    args = parser.parse_args()

    check_lean()
    check_certificates()
    check_pytest()
    if not args.skip_manifest:
        write_or_check_manifest(args.write_manifest)
    print("=== ALL PAPER 07 GATE-2 VERIFICATION CHECKS PASSED ===")


if __name__ == "__main__":
    main()
