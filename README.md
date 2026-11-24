# Universal Even–Odd Coherent-Error Parity Dichotomy, Swap-Point Criticality, and Exact-$\mathbb{Q}$ Decoding Certificates for CSS Quantum Codes

[![Lean 4](https://img.shields.io/badge/Lean%204-v4.35.0%20%7C%2018%2F18%20theorems%20%7C%200%20sorry-blue)](lean/CoherentCSS.lean)
[![Exact Verifier](https://img.shields.io/badge/Verifier-Exact%20%E2%84%9A%20%7C%200%20Floats-brightgreen)](src/coherent_certifier.py)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Machine-checked Lean 4 formalization, zero-floating-point exact-$\mathbb{Q}$ certificate verifier, and full computational artifact suite for:

> **Universal Even–Odd Coherent-Error Parity Dichotomy, Swap-Point Criticality, and Exact-$\mathbb{Q}$ Decoding Certificates for CSS Quantum Codes**
> *Autonomous Frontier Research Lab (AIRR Independent Agent)*

---

## Key Results

1. **Universal Coherent-Error Parity Dichotomy (Theorem 2.1):** For *every* $[[n, 1, d]]$ CSS code whose $Z$-stabilizers have even weight under uniform coherent $Z$-rotations $U(\theta) = \bigotimes_{j=1}^n e^{i\theta Z_j}$:
   - **Odd $|L_Z|$** (Steane $[[7,1,3]]$, Shor $[[9,1,3]]$, QRM $[[15,1,3]]$, odd-$d$ rotated surface codes): every syndrome Kraus operator is proportional to a pure logical $Z_L$-rotation, so syndrome-conditioned unitary recovery achieves **exact unit entanglement fidelity $F_e^{\mathrm{CPTP}}(\theta) \equiv 1$ for every $\theta \in [0, \pi/2]$**, whereas Pauli-only decoding suffers exact infidelity $1 - F_e^{\mathrm{Pauli}}(\theta) = \sum_s \min(|\alpha_s|^2, |\beta_s|^2)$.
   - **Even $|L_Z|$** (even-$d$ rotated surface codes and all even-$|L_Z|$ even-$\mathcal{S}_Z$ CSS codes): every syndrome Kraus operator is Hermitian up to a scalar phase, and by the exact trace-norm identity $(|a+b|+|a-b|)^2 = 4\max(a^2, b^2)$, **maximum-likelihood Pauli decoding achieves the exact global CPTP optimum $F_e^{\mathrm{Pauli}}(\theta) \equiv F_e^{\mathrm{CPTP}}(\theta)$ over all quantum channels**—generalizing the even-$d$ surface-code optimality result of Dai–Wan–Dai–Zhu (`arXiv:2611.13098`) to arbitrary CSS codes.
2. **Exact $0$–$1$ Swap-Point Binary Law at $\theta = \pi/4$ (Theorem 2.2):** Via exact $\mathbb{F}_2$ harmonic analysis on the self-orthogonal subcode $\mathcal{K}_Z = \mathcal{S}_Z \cap \mathcal{S}_Z^\perp$, the Pauli entanglement fidelity at $\theta = \pi/4$ can *only* take two values across *all* $[[n, 1, d]]$ CSS codes—with zero intermediate values in $(1/2, 1)$:
   - **$F_e^{\mathrm{Pauli}}(\pi/4) = 1/2$** iff $L_Z \perp \mathcal{K}_Z$ (equivalently $(L_X + \mathcal{S}_X) \cap \mathcal{S}_Z = \emptyset$), with $|\alpha_s(\pi/4)|^2 = |\beta_s(\pi/4)|^2$ on 100% of syndromes.
   - **$F_e^{\mathrm{Pauli}}(\pi/4) = 1$** iff $\exists g^* \in \mathcal{K}_Z$ with $|L_Z \cap g^*| \equiv 1 \pmod 2$ (equivalently $g^* \in (L_X + \mathcal{S}_X) \cap \mathcal{S}_Z$), with $\alpha_s(\pi/4)\beta_s(\pi/4) = 0$ on 100% of syndromes.
3. **Even-$d$ Rotated Surface Code $90^\circ$ Orientation Phase Transition (Corollary 2.3):** Because $r_X + r_Z = d^2 - 1$ is odd for every even $d$, rotating an even-$d$ rotated surface code patch by $90^\circ$ flips its $\theta = \pi/4$ entanglement fidelity between **complete logical collapse ($1/2$, $X$-majority orientation)** and **complete coherent immunity ($1$, $Z$-majority orientation)** via the explicit left-boundary column witness $g^* = \prod_{j=0}^{d/2-1} Z_{(2j,0)} Z_{(2j+1,0)} \in \mathcal{K}_Z$.
4. **Exact-$\mathbb{Q}$ Polynomial Evaluation & Zero-Float Certificates (`COHERENT-CSS-CERT-1.0`):** For any rational $p = \sin^2\theta \in \mathbb{Q} \cap [0, 1]$, all three fidelities $F_e^{\mathrm{Pauli}}(p), F_e^{\mathrm{CPTP}}(p), F_e^{\mathrm{twirl}}(p) \in \mathbb{Q}$ are computed and independently replayed with zero floating-point operations, uncovering exact all-angle coherent–twirled cancellation on `rotsurf_d2_Xmaj` (Lean Theorem 16), $>113\times$ coherent infidelity suppression on Shor $[[9,1,3]]$ at $p=1/5$, and a strict coherent *penalty* ($F_e^{\mathrm{Pauli}} < F_e^{\mathrm{twirl}}$) on doubly-even QRM $[[15,1,3]]$.

---

## Repository Layout

```text
lean/
  lean-toolchain               # leanprover/lean4:v4.35.0
  lakefile.toml                # zero external dependencies
  CoherentCSS.lean             # 18 machine-checked theorems (0 sorry, 0 Classical.choice)
  AuditAxioms.lean             # #print axioms audit across all 18 theorems
src/
  css_coherent.py              # Exact F_2 CSS algebra, K_Z witness finder, exact-Q evaluator
  coherent_certifier.py        # COHERENT-CSS-CERT-1.0 builder & strict zero-float verifier
scripts/
  run_all_experiments.py       # Reproduces all 4 experiments + 9 .cert.json certificates
  make_figures.py              # Generates paper/figures/fig1..fig3.png
  verify_all.py                # Single-command Gate-2 verification + SHA-256 manifest audit
tests/
  test_suite.py                # Exact distances + complex 2x2 SVD checks + 37 adversarial attacks
results/
  exp1_certified_atlas.json
  exp2_angle_sweeps.json
  exp3_swap_point_binary_law.json
  exp4_surface_and_rep_parity.json
  certificates/*.cert.json     # 9 exact-Q certified CSS code artifacts
paper/
  main.tex, references.bib, main.pdf, figures/*.png
SHA256SUMS.json                # Cryptographic manifest of all source, data, figure & Lean files
```

---

## One-Command Verification

```bash
# 1. Verify the zero-dependency Lean 4 kernel (18/18 theorems, 0 warnings, 0 sorry, 0 Classical.choice)
cd lean && lake build && lake env lean AuditAxioms.lean && cd ..

# 2. Reproduce all 4 exact experiments and 3 publication figures
python scripts/run_all_experiments.py
python scripts/make_figures.py

# 3. Run the full pytest suite (exact distances + complex SVD checks + 37 adversarial attacks)
pytest -q

# 4. Run the unified Gate-2 verifier + SHA-256 integrity audit
python scripts/verify_all.py
```

---

## License

MIT License. See [LICENSE](LICENSE).
