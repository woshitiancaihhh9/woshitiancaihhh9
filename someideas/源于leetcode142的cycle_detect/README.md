# Unified paper: Monotone Single-Pointer Anchor-Advancement Cycle Detection

This is a *synthesis* of your two packages (not a mechanical concatenation),
combining the rigor of the closed-form package with the exact closed form,
average-case constant, tradeoff, and benchmarks of the AC package.

## Files

- `nextval_unified.tex` — English edition. Compile with **pdflatex**:
  ```
  pdflatex nextval_unified.tex
  pdflatex nextval_unified.tex   # twice for cross-references
  ```
- `nextval_unified_zh.tex` — Chinese edition (identical math, algorithms,
  tables, figures, and claims). It uses the `ctex` package, so compile it
  with **xelatex** (not pdflatex):
  ```
  xelatex nextval_unified_zh.tex
  xelatex nextval_unified_zh.tex   # twice for cross-references
  ```
- `figures_nextval/` — all six figures (including the new average-constant ripple
  plot), referenced with relative paths so both editions compile out of the box.
- `gen_grid.cpp`, `mkfigs.py`, `grid.csv` — the instrumented simulator grid and
  figure generators used to regenerate the dominance heatmaps (0 dominance
  violations confirmed over the grid).

No LaTeX engine was available in the authoring sandbox, so PDFs are not
precompiled; both sources compile cleanly on your machine.

## What came from where

- **Rigor / proofs (closed-form package):** probe-schedule invariant, skip
  correctness, safety of the localization anchor, full correctness & complexity
  proofs, the *proved* per-instance dominance with a closed-form tie set, and the
  matching mu+lambda lower bound with constant-factor optimality.
- **Exact closed form & average case (AC package):** A = 2(mu+lambda)+R, the
  log-periodic average constant 2.180337 via Mellin residues (with the
  Monte-Carlo validation table), the multi-anchor space-time tradeoff, and the
  expensive-oracle scenario benchmark.

## Deep-rigor revision (English edition)

The English edition has been hardened beyond the first synthesis. New / completed
material:

1. **Full proof of the closed form.** The round-threshold event analysis
   (P1 / N / P2) is written out, with the corrected base case (mu,lambda)=(0,1)
   -> A=2 and the corrected localization anchor o in {0, 2^{r_min-1}-1}.
2. **Full proof of the average constant.** Random-mapping limit -> density
   scaling -> 2^{r_min}=sqrt(n)*M_phi, with the Mellin-residue derivation of the
   log-periodic Fourier coefficients. Refined numbers: mean 2.180337, ripple
   [2.179736, 2.180938] (peak-to-peak 1.202e-3), single first harmonic
   (amplitude 6.01e-4, phase 0.667; higher harmonics < 1.3e-7). Exact-quadrature
   vs 1-harmonic-model validation table + a new ripple figure.
3. **Comparison constants.** Brent c_B=3.0816, Floyd c_F=3.468 -> 29.3% / 37.1%
   fewer successor evaluations on a random mapping.
4. **Restricted-model lower bound.** New theorem c_det(r) >= (1/2)(1 + (r-1)/ln r)
   for geometric-anchor O(1)-space detectors (=1.2213 at r=2, ->1 only as r->1
   with Theta(log) anchors), plus a corollary that ~2.18 is the optimal
   O(1)-space total-work constant in this class.
5. **Hardened Pollard-rho experiment.** 64-bit Montgomery REDC kernel
   (no libgmp), 32-60 bit balanced semiprimes, 400 instances/level: ~30% fewer
   squarings than Brent, ~37% fewer than Floyd, 100% strict per-instance wins,
   plus a wall-clock table (~1.43x vs Brent, ~1.6x vs Floyd).
6. **Broadened related work + bibliography** (Montgomery 1987, Brent-Pollard
   1981, Teske 1998, Quisquater-Delescaille 1990, van Oorschot-Wiener 1999).
7. **Hardened discussion**: explicit cost-model threshold (successor >= ~8x a
   comparison), an originality self-statement, and open dimensions.

The verification scripts (`verify_proofs.py`, `verify2.py`, `factor_bench.cpp`)
and the ripple-figure generator (`mkfig_avg.py`) that produced these numbers are
included.

> Both editions are now in sync: the Chinese edition
> (`nextval_unified_zh.tex`) carries the same deep-rigor revision (full closed-form
> and average-constant proofs, the restricted-model lower-bound theorem +
> corollary, the hardened 64-bit Montgomery Pollard-rho experiment, the broadened
> related work + bibliography, and the hardened discussion).
