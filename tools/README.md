# tools/aao_rc_point — radiative correction at a fixed kinematic point

Computes RC = sigma_obs/sigma_Born of aao_rad's Mo-Tsai machinery at a
FIXED point (Q2, xB or W, t, phi): no event generation, no binning, no
migration. Built for point-by-point comparison with EXCLURAD delta(phi).

## Usage
```
./build_rc_point.sh
./scan_phi.py -E0 5.75 -Q2 1.125 -xb 0.1372 -t 0.12 -vcut 0.2
```
Output: rc_point_scan.txt (phi, RC). ~17 s per phi point at 2e7 tries.

## How it works
aao_rc_point.F is a surgical copy of aao_rad.F (regenerate with
make_rc_point.py after any physics change in aao_rad.F!) with:
- electron kinematics frozen via degenerate Q2/E' input ranges;
- hadronic CM angles (cos theta*, phi) frozen from point.inp;
- on every MC try TWO weights accumulate through the identical code
  path: the observed one (with radiative factors) and a Born one
  (soft branch, radiative factors set to 1; the soft sampling measure
  integrates to unity). RC = sum(obs)/sum(born): all jacobians,
  region factors and common constants cancel exactly by construction;
- event building/output is bypassed; fixed number of tries.

## Caveats
- Hadron angles are fixed at the VERTEX; EXCLURAD fixes the OBSERVED
  proton t/phi. Within the vcut window the difference is small but
  nonzero for hard photons.
- 2026-07 finding obtained with this tool: the factorized Mo-Tsai tail
  gives a phi-dependence of RC ~10x flatter than the exact
  Bardin-Shumeiko integral (EXCLURAD), and sits ~10% below it at the
  plateau of the 5.75 GeV pi0 test point. Use EXCLURAD for RC numbers;
  use aao_rad for event samples (acceptance, MM2 shapes).
