# aao_gen — exclusive single-pion electroproduction generators

This repository contains two Fortran event generators for exclusive
pseudoscalar-meson electroproduction off the proton
(`e p -> e' p pi0`, `e p -> e' n pi+`, `e p -> e' p eta`),
based on the AO and MAID physics models:

| Generator | Description |
|---|---|
| `aao_norad` | Born-level (non-radiative) cross section |
| `aao_rad`   | Same physics model **plus internal radiative effects** (vertex corrections, vacuum polarization, real-photon emission) in the Mo–Tsai formalism |

Both generators use accept–reject sampling, so the output events are
**unweighted**. At the end of the run each generator prints the
**integrated cross section** over the generation region (within all input
limits and cuts) — this number is required for absolute normalization and
for radiative-correction calculations (see
[tex/rad_corrections_phi.pdf](tex/rad_corrections_phi.pdf)).

## Building

```bash
make            # builds both generators
```

Executables appear in `aao_rad/build/aao_rad` and
`aao_norad/build/aao_norad`.

## Running

Both programs read their parameters from standard input:

```bash
cd aao_rad
./build/aao_rad   < aao_input.inp

cd aao_norad
./build/aao_norad < aao_norad_input.inp
```

Output files:

| File | Content |
|---|---|
| `aao_rad.lund` / `aao_norad.lund` | events in LUND format (input for GEMC) |
| `aao_rad.out`  / `aao_norad.out`  | run log, including the integrated cross section |
| `aao_rad.sum`  / `aao_norad.sum`  | run summary (input parameters, integrated cross section) |

---

## Input file: `aao_rad`

Parameters are read **line by line, in this exact order**
(comments after `!` are allowed):

```
1                 ! (1)  th_opt: physics model. 1=AO, 4=MAID98, 5=MAID2000
1                 ! (2)  flag_ehel: 1=polarized electron, 0=unpolarized
.20 .12 .20 .20   ! (3)  reg1..reg4: sizes of the cos(theta_k) MC integration
                  !      regions (sum must be < 0.95)
4                 ! (4)  npart: 2 = e- and charged hadron only;
                  !      4 = also write the neutral hadron and radiated photon
1                 ! (5)  epirea: final state. 1=pi0, 2=pi-, 3=pi+, 5=eta
.2                ! (6)  mm_cut: half-width of the missing-mass-squared window
                  !      (GeV^2) -- SEE THE IMPORTANT SECTION BELOW
4.0               ! (7)  t_targ: target cell length (cm), used for external
                  !      radiation in the target material
.43               ! (8)  r_targ: target radius (cm)
.082              ! (9)  vertex_x: x-coordinate of the beam position (cm)
-.449             ! (10) vertex_y: y-coordinate of the beam position (cm)
-.139             ! (11) vertex_z: z-coordinate of the beam position (cm)
10.6              ! (12) ebeam: incident electron energy (GeV)
0.3 10.6          ! (13) q2_min q2_max: Q^2 generation limits (GeV^2)
0.3 10.6          ! (14) ep_min ep_max: scattered-electron energy limits (GeV)
.005              ! (15) delta: minimum photon energy for integration (GeV).
                  !      Soft-photon cutoff: below delta, emission is
                  !      integrated analytically and exponentiated. Keep it
                  !      small (~5 MeV); physics results must not depend on it.
                  !      This is NOT an analysis cut.
100000            ! (16) nmax: number of events to generate
1.0               ! (17) fmcall: multiplication factor for sigr_max
                  !      (accept-reject ceiling). If 0, the next line is read:
[.005]            ! (18) sigr_max: explicit maximum cross section
                  !      (ONLY read if fmcall = 0)
[1.1]             ! (19) w_min_input: minimum W (GeV)
                  !      (ONLY read if th_opt > 10)
```

### The `mm_cut` parameter (line 6) — read this before computing radiative corrections

`mm_cut` is a cut applied **inside the generator, during event generation**,
on the *experimental* missing mass squared:

```fortran
if (abs(mm2 - mm_exp) .gt. mm_cut) go to 20     ! reject, regenerate
sig_tot = sig_tot + sigr                        ! only surviving events count
```

- `mm2` is computed the way an experimentalist computes it — from the
  nominal beam energy, the scattered electron, and the detected hadron
  (proton for pi0/eta, pion for pi+), **ignoring the radiated photon**.
  For the pi0 channel it is MM²(ep → e′pX). When a hard photon is emitted,
  `mm2` shifts away from its nominal value: this shift is the radiative tail.
- `mm_exp` is the nominal value: m²(π⁰) for `epirea=1`, m²(p)≈m²(n) for
  `epirea=3`, m²(η) for `epirea=5`.
- Events with |mm2 − mm_exp| > `mm_cut` are rejected **before** the cross
  section is accumulated. Therefore **the integrated cross section printed
  at the end of the run is the cross section *within* the mm_cut window.**

Consequences for a radiative-correction analysis:

1. The radiative correction is only defined *with respect to a specific
   missing-mass cut*. `mm_cut` plays exactly the role of the `v_cut`
   (inelasticity cut) in other RC codes.
2. `mm_cut` must **match, or be looser than**, the MM² cut applied to the
   data. If it is looser, apply the exact analysis cut offline (after
   detector simulation, on the reconstructed MM², same code as for data) —
   the normalization stays correct because
   dσ = σ_printed × n_pass / N_generated.
3. **Never** set `mm_cut` tighter than the analysis cut on data: the part of
   the radiative tail present in the data would then be missing from the
   Monte Carlo, and the radiative correction would be underestimated with no
   way to recover it offline.
4. `aao_norad` has no radiated photon, so at generator level its MM² is
   exact and no such parameter exists there. After detector smearing, apply
   the same reconstructed-MM² cut to both samples.

Note: for `epirea=2` (pi-) the reference value `mm_exp` is not set in the
code — use the pi0, pi+ or eta channels when the missing-mass cut matters.

### Example input files

Ready-to-use examples in `aao_rad/`: `aao_input.inp`,
`aao-pi0-1.6.inp`, `aao-pi0-2.4.inp`, `amaid-pi0-1.5.inp`,
`amaid-pi0-1.6.inp`, `amaid-pi0-2.4.inp`.

---

## Input file: `aao_norad`

All parameters are read with a single list-directed read, in this order
(free format — one value per line works):

```
5                 ! (1)  phys: physics model. 1=AO, 4=MAID98, 5=MAID2000
1                 ! (2)  flag_ehel: 1=polarized electron, 0=unpolarized
3                 ! (3)  npart: 2 = e- and charged hadron; 3 = also the
                  !      neutral hadron (e-, h+, h0)
1                 ! (4)  epirea: final state. 1=pi0, 2=pi-, 3=pi+
10.6              ! (5)  ebeam: incident electron energy (GeV)
0.2 10.6          ! (6,7) q2_min q2_max: Q^2 generation limits (GeV^2)
0.2 10.6          ! (8,9) ep_min ep_max: scattered-electron energy limits (GeV)
10000             ! (10) nmax: number of events to generate
1.0               ! (11) fmcall: scale factor for the accept-reject ceiling
1                 ! (12) boso: 1=BOS output, 0=no BOS output
0                 ! (13) seed_source: 0 = seed from machine time;
                  !      any other value is used as the random seed
```

---

## Radiative corrections

To compute the radiative correction factor as a function of a kinematic
variable (e.g. φ):

1. Generate N events with each generator using **identical generation
   limits** (Q², E′ ranges) and identical analysis cuts.
2. Record the two integrated cross sections σ_rad and σ_norad printed at
   the end of the runs.
3. Fill the norad histogram with weight 1 and the rad histogram with the
   constant weight w = σ_rad/σ_norad (for equal N); then

   RC(φ) = (σ_rad/σ_norad) · h_rad(φ)/h_norad(φ),
   σ_Born = σ_exp / RC.

For radiated events, reconstruct the kinematics from the scattered electron
(as in data), not from the true hadronic-vertex values. Full details,
formulas and caveats: [tex/rad_corrections_phi.pdf](tex/rad_corrections_phi.pdf).

## Python wrapper

`aao_gen` / `gen_wrapper/` contain a Python wrapper (`aao_gen.py`) that
drives `aao_norad` with command-line options and additional kinematic
filtering (xB, W², t limits). Run `./aao_gen.py -h` for options. The
radiative generator is **not** supported by the wrapper; run
`aao_rad` directly as described above.
