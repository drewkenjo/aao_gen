#!/usr/bin/env python3
"""
scan_phi.py - radiative correction RC = sigma_obs/sigma_Born at a FIXED
kinematic point (Q2, xB, t) as a function of phi, using aao_rc_point
(the Mo-Tsai machinery of aao_rad with frozen kinematics).

    ./scan_phi.py -E0 5.75 -Q2 1.125 -xb 0.1372 -t 0.12 -vcut 0.2
                  [-nphi 25] [-ntries 20000000] [-model 5] [-out rc_point_scan.txt]

Output: two columns (phi_deg, RC). Directly comparable point-by-point
with the EXCLURAD delta(phi) at the same (Q2, xB, t, vcut).

Caveat: the hadronic angles (cos theta*, phi) are fixed at the VERTEX,
while EXCLURAD fixes the OBSERVED proton t/phi; for hard photons inside
the vcut window the two differ slightly.
"""
import os, sys, argparse, subprocess

MP, MPI0 = 0.938272, 0.134976
HERE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("-E0",  type=float, required=True, help="beam energy (GeV)")
ap.add_argument("-Q2",  type=float, required=True)
ap.add_argument("-xb",  type=float, required=True)
ap.add_argument("-t",   type=float, required=True, help="|t| (GeV^2, positive)")
ap.add_argument("-vcut",type=float, default=0.2, help="mm_cut = half-width of the MM^2 window (GeV^2)")
ap.add_argument("-nphi",type=int,   default=25)
ap.add_argument("-ntries", type=int, default=20000000, help="MC tries per phi point")
ap.add_argument("-model",  type=int, default=5, help="physics model (5=MAID+DVMP)")
ap.add_argument("-out", type=str, default="rc_point_scan.txt")
a = ap.parse_args()

# point kinematics
W2 = MP*MP + a.Q2*(1.0/a.xb - 1.0)
W  = W2**0.5
nu = (W2 - MP*MP + a.Q2)/(2*MP)
Ep = a.E0 - nu
if Ep <= 0: sys.exit(f"ERROR: E' = {Ep:.3f} <= 0 at this point")
nu_cm  = (W2 - MP*MP - a.Q2)/(2*W)
E_pi   = (W2 + MPI0**2 - MP*MP)/(2*W)
p_pi   = (E_pi**2 - MPI0**2)**0.5
qv     = (((W2 + a.Q2 + MP*MP)/(2*W))**2 - MP*MP)**0.5
cscm   = (a.Q2 - MPI0**2 + 2*nu_cm*E_pi - a.t)/(2*qv*p_pi)
if abs(cscm) > 1: sys.exit(f"ERROR: cos(theta*) = {cscm:.3f} - point outside phase space")
print(f"point: W={W:.4f}  E'={Ep:.4f}  cos(theta*)={cscm:.5f}")

gen_inp = (f"{a.model}\n0\n.20 .12 .20 .20\n2\n1\n{a.vcut}\n0.0\n.43\n0.\n0.\n0.\n"
           f"{a.E0}\n{a.Q2} {a.Q2}\n{Ep:.5f} {Ep:.5f}\n.005\n1000\n0.\n1.0\n")
open("aao_point.inp","w").write(gen_inp)

with open(a.out,"w") as fo:
    fo.write(f"# aao_rc_point scan: E0={a.E0} Q2={a.Q2} xB={a.xb} t={a.t} vcut={a.vcut} W={W:.4f}\n")
    for i in range(a.nphi):
        phi = (i+0.5)*360.0/a.nphi
        open("point.inp","w").write(f"{cscm:.6f}\n{phi:.2f}\n{a.ntries}\n")
        r = subprocess.run([os.path.join(HERE,"aao_rc_point")],
                           stdin=open("aao_point.inp"), capture_output=True, text=True)
        rc = None
        for line in r.stdout.splitlines():
            if "RC_POINT_RESULT" in line:
                rc = float(line.split()[-1])
        if rc is None: sys.exit(f"ERROR at phi={phi}: no result\n{r.stdout[-500:]}")
        print(f"phi {phi:6.1f}  RC = {rc:.4f}")
        fo.write(f"{phi:8.2f} {rc:10.5f}\n")
print(f"-> {a.out}")
