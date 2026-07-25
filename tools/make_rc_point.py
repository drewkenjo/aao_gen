#!/usr/bin/env python3
"""
make_rc_point.py - regenerate aao_rc_point.F from ../aao_rad/aao_rad.F.

Run this after ANY physics change in aao_rad.F, then ./build_rc_point.sh.
The tool source is a surgical copy of the generator with six patches:
frozen kinematics, a parallel Born weight, fixed-tries termination and
a clean result printout. If an anchor is not found the script aborts -
that means aao_rad.F changed in a patched region and the patch below
must be reviewed by hand.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(HERE, "..", "aao_rad", "aao_rad.F")
c = open(src).read()

def rep(old, new, n=1):
    global c
    assert c.count(old) == n, f"anchor found {c.count(old)}x (expected {n}): {old[:60]!r}"
    c = c.replace(old, new)

rep("      real*8 sig_tot,sig_sum",
    """      real*8 sig_tot,sig_sum
      real*8 born_tot,sigb
      integer*8 ntry_max
      real csthcm_fix,phicm_fix""")

rep("      sig_tot	= 0.",
    """c     RC-POINT MODE: fixed hadronic CM point, fixed number of tries
      open(unit=15,file='point.inp',status='old')
      read(15,*) csthcm_fix
      read(15,*) phicm_fix
      read(15,*) ntry_max
      close(15)
      born_tot	= 0.d0
      write(6,*)' RC POINT MODE: csthcm,phicm,ntry =',
     &  csthcm_fix,phicm_fix,ntry_max

      sig_tot	= 0.""")

rep("""      csthcm	= -1.+2.*myran()
      phicm	= 360.*myran()""",
    """      csthcm	= csthcm_fix
      phicm	= phicm_fix""")

rep("""      ntries	= ntries+1""",
    """      ntries	= ntries+1
      sigb	= 0.d0
      if (ntries .gt. ntry_max) go to 50""")

rep("""	sigr	= sigr1/delta/4./pi""",
    """	sigr	= sigr1/delta/4./pi
	sigb	= signr/delta/4./pi""")

rep("""      sigr	= mcfac*mpfac*sigr
      sig_ratio	= sigr/sigr_max
      sig_tot	= sig_tot + sigr""",
    """      sigr	= mcfac*mpfac*sigr
      sig_ratio	= sigr/sigr_max
      sig_tot	= sig_tot + sigr
      sigb	= sigb*jacob*mcfac*mpfac
      born_tot	= born_tot + sigb
      go to 20""")

rep(""" 50   continue""",
    """ 50   continue
      write(6,*)' RC_POINT obs,born,ntries:',sig_tot,born_tot,ntries
      if (born_tot .gt. 0.d0) then
        write(6,*)' RC_POINT_RESULT ',sig_tot/born_tot
      endif
      stop""")

out = os.path.join(HERE, "aao_rc_point.F")
open(out, "w").write(c)
print(f"regenerated {out} from aao_rad.F")
