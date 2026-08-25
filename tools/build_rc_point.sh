#!/bin/bash
# Build aao_rc_point: links the tool source against all aao_rad physics
# sources except aao_rad.F itself (aao_rc_point.F is a modified copy of it).
cd "$(dirname "$0")"
R=../aao_rad
REPO=$(cd .. && pwd)
gfortran -fno-automatic -ffixed-line-length-none -fno-second-underscore \
  -funroll-loops -I$R -DMAID_TBL_DIR="\"$REPO\"" \
  aao_rc_point.F $R/aao.F $R/cgln_amps.F $R/daresbury.F $R/dsigma.F \
  $R/dvmpw.F $R/dvmpx.F $R/fint.F $R/helicity_amps.F $R/interp.F \
  $R/legendre.F $R/maid_lee.F $R/multipole_amps.F $R/read_sf_file.F \
  $R/sigmao.F $R/splie2.F $R/splin2.F $R/spline.F $R/splint.F \
  $R/strlen.F $R/timex.F $R/xsection.F $R/unixtime.c -o aao_rc_point \
  && echo "built: tools/aao_rc_point"
