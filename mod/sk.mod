: Small-conductance Ca-activated K current (Kohler et al. 1996, as in Hay et al. 2011)
NEURON { SUFFIX sk  USEION k READ ek WRITE ik  USEION ca READ cai  RANGE gbar }
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) (mM) = (milli/liter) }
PARAMETER { gbar = 0.001 (S/cm2) ztau = 1 (ms) }
STATE { z }
ASSIGNED { v (mV) ek (mV) ik (mA/cm2) cai (mM) zinf }
BREAKPOINT { SOLVE states METHOD cnexp  ik = gbar*z*(v - ek) }
INITIAL { rates(cai)  z = zinf }
DERIVATIVE states { rates(cai)  z' = (zinf - z)/ztau }
UNITSOFF
PROCEDURE rates(ca(mM)) {
  if (ca < 1e-7) { ca = 1e-7 }
  zinf = 1/(1 + (0.00043/ca)^4.8)
}
UNITSON
