: High-voltage-activated Ca current (Reuveni et al. 1993), used in the L5 PT apical hot zone
NEURON { SUFFIX cahva  USEION ca READ eca WRITE ica  RANGE gbar }
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) }
PARAMETER { gbar = 0.001 (S/cm2) }
STATE { m h }
ASSIGNED { v (mV) eca (mV) ica (mA/cm2) ma mb ha hb }
BREAKPOINT { SOLVE states METHOD cnexp  ica = gbar*m*m*h*(v - eca) }
INITIAL { rates(v)  m = ma/(ma+mb)  h = ha/(ha+hb) }
DERIVATIVE states { rates(v)  m' = ma*(1-m) - mb*m  h' = ha*(1-h) - hb*h }
FUNCTION efun(x) { if (fabs(x) < 1e-6) { efun = 1 - x/2 } else { efun = x/(exp(x) - 1) } }
UNITSOFF
PROCEDURE rates(v(mV)) {
  ma = 0.055*3.8*efun((-27 - v)/3.8)
  mb = 0.94*exp((-75 - v)/17)
  ha = 0.000457*exp((-13 - v)/50)
  hb = 0.0065/(exp((-v - 15)/28) + 1)
}
UNITSON
