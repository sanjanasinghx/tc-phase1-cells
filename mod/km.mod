: Slow non-inactivating K (M) current, spike-frequency adaptation (Yamada et al. 1989; Pospischil 2008)
NEURON { SUFFIX km  USEION k READ ek WRITE ik  RANGE gbar, taumax }
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) }
PARAMETER { gbar = 7e-5 (S/cm2) taumax = 1000 (ms) }
STATE { p }
ASSIGNED { v (mV) ek (mV) ik (mA/cm2) p_inf tau_p (ms) }
BREAKPOINT { SOLVE states METHOD cnexp  ik = gbar*p*(v - ek) }
INITIAL { rates(v)  p = p_inf }
DERIVATIVE states { rates(v)  p' = (p_inf - p)/tau_p }
UNITSOFF
PROCEDURE rates(v(mV)) {
  p_inf = 1/(1 + exp(-(v + 35)/10))
  tau_p = taumax/(3.3*exp((v + 35)/20) + exp(-(v + 35)/20))
}
UNITSON
