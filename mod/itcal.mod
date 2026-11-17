: T-type Ca current (Huguenard & McCormick 1992), ohmic form, with voltage shift for TRN
NEURON {
  SUFFIX itcal
  USEION ca READ eca WRITE ica
  RANGE gmax, vshift, m_inf, h_inf, tau_m, tau_h
}
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) }
PARAMETER { gmax = 0.002 (S/cm2) vshift = 0 (mV) celsius (degC) }
STATE { m h }
ASSIGNED { v (mV) eca (mV) ica (mA/cm2) m_inf h_inf tau_m (ms) tau_h (ms) phi_m phi_h }
BREAKPOINT { SOLVE states METHOD cnexp  ica = gmax*m*m*h*(v - eca) }
INITIAL { phi_m = 3.55^((celsius-24)/10)  phi_h = 2.8^((celsius-24)/10)  rates(v)  m = m_inf  h = h_inf }
DERIVATIVE states { rates(v)  m' = (m_inf - m)/tau_m  h' = (h_inf - h)/tau_h }
UNITSOFF
PROCEDURE rates(v(mV)) { LOCAL vv
  vv = v - vshift
  m_inf = 1/(1 + exp(-(vv + 57)/6.2))
  h_inf = 1/(1 + exp((vv + 81)/4))
  tau_m = (0.612 + 1/(exp(-(vv + 132)/16.7) + exp((vv + 16.8)/18.2)))/phi_m
  if (vv < -80) { tau_h = exp((vv + 467)/66.6)/phi_h } else { tau_h = (28 + exp(-(vv + 22)/10.5))/phi_h }
}
UNITSON
