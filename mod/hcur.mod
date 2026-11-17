: Hyperpolarization-activated cation current Ih (McCormick & Pape 1990; Destexhe et al. 1993)
NEURON { SUFFIX hcur  NONSPECIFIC_CURRENT ih  RANGE ghbar, eh, m_inf, tau_m }
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) }
PARAMETER { ghbar = 2e-5 (S/cm2) eh = -40 (mV) celsius (degC) }
STATE { m }
ASSIGNED { v (mV) ih (mA/cm2) m_inf tau_m (ms) tadj }
BREAKPOINT { SOLVE states METHOD cnexp  ih = ghbar*m*(v - eh) }
INITIAL { tadj = 3^((celsius-36)/10)  rates(v)  m = m_inf }
DERIVATIVE states { rates(v)  m' = (m_inf - m)/tau_m }
UNITSOFF
PROCEDURE rates(v(mV)) {
  m_inf = 1/(1 + exp((v + 75)/5.5))
  tau_m = (1/(exp(-14.59 - 0.086*v) + exp(-1.87 + 0.0701*v)))/tadj
}
UNITSON
