: Submembrane Ca accumulation with first-order decay (Destexhe et al. 1994)
NEURON { SUFFIX cadyn  USEION ca READ ica WRITE cai  RANGE depth, taur, cainf, gamma }
UNITS { (mA) = (milliamp) (mM) = (milli/liter) (um) = (micron) FARADAY = (faraday) (coulomb) }
PARAMETER { depth = 0.1 (um) taur = 80 (ms) cainf = 1e-4 (mM) gamma = 0.05 }
STATE { cai (mM) }
ASSIGNED { ica (mA/cm2) }
INITIAL { cai = cainf }
BREAKPOINT { SOLVE state METHOD cnexp }
DERIVATIVE state { cai' = -(10000)*gamma*ica/(2*FARADAY*depth) - (cai - cainf)/taur }
