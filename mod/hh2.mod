: Traub-Miles fast Na / delayed-rectifier K (Destexhe & Pare 1999; Pospischil et al. 2008)
NEURON {
  SUFFIX hh2
  USEION na READ ena WRITE ina
  USEION k READ ek WRITE ik
  RANGE gnabar, gkbar, vtraub
}
UNITS { (mA) = (milliamp) (mV) = (millivolt) (S) = (siemens) }
PARAMETER {
  gnabar = 0.05 (S/cm2)
  gkbar  = 0.005 (S/cm2)
  vtraub = -60 (mV)
}
STATE { m h n }
ASSIGNED { v (mV) ena (mV) ek (mV) ina (mA/cm2) ik (mA/cm2)
  am bm ah bh an bn }
BREAKPOINT {
  SOLVE states METHOD cnexp
  ina = gnabar*m*m*m*h*(v - ena)
  ik  = gkbar*n*n*n*n*(v - ek)
}
INITIAL { rates(v) m = am/(am+bm) h = ah/(ah+bh) n = an/(an+bn) }
DERIVATIVE states {
  rates(v)
  m' = am*(1-m) - bm*m
  h' = ah*(1-h) - bh*h
  n' = an*(1-n) - bn*n
}
FUNCTION efun(x) { if (fabs(x) < 1e-6) { efun = 1 - x/2 } else { efun = x/(exp(x) - 1) } }
UNITSOFF
PROCEDURE rates(v(mV)) { LOCAL v2
  v2 = v - vtraub
  am = 0.32*4*efun((13 - v2)/4)
  bm = 0.28*5*efun((v2 - 40)/5)
  ah = 0.128*exp((17 - v2)/18)
  bh = 4/(1 + exp((40 - v2)/5))
  an = 0.032*5*efun((15 - v2)/5)
  bn = 0.5*exp((10 - v2)/40)
}
UNITSON
