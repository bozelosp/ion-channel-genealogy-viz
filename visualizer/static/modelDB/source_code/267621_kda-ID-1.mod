TITLE Kd current from C-S model, adaptable with shifts and scales

COMMENT
  Assumes the temp is 10 degC, and rates are corrected from C-S accordingly.
ENDCOMMENT

UNITS
  {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (millimolar)
  (S) = (mho)
  (uS) = (microS)
  }

NEURON  : this block manages the interface with NEURON
  {
  SUFFIX kda
  USEION k READ ek WRITE ik
  RANGE gbar, vshift, g, ik,vshiftm,vscalem,tauscalem
  }

 
: All var names used must appear in some variable-declaration block,
: even if they're mentioned in the NEURON block already

PARAMETER  : this is a variable-declaration block, for params settable
           : from the interface
  {
  gbar = 0 (uS/mm2)
  vshiftm =     2.7706100000000000e-001  (mV)     : per params_07 in ours-na-kd-v8 
  vscalem =     8.5349500000000000e-001
  tauscalem =   2.0000000000000000e+000
  }

ASSIGNED
  {
  v (mV)
  ek (mV)
  minf
  taum (ms)
  g (uS/mm2)
  ik (mA/cm2)
  }

STATE  : also a var-decl block, where one declares state vars local to
       : this particular mechanism
  {
  m
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m^4
  ik = (1e-4)*(g*(v-ek))
  }

DERIVATIVE state_change
  {
  rates(v)  : Calculate minf, taum
  m' = (minf-m)/taum
  }

INITIAL 
  {
  rates(v)  : Calculate minf, taum
  m = minf
  }

PROCEDURE rates(v(mV)) 
  {
  :TABLE minf, taum DEPEND vshift FROM -300 TO 300 WITH 601
  LOCAL Q, vscaledm, alpham, betam
  Q=1.501533 
  vscaledm=v/vscalem-vshiftm
  alpham=Q*(0.05(1/ms))*linoid((vscaledm-(-50.0(mV)+4.3(mV)))/(10.0(mV)))/tauscalem
  betam=Q*(0.0625(1/ms))*exp(-(vscaledm-(-60.0(mV)+4.3(mV)))/(80.0(mV)))/tauscalem
  minf=alpham/(alpham+betam)
  taum=1.0/(alpham+betam)
  }

FUNCTION linoid(x)
  {
  if ( (-1e-6<x) && (x<1e-6) )
    {
    linoid=1+x/2
    }
  else
    {
    linoid=x/(1-exp(-x))
    }
  }
