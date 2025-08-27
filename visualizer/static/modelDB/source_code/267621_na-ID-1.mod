TITLE Fast Na current from C-S model, adaptable with shifts and scales

COMMENT
  Assumes the temp is 10 degC, wth rates corrected from C-S accordingly.
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
  SUFFIX na
  USEION na READ ena WRITE ina

  RANGE gbar, g, ina
  RANGE minf, hinf, taum, tauh,vshiftm,vscalem,tauscalem,vshifth,vscaleh,tauscaleh
  }

PARAMETER 
  {
  gbar =        0      (uS/mm2)
  vshiftm =    1.4763930000000001e+000  (mV)     : per params_07 in ours-na-kd-v8 
  vscalem =     9.1553499999999999e-001 
  tauscalem =   5.0015500000000002e-001
  vshifth =     -9.5136900000000004e+000  (mV)
  vscaleh =     5.0000000000000000e-001    
  tauscaleh =   2.0000000000000000e+000
  }

ASSIGNED
  {
  v (mV)
  ena (mV)
  minf
  taum (ms)
  hinf
  tauh (ms)
  g (uS/mm2)
  ina (mA/cm2)
  }

STATE
  {
  m
  h
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m^3*h
  ina = (1e-4)*(g*(v-ena))
  }

DERIVATIVE state_change
  {
  rates(v)  : Calculate minf, taum, hinf, tauh
  m' = (minf-m)/taum
  h' = (hinf-h)/tauh
  }

INITIAL 
  {
  rates(v)  : Calculate minf, taum, hinf, tauh
  m = minf
  h = hinf
  }

PROCEDURE rates(v(mV)) 
  {
  LOCAL Q, vscaledm, alpham, betam, vscaledh, alphah, betah
  Q=1.501533
  vscaledm=v/vscalem-vshiftm
  alpham=Q*1.0(1/ms)*linoid((vscaledm-(-35.0(mV)+5.3(mV)))/10.0(mV))/tauscalem
  betam=Q*4.0(1/ms)*exp(-(vscaledm-(-60.0(mV)+5.3(mV)))/18.0(mV))/tauscalem
  minf=alpham/(alpham+betam)
  taum=1.0/(alpham+betam)
  vscaledh=v/vscaleh-vshifth
  alphah=Q*0.07(1/ms)*exp(-(vscaledh-(-60.0(mV)+12.0(mV)))/20.0(mV))/tauscaleh
  betah=Q*1.0(1/ms)*lgc((vscaledh-(-30.0(mV)+12.0(mV)))/10.0(mV))/tauscaleh
  hinf=alphah/(alphah+betah)
  tauh=1.0/(alphah+betah)
  }

FUNCTION lgc(x) 
  {
  lgc = 1/(1+exp(-x))
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
  
