TITLE I_A current from Connor-Stevens model, adaptable with shifts and scales

COMMENT
  Assumes the temp is 10 degC, and rates are corrected from CS accordingly.
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
NEURON
  {
  SUFFIX aa
  USEION k READ ek WRITE ik
  RANGE gbar, ik, minf, hinf, taum, tauh, vshiftm,vscalem,tauscalem,vshifth,vscaleh,tauscaleh
  }

PARAMETER  : this is a variable-declaration block, for params settable
           : from the interface
  {
  gbar =        0      (uS/mm2)
  vshiftm =     -4.6005120000000002e+000  (mV)     : per params_07 in ours-na-kd-v8 
  vscalem =     1.0363050000000000e+000 
  tauscalem =   1.8606229999999999e+000
  vshifth =    -8.6378999999999995e-001  (mV)
  vscaleh =     1.0326500000000001e+000    
  tauscaleh =   2.0000000000000000e+000
  }

ASSIGNED
  {
  v (mV)
  ek (mV)
  minf
  taum (ms)
  hinf
  tauh (ms)
  g (uS/mm2)
  ik (mA/cm2)
  }

STATE  : also a var-decl block, where one declares state vars local to
       : this particular mechanism
  {
  m
  h
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m^3*h
  ik = (1e-4)*(g*(v-ek))
  }

DERIVATIVE state_change
  {
  rates(v)  : Calculate minf, taum, hinf, tauh
  m' = (minf-m)/taum
  h' = (hinf-h)/tauh
  }

INITIAL 
  {
  rates(v)
  m = minf
  h = hinf
  }

PROCEDURE rates(v(mV)) 
  {
  LOCAL Q, vscaledm, vscaledh
  Q=1.501533
  vscaledm=v/vscalem-vshiftm
  minf=(0.0761*exp((vscaledm-(-94.22(mV)))/(31.84(mV)))*lgc(-(vscaledm-(-1.17(mV)))/(28.93(mV))))^(1/3)
  taum=tauscalem*(0.3632(ms)+1.158(ms)*lgc(-(vscaledm-(-55.96(mV)))/(20.12(mV))))*(3.8/Q)
  vscaledh=v/vscaleh-vshifth
  hinf=(lgc(-(vscaledh-(-53.3(mV)))/(14.54(mV))))^4
  tauh=tauscaleh*(1.24(ms)+2.678(ms)*lgc(-(vscaledh-(-50.0(mV)))/(16.027(mV))))*(3.8/Q)
  }

FUNCTION lgc(x) 
  {
  lgc = 1/(1+exp(-x))
  }
