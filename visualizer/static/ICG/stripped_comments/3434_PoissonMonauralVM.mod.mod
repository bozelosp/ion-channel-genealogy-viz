NEURON {

 
   POINT_PROCESS PoissonMonauralVM
 

 RANGE 

 
  
   vs, freq, stimPhaseIpsi, stimPhaseContra, stimRate,
   genericParam1, genericParam2,
  
   angfreq, stimProb, kappa, norm, stim,
 
 
  poissDt, tickSize, tickVal, tickOnce 
 POINTER probIpsi, probContra
}

UNITS {

 
  (kHz) = (kiloHz)
 

 PI = (pi) (1)
 (tick) = (1)
 (prob) = (1)
 (angle) = (1) 
}

PARAMETER {

 
  vs = 0.43 (1)               <0,0.999> 
  freq = 1000 (Hz)
  stimPhaseIpsi = 0 (angle)
  stimPhaseContra = 0 (angle)
  stimRate = 3.46 (/ms)       <1e-9,1e9>
  genericParam1 = 0 (1)       <0,1e9> 
  genericParam2 = 0     
 

 poissDt = 0.025 (ms)        <1e-9,1e9>
 tickSize = 1 (tick)         <1e-9,1e9>
}

ASSIGNED {
 
  angfreq (/ms)
  stimProb (prob)
  kappa (1)
  norm (1)
 

 tickVal (tick)
 tickOnce (tick)
 stim[2]  (prob) 
 probIpsi (prob)
 probContra (prob)

}

INITIAL {

 
  angfreq = 2*PI*freq*(0.001) 
  stimProb = stimRate*poissDt
  if (vs < 0) {vs = 0} 
  if (vs >= 0.999) {vs = 0.999} 
  kappa =  A1Inv(vs)
  norm = 1/(2*PI*IO(kappa))
 

 tickVal = 0
 tickOnce = 0

 stim[0] = calcStimIpsi()
 stim[1] = calcStimContra()
 if (pointersValid()) {
  probIpsi = stimProb*stim[0]
  probContra = stimProb*stim[1]  
 }
}


 
 FUNCTION calcStimIpsi() {

  

  calcStimIpsi = norm*exp(kappa*sin(angfreq*t+stimPhaseIpsi))

 }
 
 FUNCTION calcStimContra() {

  

  calcStimContra = genericParam1

 }
 



 FUNCTION IO(x) { 
 
 
  LOCAL w, w2, w4, w6, w8 
  w = x/3.75
  w2 = w*w
  w4 = w2*w2
  w6 = w4*w2
  w8 = w4*w4
  if (w < 1) {
   IO = (1 + 3.5156229*w2 + 3.0899424*w4 + 1.2067492*w6 + 0.2659732*w8 + 
             0.0360768*w8*w2 + 0.0045813*w8*w4)
  } else {
   IO = (0.39894228 + 0.01328592/w + 0.00225319/w2 - 0.00157565/(w2*w) + 
         0.00916281/w4 - 0.02057706/(w4*w) + 0.02635537/w6 - 
         0.01647633/(w6*w) + 0.00392377/w8)*exp(x)/sqrt(x)
  }
 }
 
 FUNCTION A1Inv(x) { 
 
 
  if (x < 0.53) {
   A1Inv = (2*x + x*x*x + 5/6*x*x*x*x*x)
  } else { 
   if (x < 0.85) {
    A1Inv = (-0.4 + 1.39*x + 0.43/(1-x))
   } else {
    A1Inv = (1/(x*x*x - 4*x*x + 3*x))
   }
  }
 }





PROCEDURE prerun() {
 tickVal  = tickSize
 tickOnce = tickSize
}

NET_RECEIVE (weight) {
 if(tickVal) {
  tickVal = 0
 } else {
  tickVal = tickSize
  stim[0] = calcStimIpsi()
  stim[1] = calcStimContra()
  probIpsi   = stimProb*stim[0]
  probContra = stimProb*stim[1]
 }
 net_send(poissDt,1)
}

FUNCTION pointersValid() {
 if ((!nrn_pointing(probIpsi))||(!nrn_pointing(probContra))) {
  printf(
   "Poisson Process
  pointersValid = 0
 } else {
  pointersValid = 1
 }
}