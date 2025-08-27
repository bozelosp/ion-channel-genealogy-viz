NEURON
 {
	 
  SUFFIX kaDasoma  USEION k READ ek WRITE ik 
 
RANGE gbar, g, i  

RANGE atau, btau

RANGE ainf, binf

RANGE taurecov


 }

UNITS {

(pS) =(picosiemens)
(mV) = (millivolt)
(mA) = (milliamp)

}

PARAMETER { 
  gbar = 150 (pS/microm2)
  Vmid_ac =-30 (mV)
  k_ac = 7 (mV) 
  Vmid_ina = -75 (mV) 
  k_ina = -7 (mV)
   taurecov=25 (ms)
  Vshift=-90 (mV) 
   m=1         
  h=1        

}

 ASSIGNED {
  v	(mV)
  ik 	(mA/cm2)
  i 	(mA/cm2)
  g	(pS/microm2)
  atau (ms)
  btau (ms)
   
  ainf (1)
  binf (1)
 
ek (mV)
 }


STATE {a b b2}

BREAKPOINT {
  SOLVE states METHOD cnexp

 
  g = gbar*(a^m)*(b^h)
  i = (0.0001)*g*(v-ek)
  ik = i
}

INITIAL {

rates(v)
a= ainf
b=binf


}

DERIVATIVE states {
 rates(v)
  a' = (ainf-a)/atau
  b' = (binf-b)/btau
  
}


FUNCTION a_inf (V (mV)) () {

  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  
}

FUNCTION b_inf (V (mV)) () {
  b_inf = 1/(1+exp(-(V-Vmid_ina)/k_ina)) 
}


FUNCTION a_tau (V (mV)) (ms) {
UNITSOFF

a_tau= 1.029 + (4.83/(1+exp((V+57)/6.22)))


UNITSON
}

FUNCTION b_tau (V (mV)) (ms) {

UNITSOFF

  b_tau = taurecov + (80+(78.4/(1+exp(V+68.5)/5.95))-taurecov)/(1+exp((-V+Vshift)*5))
 


UNITSON 

}



PROCEDURE rates(V (mV)) {





atau=a_tau(V)

ainf=a_inf(V)

btau=b_tau(V)

binf=b_inf(V)



}