TITLE hh_mod.mod   squid sodium, potassium, and leak channels
 
 
: Modified from NEURON's built-in hh.mod 
: 2021, Janos Brunner
:
: relates to the publication:
: "Small size of recorded neuronal structures confines the accuracy in direct axonal voltage measurements" 
: by Viktor Janos Olah, Gergely Tarcsay and Janos Brunner 
: ENEURO.0059-21.2021


: scale_na_act, scale_na_inact, scale_k, v_shift_na, v_shift_k are introduced
 
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT


	



 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hh_mod
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, scale_na_act,scale_na_inact,scale_k, v_shift_na,v_shift_k
        RANGE minf, hinf, ninf, mtau, htau, ntau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        el = -80 (mV)
		scale_na_act=1	<0.001,1e3>
		scale_na_inact=1	<0.001,1e3>
		scale_k=1	<0.001,1e3>
		v_shift_na=0
		v_shift_k=0
		
		
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 



? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       

UNITSOFF
        alpha = .1 * vtrap(-((v-v_shift_na)+40),10)
        beta =  4 * exp(-((v-v_shift_na)+65)/18)
        sum = alpha + beta
	mtau = 1/(scale_na_act*sum)
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .07 * exp(-((v-v_shift_na)+65)/20)
        beta = 1 / (exp(-((v-v_shift_na)+35)/10) + 1)
        sum = alpha + beta
	htau = 1/(scale_na_inact*sum)
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .01*vtrap(-((v-v_shift_k)+55),10) 
        beta = .125*exp(-((v-v_shift_k)+65)/80)
	sum = alpha + beta
        ntau = 1/(scale_k*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
