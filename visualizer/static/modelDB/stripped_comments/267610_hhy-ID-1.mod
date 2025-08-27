UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)

}

NEURON {
        SUFFIX hhy

        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk,a2,s_init,smin,tslow

        GLOBAL minf, hinf, ninf, mtau, htau, ntau

}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
        
        gamma = 0.0000025 (/ms)
        delta = 0.25 (/ms)
        
    	vhalfs=-60	(mV)		
        a0s=0.0003	(ms)		
        smin=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        a2=1		(1)		
        s_init = 1
        tslow = 0
}
 
STATE {
        m h n s
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

		gna (S/cm2)
		gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf sinf
		mtau (ms)
		htau (ms)
		ntau (ms)
		stau (ms)
}


BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h*s
		ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
		ik = gk*(v - ek)      
        il = gl*(v - el)
}


DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		s' = (sinf-s)/stau

}

UNITSOFF

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	s = s_init       
}

PROCEDURE rates(v(mV)) { 

        LOCAL  alpha, beta, sum, q10,cs,avail

        q10 = 3^((celsius - 6.3)/10)
                
        alpha = .1 * vtrap(-(v+40),10)
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
		mtau = 1/(q10*sum)
        minf = alpha/sum
                
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
		htau = 1/(q10*sum)
        hinf = alpha/sum
                
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
		sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
 		
        cs = 1/(1+exp((v-vvh)/vvs))
        avail = 1
        if (t>2000) {avail = a2}
        sinf = cs+avail*(1-cs)
        if (t<2000) {sinf = s_init}
        
        alpha = exp(0.45*(v-vhalfs))
        beta = exp(0.09*(v-vhalfs))
        stau = beta/(a0s*(1+alpha))
        
        if (stau<smin) {stau=smin}
}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON