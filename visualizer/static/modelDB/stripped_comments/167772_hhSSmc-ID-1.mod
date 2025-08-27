UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
	SUFFIX hhSSmc
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, NNa, NK, se
  RANGE nextRNa, nextRK, prev_evNa, prev_evK, next_evNa, next_evK,NaT, K_3_4, K_4_3
  RANGE sumrtNa, sumrtK, mh0, n0, N
  RANGE an, bn, am, bm
}
 
PARAMETER {
	se = -1
	gnabar = .12 (S/cm2)	<0,1e9>
	gkbar = .036 (S/cm2)	<0,1e9>
	gl = .0003 (S/cm2)	<0,1e9>
	el = -54.3 (mV)
	NNa = 5000
	NK = 1500 
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)
	dt (ms)
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	am	(/ms)
	ah	(/ms)
	an	(/ms)
	bm	(/ms)
	bh	(/ms)
	bn	(/ms)
	stsum
	M
	N
	H
	NaT[4]	(/ms)
	K_4_3	(/ms)
	K_3_4	(/ms)
	cumsumNa[4]	(/ms)
	prev_evNa	(ms)
	next_evNa	(ms)
	prev_evK	(ms)
	next_evK	(ms)
	nextRNa
	sumrtNa  	(/ms)
	nextRK
	sumrtK		(/ms)
	ev
	mh0    
            
  n0
	
}
 
STATE {	
	mh1
	mh2
	mh3
	mh4
	mh5
	mh6
	mh7
	n1
	n2
	n3
	n4
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*mh7*(v - ena)/NNa
	ik = gkbar*n4*(v - ek)/NK
	il = gl*(v - el)
}
 
INITIAL {
	if (se>0) {set_seed(se)}
	ratesNa(v)
	ratesK(v)
	M=am/bm
	H=ah/bh
	N=an/bn
    
    
	stsum=(1+H)*(1+M)^3
	mh0=NNa/stsum
	mh1=NNa*3*M/stsum
	mh2=NNa*3*M^2/stsum
	mh3=NNa*M^3/stsum
	mh4=NNa*H/stsum
	mh5=NNa*H*3*M/stsum
	mh6=NNa*H*3*M^2/stsum
	mh7=NNa*H*M^3/stsum
	
	stsum=(1+N)^4
	n0=NK/stsum
	n1=NK*4*N/stsum
	n2=NK*6*N^2/stsum
	n3=NK*4*N^3/stsum
	n4=NK*N^4/stsum
	ratesNa(v)
	ratesK(v)
    
    
    
  nextRNa=log(scop_random())
  nextRK=log(scop_random())
	prev_evK=0
	prev_evNa=0
  
}

DERIVATIVE states {  
	ratesNa(v)
	
	mh1' = (-2*am-bm-ah)*mh1 + 3*am*mh0 + 2*bm*mh2 + bh*mh5 
	mh2' = (-am-2*bm-ah)*mh2 + 2*am*mh1 + 3*bm*mh3 + bh*mh6 
	mh3' = (-3*bm)*mh3 + am*mh2  			
	mh4' = (-3*am-bh)*mh4 + bm*mh5 + ah*mh0 
    mh5' = (-2*am-bm-bh)*mh5 + 3*am*mh4 + 2*bm*mh6 + ah*mh1 
   	mh6' = (-2*bm-bh)*mh6 + 2*am*mh5 + ah*mh2    	
    
    
    
	next_evNa = prev_evNa - nextRNa/sumrtNa  
	while (t>= next_evNa){
		transNa()                 
		prev_evNa = next_evNa     
		next_evNa = prev_evNa - nextRNa/sumrtNa  
	}
    
	mh0 = NNa-mh1-mh2-mh3-mh4-mh5-mh6-mh7
	
	ratesK(v)		
	
	n1' = (-3*an-bn)*n1 + 4*an*n0 + 2*bn*n2 
	n2' = (-2*an-2*bn)*n2 + 3*an*n1 + 3*bn*n3 
	n3' = (-3*bn)*n3 + 2*an*n2
    
	next_evK = prev_evK - nextRK/sumrtK
	while (t>= next_evK){
		transK()
        prev_evK = next_evK
		next_evK = prev_evK - nextRK/sumrtK
	}
	n0 = NK-n1-n2-n3-n4
}
 
LOCAL q10

PROCEDURE ratesNa(v(mV)) {  

	UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
	am = q10*0.1*(v+40)/(1-exp(-(v+40)/10))
	bm = q10*4*exp(-(v+65)/18)
	ah = q10*0.07*exp(-(v+65)/20) 
	bh = q10/(1+exp(-(v+35)/10))

    
	NaT[0]=ah*mh3  
	NaT[1]=bh*mh7	
	NaT[2]=am*mh6	
	NaT[3]=3*bm*mh7	
	sumrtNa=NaT[0] + NaT[1] + NaT[2] + NaT[3]
}
		
PROCEDURE ratesK(v(mV)) {  
	q10 = 3^((celsius - 6.3)/10)
	an = q10*0.01*(v+55)/(1-exp(-(v+55)/10))
	bn = q10*0.125*exp(-(v+65)/80)
	K_3_4 = an*n3     
	K_4_3 = 4*bn*n4   
	sumrtK = K_3_4 + K_4_3

	UNITSON 
}


PROCEDURE transNa() {
	ratesNa(v)
	UNITSOFF
	sumrtNa=0
    
	FROM ii=0 TO 3 {
		sumrtNa = sumrtNa + NaT[ii]
		cumsumNa[ii] = sumrtNa
	}
    
	FROM ii=0 TO 3 {cumsumNa[ii] = cumsumNa[ii] / sumrtNa}

    
	ev = scop_random()
	if (ev <= cumsumNa[0]) { 
		mh3=mh3-1
		mh7=mh7+1
	}
	if (cumsumNa[0] < ev && ev <= cumsumNa[1]) {  
		mh3=mh3+1
		mh7=mh7-1
	}	
	if (cumsumNa[1] < ev && ev <= cumsumNa[2]) {  
		mh6=mh6-1
		mh7=mh7+1
	}
	if (cumsumNa[2] < ev && ev <= cumsumNa[3]) {  
		mh6=mh6+1
		mh7=mh7-1
	}
	UNITSON
    
    
	if (mh3>NNa){mh3=NNa}
	if (mh6>NNa){mh6=NNa}
	if (mh7>NNa){mh7=NNa}
	if (mh3<0){mh3=0}
	if (mh6<0){mh6=0}
	if (mh7<0){mh7=0}

	nextRNa=log(scop_random())
}

PROCEDURE transK() {
	ratesK(v)
   	sumrtK=0
	sumrtK = K_3_4 + K_4_3
	ev = scop_random()
    
    
	if (ev <= (K_3_4 / sumrtK)) {
		n3=n3-1
		n4=n4+1
	} else {
		n3=n3+1
		n4=n4-1
	}
    
	if (n3>NK){n3=NK}
	if (n4>NK){n4=NK}
	if (n3<0){n3=0}
	if (n4<0){n4=0}
	nextRK = log(scop_random())
}