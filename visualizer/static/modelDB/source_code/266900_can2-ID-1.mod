TITLE n-calcium channel
: n-type calcium channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER { 
   
   curr
	v (mV)
     sh=0 (mV)
	celsius 		(degC)
	gcanbar=.0003 (mho/cm2)
	ki=.001 (mM)
	cai=50.e-6 (mM)
	cao = 2  (mM)
	q10=5
	mmin = 0.2
	hmin = 3
	a0m =0.03
	zetam = 2
	vhalfm = -14
	gmm=0.1	
	time1=600
	time0=100
	alphash0=0
	alphash1=0.15
	sh2
	count=1
	vrun (mV)
	delta=0
	vinit=-76.2 
	:vvrun=0
	 timestep=1000
	vrun2
	v0
	dv0
	ddv
	flag=0
	FCa = 2
	PCa = 1
	BCa = 2
	CCa = 50
	stim_moltCa=1
	eca2=140
}


NEURON {
	SUFFIX can
	POINTER stim_i
	USEION ca READ cai,cao WRITE ica
        RANGE flag, curr, gcanbar, ica, gcan,sh , count,delta2,vrun2, stim_moltCa , eca2, vvrun
        GLOBAL  hinf,minf,taum,tauh
}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcan  (mho/cm2) 
        minf
        hinf
        taum
        tauh
		stim_i
		vvrun
}

INITIAL {
        rates(v,sh2)
        m = minf
        h = hinf
		vrun=0
	  vvrun=vrun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcan = gcanbar*m*m*h*h2(cai)
	ica = gcan*ghk(v,cai,cao)

}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
	alph = 1.6e-4*exp(-v/48.4)
}

FUNCTION beth(v(mV)) {
	beth = 1/(exp((-v+39.0)/10.)+1.)
}

FUNCTION alpm(v(mV)) {
	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
	betm = 0.046*exp(-v/20.73)
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm)) 
}

FUNCTION betmt(v(mV)) {
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm)) 
}

UNITSON

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v,sh2)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}


BEFORE STEP { LOCAL i
       
	  if(stim_i==0 && flag==0){ 
		  vrun=0
		  vvrun=0
		  
	    }else{
		 flag=1
		             		  
		delta=v-vinit
		if (count<timestep+1){
		   vrun= (delta-vrun)*(FCa/(count+1))+vrun
	       vrun2=vrun 
		 }else{

		vrun2= (delta)*(FCa/(timestep+1))+vrun2*pow((1-FCa/(timestep+1)),PCa)
			
			}
		  	   
	   vvrun=(BCa*vrun2/(1+vrun2/CCa))
	    
		count=count+1   
        }
								
		 sh2=sh+alphash1*vvrun
	
}
   
PROCEDURE rates(v (mV),sh2) { :callable from hoc
        LOCAL a, b, qt,i
		
			
        a = alpm(v-sh2)
        b = 1/(a + betm(v-sh2))
        minf = a*b
	taum = betmt(v-sh2)/(qt*a0m*(1+alpmt(v-sh2)))
	if (taum<mmin/qt) {taum=mmin/qt}
        a = alph(v-sh2)
        b = 1/(a + beth(v-sh2))
        hinf = a*b
	tauh=b/qt
	:tauh= 80
	if (tauh<hmin) {tauh=hmin}
}
