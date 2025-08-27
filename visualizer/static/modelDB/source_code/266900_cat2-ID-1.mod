TITLE T-calcium channel
: T-type calcium channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER { 
   
    curr
	v (mV)
     sh=0 (mV) 
	celsius = 34	(degC)
	gcatbar=.003 (mho/cm2)
	cai = 50.e-6 (mM)
	cao = 2 (mM)
	q10 = 5
	mmin=0.2
	hmin=10
	a0h =0.015
	zetah = 3.5
	vhalfh = -75
	gmh=0.6	
	a0m =0.04  :0.04
	zetam = 2
	vhalfm = -28
	gmm=0.1	
     ki=.001 (mM)
	time1=600
	time0=100
	alphash0=0
	alphash1=0.15
	sh2
	count=1
	vrun (mV)
	delta=0
	vinit=-76.2
	vvrun=0
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
}


NEURON {
	SUFFIX cat
	POINTER stim_i
	USEION ca READ cai,cao WRITE ica
        RANGE flag, curr, gcatbar, ica, gcat,sh,count,delta2,vrun2, stim_moltCa
        GLOBAL  hinf,minf,mtau,htau
}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcat (mho/cm2)
	hinf
	htau
	minf
	mtau
	stim_i
	
}

INITIAL {
	rates(v,sh2)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gcatbar*m*m*h
	ica = gcat*ghk(v,cai,cao)

}

DERIVATIVE states {	: exact when v held constant
	rates(v,sh2)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
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

FUNCTION KTF(celsius (DegC)) (mV) {
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
  alph = exp(0.0378*zetah*(v-vhalfh)) 
}

FUNCTION beth(v(mV)) {
  beth = exp(0.0378*zetah*gmh*(v-vhalfh)) 
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm)) 
}

FUNCTION betmt(v(mV)) {
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm)) 
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
	LOCAL a,b, qt,i
        qt=q10^((celsius-25)/10)
						

	a = 0.2*(-1.0*v+19.26+sh2)/(exp((-1.0*v+19.26+sh2)/10.0)-1.0)
	b = 0.009*exp(-v+sh2/22.03)
	minf = a/(a+b)
	mtau = betmt(v-sh2)/(qt*a0m*(1+alpmt(v-sh2)))
	if (mtau<mmin) {mtau=mmin}

	a = 1.e-6*exp(-v+sh2/16.26)
	b = 1/(exp((-v+29.79+sh2)/10.)+1.)
	hinf = a/(a+b)
	htau = beth(v-sh2)/(a0h*(1+alph(v-sh2)))
	if (htau<hmin) {htau=hmin}
}

