NEURON {
	SUFFIX saccum

	USEION na READ ina, nao, nai WRITE nai, nao, ina
	USEION k  READ ik, ko, ki WRITE ki, ko, ik
	USEION ca READ ica, cai WRITE cai, cao, ica VALENCE 2
	USEION cl READ icl, cli, clo WRITE cli, clo, icl VALENCE -1
	USEION hco3  READ hco3o, hco3i VALENCE -1
	USEION a  READ ia, ao, ai WRITE ai, ao VALENCE -1
	
	RANGE volin, volout, delta
	POINTER ko2, ko3, nao2, nao3
	GLOBAL TotalBuffer
}

UNITS	{
	(mV)		= (millivolt)
	(mA)		= (milliamp)
	FARADAY		= 96485.309 (coul/mole)
	(molar)		= (1/liter)
	(mM)		= (millimolar)
	PI		    = (pi) (1)
	R 	    	= (k-mole) (joule/degC)
}

PARAMETER {
	
	Difna = 1.33	(um2/ms)   
	Difk  = 1.96	(um2/ms)   
	Difca = 0.6	    (um2/ms)   
	Difcl = 2.03	(um2/ms)   
	Difa = 0	    (um2/ms)

	
	k1buf = 0.0008 (/ms)       
	TotalBuffer = 1100 (mM)    
	x0 = 16 (mM)               
    q = -1.25 (mM)             
    bf = 1
    
	
	nab = 140 (mM) 
	kb  = 3.5 (mM)
	clb = 135 (mM)
	
	s = 4.4e4
	sk = 4.4e4
    bath = 1
	nadif = 1
    rad = 1
    narad = 1

	
	setvolin = 1	
	setvolout = 1   

	minvol = 0.04      
	maxvol = 0.9       
	tau = 250 (ms)     
	swell = 1          
}

ASSIGNED {
	
	ina		(mA/cm2)
	ik		(mA/cm2)
	ica		(mA/cm2)
	icl		(mA/cm2)
	ia		(mA/cm2)

	
	diam		(um)

	
	naflux[2]
	kflux[2]
	caflux[2]
	clflux[2]
	aflux[2]

	
	dif[3]

	
	shell[4]
	
	nai ki cai cli hco3i ai		
	nao ko cao clo hco3o ao		

	volin     	
	volout		

	
	Kd (/mM)
	B0 (mM)

	
	delta (mM)

	
	ko2 (mM)
	ko3 (mM)
	nao2 (mM)
	nao3 (mM)
}

STATE { na[2] k[2] ca[2] cl[2] a[2] kbuf Buffer KBuffer (mM) vol[2] <1e-4> }

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL {
	na[0]=nai
	na[1]=nao
	k[0]=ki
	k[1]=ko
	ca[0]=cai
	ca[1]=cao
	cl[0]=cli
	cl[1]=clo
	a[0]=ai
	a[1]=ao
	
	
	Kd = k2buf(k[1])/k1buf
    kbuf = 0
	B0 = TotalBuffer/(1+Kd*k[1])
	Buffer = B0
	KBuffer = TotalBuffer - B0
	
	vol[0]=setvolin
	vol[1]=setvolout

	volin=vol[0]
	volout=vol[1]
}

KINETIC state {

	

	delta = ((nai + ki + cli + cai + hco3i + ai) - (nao + ko + clo + cao + hco3o + ao))/tau	

	IF (vol[1] <= minvol && delta > 0 ) {
	  delta = 0
	} ELSE {
	  IF (vol[0] <= maxvol && delta < 0 ) {
	    delta = 0
	  }
    ~  vol[0] << (swell*delta/(diam*diam*PI/4)) 
	~  vol[1] << (-swell*delta/(diam*diam*PI/4)) 
	}

	COMPARTMENT i, vol[i]*PI/4*diam*diam { na k ca cl a }

	
	naflux[0] = -delta*na[0] -(ina*diam)*PI*(1e4)/FARADAY
    kflux[0]  = -delta*k[0] -(ik*diam)*PI*(1e4)/FARADAY
	caflux[0] = -delta*ca[0] -(ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[0] = -delta*cl[0] -(icl*diam)*PI*(1e4)/(FARADAY*-1)
    aflux[0]  = -delta*a[0]
    	
	
	naflux[1] = delta*na[1] + (ina*diam)*PI*(1e4)/FARADAY
    kflux[1]  = delta*k[1]  + (ik*diam)*PI*(1e4)/FARADAY
	caflux[1] = delta*ca[1] + (ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[1] = delta*cl[1] + (icl*diam)*PI*(1e4)/(FARADAY*-1)
	aflux[1]  = delta*a[1]

	
	dif[0]    = nadif*bath*Difna*((geom(diam)/2)/2)*((nab-na[1])/s)*(PI*(diam+dr(diam)))
	dif[1]    = bath*Difk*((geom(diam)/2)/2)*((kb-k[1]-kbuf)/sk)*(PI*(diam+dr(diam)))
	dif[2]    = bath*Difcl*((geom(diam)/2)/2)*((clb-cl[1])/s)*(PI*(diam+dr(diam)))

	
	shell[0]  = rad*Difk*(geom(diam)/2)*surf(diam)*(ko2-k[1]-kbuf)
	shell[1]  = rad*Difk*(geom(diam)/2)*surf(diam)*(ko3-k[1]-kbuf)
	shell[2]  = narad*Difna*(geom(diam)/2)*surf(diam)*(nao2-na[1])
	shell[3]  = narad*Difna*(geom(diam)/2)*surf(diam)*(nao3-na[1])
	
	
	IF (bf == 0 ) {
	  kbuf = 0
      Buffer = 0
      KBuffer = 0
	} ELSE {
	~  kbuf  << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  Buffer << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  KBuffer << (k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer-k1buf*KBuffer)
        }

	
	~  na[0] << (naflux[0])
    ~  k[0]  << (kflux[0])
	~  ca[0] << (caflux[0])
	~  cl[0] << (clflux[0])
	~  a[0] << (aflux[0])
	~  na[1] << (naflux[1]+dif[0]+shell[2]+shell[3])
	~  k[1]  << (kflux[1]+dif[1]+shell[0]+shell[1])
	~  ca[1] << (caflux[1])
	~  cl[1] << (clflux[1]+dif[2])
	~  a[1] << (aflux[1])

	
	volin=vol[0]
	volout=vol[1]
	nai=na[0] nao=na[1]
	ki=k[0]   ko=k[1]+kbuf
	cai=ca[0] cao=ca[1]
	cli=cl[0] clo=cl[1]
	ai=a[0]   ao=a[1]
}



FUNCTION k2buf(x(mM)) { 	         
	TABLE FROM 0 TO 199 WITH 200
	k2buf = k1buf/(1+exp((x-x0)/(q)))
}

FUNCTION dr(x(um)) { 			 
	TABLE FROM 0 TO 15 WITH 10
	dr = x*(sqrt(1+setvolout)-1)
}

FUNCTION geom(x(um)) {
	TABLE FROM 0 TO 15 WITH 10
	geom = 1/dr(x)
}

FUNCTION surf(x(um)) { 			 
	TABLE FROM 0 TO 15 WITH 10
	surf = (PI-2*atan2(x+2*dr(x/2),(x/2)))*((x+2*dr(x/2))/sin(atan2(x+2*dr(x/2),(x/2))))
}