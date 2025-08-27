NEURON{
	SUFFIX cadyn
	USEION ca READ ica, cai, cao WRITE cai, ica VALENCE 2 		
	USEION na READ nai VALENCE 1								
	USEION caer READ caeri WRITE caeri VALENCE 2 				
	USEION camt READ camti WRITE camti VALENCE 2 				
	USEION caip3r READ caip3ri WRITE caip3ri VALENCE 2			
	USEION ip3 READ ip3i VALENCE 1								
   
	GLOBAL DCa, cai0, caeri0, camti0							
	RANGE k1, k2, k3, k4, ica_pmp, pump0 			 		 	
	GLOBAL bbr													
	
	
	RANGE vmaxsr, kpsr		                        			
	RANGE kactip3, konip3, kinhip3,  kip3, jmaxsr			  	
	RANGE ktcicr, kcicr, vcicr                           		
    RANGE vmcu, kmcu, nmcu, vncx, kna, kncx			    	   	
	RANGE jer, jserca, jip3, jcicr								
	RANGE jmcu, jmncx, jmito

	RANGE Kmmt, Bmmt, fmmt										
	RANGE Kmer, Bmer, fmer										
	
    GLOBAL vol
	THREADSAFE
}
DEFINE NANN  12    

UNITS{
	(mV)	= (millivolt)
	(um)    = (micron)
	(mM)    = (milli/liter)
    (nM)    = (nano/liter)
	(mA)    = (milliamp)
	F       = (faraday) (coulombs)
	PI      = (pi) (1)
	R 		= (k-mole) (joule/degC)
    (mol)   = (1) 
}

PARAMETER{
	diam	(um)
	L		(um)

	
	cai0 	= 136e-6	(mM)
	cao0 	= 2 		(mM)
	caeri0 	= 0.4		(mM)
   	camti0 	= 2e-4 		(mM)
	DCa 	= 0.6		(um2/ms)     

	
	bbr = 370 		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	vmaxsr 	= 0.00027	(mM/ms)
	kpsr 	= 3.75e-6 	(mM) 
	

	
    jmaxsr 	= 3.5e-6	(mM/ms)
    kip3 	= 0.0008    (mM) 	
	kactip3 = 0.0003	(mM) 	
	konip3 	= 2.7		(/mM-ms)
	kinhip3 = 0.0002	(mM)    


	
	
	kcicr 	= 0.00198	(mM) 
	ktcicr 	= 0.0006  	(mM)
	vcicr 	= 5e-7     	(/ms)

	
	Kmer = 0.5 (mM)
	Bmer = 10 (mM)
	


	
	
	vmcu = 1.4468e-6	(mM/ms)  
	kmcu = 606e-6 	(mM) 	     
    nmcu = 2.3 (1) 				 
	
	
	vncx = 6e-5 	(mM/ms)
	kna = 	8		(mM)	
	kncx = 	35e-3	(mM)	


	
	Kmmt = 0.01e-3 (mM) 
	Bmmt = 0.065 (mM)
	

	
	k1 = 3.74e7        (/mM-s) 
	k2 = .25e6      (/s)
	k3 = .5e3       (/s)
	k4 = 5e0        (/mM-s)
	pump0 = 1.3725e-13 (mol/cm2)  
}

ASSIGNED{
	celsius		(degC)
	ica			(mA/cm2)
	
	cai			(mM)
	cao			(mM)
	caeri		(mM)
	camti 		(mM)
	
	vol[NANN]	(1)

	
	ica_pmp (mA/cm2)
	last_ica_pmp (mA/cm2)
	parea    (um)
	c1      (1+8 um4/ms)
	c2      (1-10 um/ms)
	c3      (1-10 um/ms)
	c4      (1+8 um4/ms)


	
	ip3i		(mM)


	
	nai			(mM)

	
	jmcu[NANN] 			(mM/ms)
	jmncx[NANN]			(mM /ms)
	jmito[NANN]  		(mM /ms)    

	
	jer[NANN]		(mM /ms)
	jip3[NANN]		(mM /ms)
	jserca[NANN]	(mM /ms)
	jcicr[NANN]     (mM /ms) 
	
	
	
	
	
	fmer[NANN]
	fmmt[NANN]
}

CONSTANT{
volo = 1e10 (um2) }


STATE{
	
	pump            (mol/cm2) <1e-16> 
	pumpca          (mol/cm2) <1e-16>

	ca[NANN]			(mM)     <1e-8>
	caer[NANN]			(mM)
	camt[NANN]			(mM)
	caip3ri				(mM)

	hc[NANN]        (1)		
	ho[NANN]		(1)     

	Ln[NANN]		(mM/ms)	
}

LOCAL factors_done

INITIAL{ LOCAL total
	
	if (factors_done==0) {
		factors_done= 1
		factors()
	}

	
	cai = cai0	

	
	parms() 
	parea = PI*diam
	pump = pump0
	pumpca = cai*pump*k1/k2
	total = pumpca + pump
	if (total > 1e-9) {
		pump = pump*(pump/total)
		pumpca = pumpca*(pump/total)
	}
	ica_pmp = 0
	last_ica_pmp = 0



	FROM i=0 TO NANN-1{
		ca[i]=cai0          
		caer[i] = caeri0	
		camt[i] = camti0	
		caip3ri = cai0
		
		
		jserca[i] = 0
		jip3[i] = 0
		jcicr[i] = 0
		
		
		jmcu[i] = 0
		jmncx[i] = 0
		
		
		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)
		
		
		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2)
	}
	caeri= caer[0]
	camti = camt[0]

	
	FROM i=0 TO NANN-1 {
    		 ho[i] = kinhip3/(ca[i]+kinhip3) 	        
    		 hc[i] = 1-ho[i] 	   			            

			   jserca[i] = (-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2))
			   
			   jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )
			   
			   if(ca[i] > ktcicr){			
					jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) ) 
				} else {
					jcicr[i] = 0
				} 			
		
			   
			   jer[i] = jserca[i]+jip3[i]+jcicr[i]
			   
			   UNITSOFF
			   jmcu[i] = (-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))
			   UNITSON
			   jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))
			   
			   jmito[i] = jmcu[i]+jmncx[i]
			   
			   Ln[i] = -(jserca[i]+jip3[i]+jcicr[i])/(1 - (ca[i]/caeri0))
    		 }
}

BREAKPOINT{
	SOLVE state METHOD sparse
	 last_ica_pmp = ica_pmp
     ica = ica_pmp
}

LOCAL frat[NANN]

PROCEDURE factors(){
	LOCAL r, dr2
	r = 1/2		        
	dr2 = r/(NANN-1)/2	
	vol[0] = 0
	frat[0] = 2*r
	FROM i=0 TO NANN-2{
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2 
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)		
                                        
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2   
		}
}

LOCAL dsq, dsqvol,dsqvolmt,dsqvoler

KINETIC state {
	COMPARTMENT ii, (1+bbr)*diam*diam*vol[ii]*0.81 {ca} 
	
	COMPARTMENT     (1+bbr)*diam*diam*vol[0]*0.81 {caip3ri}
	
	
	COMPARTMENT jj,	(1/fmer[jj])*diam*diam*vol[jj]*0.12 {caer} 
	COMPARTMENT kk, (1/fmmt[kk])*diam*diam*vol[kk]*0.07 {camt} 
	COMPARTMENT (1e10)*parea {pump pumpca}  
	COMPARTMENT volo {cao}


	~ ca[0] << (-(ica - last_ica_pmp)*PI*diam*(1e4)*frat[0]/(2*F))
	
	
	
	 ~ ca[0] + pump <-> pumpca  (c1,c2)  
	 ~ pumpca <-> pump + cao    (c3,c4)
	  
	ica_pmp = (1e-4) * 2*F*(f_flux - b_flux)/parea 
	  
   
	
	 FROM i=0 TO NANN-2{
		
		 ~ ca[i] <-> ca[i+1]	(DCa*frat[i+1], DCa*frat[i+1])
}
         
	 dsq = diam*diam
     	 
	 FROM i=0 TO NANN-1{
		 dsqvol = dsq*vol[i]*0.81
		 dsqvoler = dsq*vol[i]*0.12

		
	
		
		jserca[i] = ((-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2)))
		
		~ ca[i] << (dsqvol*jserca[i])
		~ caer[i] << (-dsqvoler*jserca[i])
		
		
		~ hc[i] <-> ho[i]  (konip3*kinhip3, konip3*ca[i])
		jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )
        
		~ ca[i] << (dsqvol*jip3[i])
		~ caer[i] << (-dsqvoler*jip3[i])
		
		
		if (i==0) {
			~ caip3ri << (dsqvol*jip3[0])
		}
		
		
		if(ca[i] > ktcicr){			
			jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) )
			~ ca[i] << (dsqvol * jcicr[i])
			~ caer[i] << (-dsqvoler * jcicr[i])
		} else {
			jcicr[i] = 0
			~ ca[i] << (dsqvol*jcicr[i])
			~ caer[i] << (-dsqvoler*jcicr[i])
		} 			
		
		
		~ ca[i] << (Ln[i]*(1-ca[i]/caeri0)*dsqvol)
		~ caer[i] << (-Ln[i]*(1-ca[i]/caeri0)*dsqvoler)
		
		jer[i] = jserca[i]+jip3[i]+jcicr[i]
		
		
		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2) 

		
		dsqvolmt = dsq*vol[i]*0.07
		
		
		UNITSOFF
		jmcu[i] = ((-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))*1/(camt[i]*1e3))
		UNITSON
		~ ca[i] << (dsqvol*jmcu[i])
		
		
		jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))
		~ ca[i] << (jmncx[i]*dsqvol)
		
		
		jmito[i] = jmcu[i]+jmncx[i]
		
		
		~ camt[i] <<  (-(jmncx[i]+jmcu[i])*dsqvolmt) 
		
		
		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)
		
		
		
	}
	
	cai = ca[0]
	caeri= caer[0]
	camti = camt[0]
}

PROCEDURE parms() {
	parea = 2*PI*(diam/2)
        c1 = (1e7)*parea * k1
        c2 = (1e7)*parea * k2
        c3 = (1e7)*parea * k3
        c4 = (1e7)*parea * k4
}