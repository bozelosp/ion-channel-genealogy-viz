NEURON {
	SUFFIX cadifus2
	USEION ca READ cao, cai, ica WRITE cai, ica
	GLOBAL vol, Buffer0
	RANGE ipump
}
DEFINE NANN  4

UNITS {
        (mol)   = (1)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	FARADAY = (faraday)	 (10000 coulomb)
	PI	= (pi) (1)
}

PARAMETER {
	DFree = .6	(um2/ms)
	diam		(um)
	cao		(mM)
	ica		(mA/cm2)
	k1buf = 50	(/mM-ms)
	k2buf = 5	(/ms)
        k1=1.e10            (um3/s)
        k2=50.e7            (/s)	
        k3=1.e10            (/s)	
        k4=5.e6	            (um3/s)	
	area		(um2)
} 
CONSTANT { volo=1  (liter)}

ASSIGNED {
	cai		(mM)
	vol[NANN]	(1)	
	ipump           (mA/cm2)
	last_ipump           (mA/cm2)

}

STATE {
	ca[NANN]	(mM) <1.e-5> 
	CaBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
        pump            (mol/cm2) <1.e-15>
        pumpca          (mol/cm2) <1.e-15>

}

LOCAL totpump, kd,totbuf

INITIAL {
           totpump=0.2
           pump=totpump/(1+1.e-18*k4*cao/k3)
           pumpca=2.e-22
	   ipump=0

           totbuf=1.2
           kd=k2buf/k1buf
           FROM i=0 TO NANN-1 {
                ca[i] = cai
		CaBuffer[i] =(totbuf*ca[i])/(kd+ca[i])
		Buffer[i] = totbuf - CaBuffer[i]
                }

}

BREAKPOINT {
	SOLVE state METHOD sparse
	last_ipump=ipump
	ica = ipump
}

LOCAL coord_done

INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	
	
	FROM i=0 TO NANN-1 {
		ca[i] = cai
	}
}

LOCAL frat[NANN] 	

PROCEDURE coord() {
	LOCAL r, dr2
	
	
	
	
	
	
	r = 1/2					
	dr2 = r/(NANN-1)/2			
	vol[0] = 0
	frat[0] = 2*r
	FROM i=0 TO NANN-2 {
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2	
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	
					
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2	
	}
}

LOCAL dsq, dsqvol 
		
KINETIC state {
	COMPARTMENT i, diam*diam*vol[i]*1(um) {ca CaBuffer Buffer}
        COMPARTMENT (1.e10)*area {pump pumpca}
        COMPARTMENT (1.e15)*volo {cao}

	~ ca[0] << (-(ica-last_ipump)*PI*diam*frat[0]*1(um)/(2*FARADAY))
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DFree*frat[i+1]*1(um), DFree*frat[i+1]*1(um))
	}
	dsq = diam*diam*1(um)
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol,k2buf*dsqvol)
	}
        ~ca[0] + pump <-> pumpca ((1.e-11)*k1*area, (1.e7)*k2*area)
        ~pumpca       <-> pump + cao ((1.e7)*k3*area, (1.e-11)*k4*area)

        ipump = 2*FARADAY*(f_flux-b_flux)/area

	cai = ca[0]
}