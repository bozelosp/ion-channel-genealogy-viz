INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica,cai WRITE cai
	RANGE depth,kt,kd,cainf,taur,icaadjust,camolflux,cabuff,capump 

}

UNITS {
	(molar) = (1/liter)      
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)   = (ms mM)
}

CONSTANT{
	FARADAY = 96489 (coul)   
}

PARAMETER {
	depth = .1	(um)     
	taur  = 1e10    (ms)     
	cainf = 2.4e-4	(mM)
	kt    = 1e-4	(mM/ms)
	kd    = 1e-4	(mM)
	icaadjust = 1     
	cainit = 1e-4      (mM)  
	capump = 1
	cabuff = 1
	camolflux = 1 
}

STATE {
	cai   (mM)
}

INITIAL {
	cai = cainit
}

ASSIGNED{
	ica		(mA/cm2)
	drive_channel   (mM/ms)
	drive_pump	(mM/ms)
}

BREAKPOINT{
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {

	drive_channel = -(10000)*(ica*icaadjust)/(2*FARADAY*depth)

	if(drive_channel <= 0.) {drive_channel = 0.}

	drive_pump = -kt*cai/(cai+kd)  

	cai' = (camolflux * drive_channel) + (capump * drive_pump) + (cabuff * ((cainf-cai)/taur))
}