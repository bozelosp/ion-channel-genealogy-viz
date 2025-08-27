INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	(um) = (micron)
	(mol) = (1)
	PI = (pi) (1)
	FARADAY = (faraday) (coulomb)
}

NEURON {
	SUFFIX capump
	USEION ca READ cao, cai WRITE cai, ica
	GLOBAL k1, k2, k3, k4
}

STATE {
	pump	(mol/cm2)
	pumpca	(mol/cm2)
	cai	(mM)
}

PARAMETER {
	cao = 10	(mM)
	diam = 2	(um)

	k1 = 5e8	(/mM-s)
	k2 = .25e6	(/s)
	k3 = .5e3	(/s)
	k4 = 5e0	(/mM-s)
}

CONSTANT {
	volo = 1 (liter)
}

ASSIGNED {
	ica (mA/cm2)
	ipump (mA/cm2)
	voli	(um3)
	area1	(um2)
	c1	(1+8 um5/ms)
	c2	(1-10 um2/ms)
	c3	(1-10 um2/ms)
	c4	(1+8 um5/ms)
}

BREAKPOINT {
	if (t == 0) {parms()}
	SOLVE pmp METHOD sparse
	ica = ipump
}

KINETIC pmp {

	COMPARTMENT voli {cai}
	COMPARTMENT (1e10)*area1 {pump pumpca}
	COMPARTMENT volo*(1e15) {cao}

	~ cai + pump <-> pumpca		(c1,c2)
	~ pumpca     <-> pump + cao	(c3,c4)

	
	ipump = (1e-4)*2*FARADAY*(f_flux - b_flux)/area1
}

INITIAL {
	
	
	VERBATIM
	cai = _ion_cai;
	ENDVERBATIM
}

PROCEDURE parms() { 
	voli = PI*(diam/2)^2 * 1(um)
	area1 = 2*PI*(diam/2) * 1(um)
	c1 = (1e7)*area1 * k1
	c2 = (1e7)*area1 * k2
	c3 = (1e7)*area1 * k3
	c4 = (1e7)*area1 * k4
}

FUNCTION ss() (mM) {	
	SOLVE pmp STEADYSTATE sparse
	ss = cai
	
}