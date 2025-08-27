NEURON {

POINT_PROCESS Dipole
RANGE ri, ia, Q, ztan
POINTER pv


RANGE Qsum 
}

UNITS {
	(nA) = (nanoamp)
	(mV) =(millivolt)
	(Mohm) = (megaohm)
	(um) = (micrometer)
	(Am) = (amp meter)
	(fAm) = (femto amp meter)
}

ASSIGNED {
	ia (nA)
	ri (Mohm)
	pv (mV)
	v (mV)
	ztan (um)
	Q  (fAm)
	Qsum (fAm)
}

AFTER SOLVE {     	
	ia=(pv-v)/ri
	Q=ia*ztan
	Qsum = Qsum + Q
}
	
AFTER INITIAL {
	ia=(pv-v)/ri
	Q=ia*ztan
	Qsum = Qsum + Q
}


 BEFORE INITIAL {
	Qsum = 0
 }
 BEFORE BREAKPOINT {
	Qsum = 0
 }