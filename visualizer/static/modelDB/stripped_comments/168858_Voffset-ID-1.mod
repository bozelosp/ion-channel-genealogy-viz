NEURON {
  SUFFIX offst
  RANGE Voff, Vshift, Vsum
}

PARAMETER { 
	Vraise = -58.5 	(mV)	
	on = 215 	(ms)	
	W = 50		(ms)	
	we = 10		(ms)	
	didAvg = 0	(1)	
}

UNITS { 
	(mV) = (millivolt) 
}

ASSIGNED {
  v (millivolt)
  
	Vsum		(mV)
	npts		(1)
	Voff		(mV)
	Vshift 		(mV)
}

INITIAL {
	Vsum = 0
	Vshift = 0
	npts = 0
	Voff = v - Voff + Vraise
        didAvg = 0

}

BREAKPOINT {
    if( t >= on-we-W && t <= on-we ) {
        Vsum = Vsum + v
        npts = npts + 1

    } 
    if( didAvg == 0 && t > on-we ) {		
        didAvg = 1
        Voff = Vsum / npts
printf("***TIME %g\tFound offset %g = %g / %g\n",t,Voff,Vsum,npts)
    }
    Vshift = v - Voff + Vraise



}