NEURON {
	SUFFIX monx


	RANGE vmax, dist, surf, temp, zin, gesyn, gisyn
}

ASSIGNED {





	area (micron2)
	surf (micron2)
	v (millivolt)
	vmax (millivolt)
	tvmax (ms)
	dist (micron)
	gesyn (microsiemens)
	gisyn (microsiemens)
	zin (megohm)
	temp (1) 
}



INITIAL {




	vmax = 0
	tvmax = 0
	surf = area
	
}










AFTER SOLVE {








	if (v>vmax) {
		vmax = v
		tvmax = t
	}
}