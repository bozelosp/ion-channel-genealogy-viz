NEURON {
	SUFFIX body
	POINTER a0Pointer, a1Pointer, a2Pointer
}

PARAMETER {
	tau_m = 2.45  
	umax = 1  
	br = 0.4  
	fsw = 0.0 
}

ASSIGNED { a0Pointer a1Pointer a2Pointer c0 c1 w0 w1 grasperstate }

STATE { u0 u1 xr sw }

BREAKPOINT {
	SOLVE states METHOD euler
}

INITIAL {
	u0 = 0.747647099749367
	u1 = 0.246345045901938
	xr = 0.649984712236374
	sw = 0.0
}

DERIVATIVE states {
	u0' = ((a0Pointer + a1Pointer) * umax - u0) / tau_m
	u1' = ((a2Pointer * umax - u1) / tau_m)
	xr' = (fmusc(u0,u1,xr) + fsw) / br
	grasperstate = (a1Pointer + a2Pointer >= 0.5)
	sw' = -grasperstate * ((fmusc(u0,u1,xr) + fsw * grasperstate) / br)
}

FUNCTION phi(x) {
	
	phi = x - 2.598076211353316 * x * (x * x - 1)
}

FUNCTION fmusc(u0,u1,xr) {
	c0 = 1.0 
	c1 = 1.1 
	w0 = 2 
	w1 = 1.1 
	fmusc = phi((c0 - xr) / w0) * u0 - phi((c1 - xr) / w1) * u1
}