NEURON {
	ARTIFICIAL_CELL IF3
	RANGE taue, taui, taum, e, i, m, b
	RANGE nself, nexcite, ninhibit
	GLOBAL eps, refrac
}

PARAMETER {
	
	taue = 3 (ms) < 1e-9, 1e9 >
	taui = 10 (ms) < 1e-9, 1e9 >
	taum = 30 (ms) < 1e-9, 1e9 >
	b = 0
	eps = 1e-6
	refrac = 5 (ms)
}

ASSIGNED {
	e i m
	enew inew mnew
	t0 (ms)
	nself nexcite ninhibit

	ae ai
	be bi
	on
}


PROCEDURE newstates(d(ms)) {
	LOCAL ee, ei, em

	ee = exp(-d/taue)
	ei = exp(-d/taui)
	em = exp(-d/taum)

	enew = e*ee
	inew = i*ei
	mnew = b + (m - b)*em
		+ ae*e*(taue/(taum - taue))*(em - ee)
		+ ai*i*(taui/(taum - taui))*(em - ei)
}

FUNCTION M() {
	if (on) {
		newstates(t - t0)
		M = mnew
	}else{
		M = 0
	}
}

FUNCTION E() {
	newstates(t - t0)
	E = enew
}

FUNCTION I() {
	newstates(t - t0)
	I = inew
}

PROCEDURE update() {
	e = enew
	i = inew
	m = mnew
	t0 = t
}

PROCEDURE factors() {
	LOCAL tp
	
	
	
	if (taue >= taui) { taui = taue + 0.01 }
	if (taum == taue) { taum = taue + 0.01 }
	if (taum == taui) { taum = taui + 0.01 }

	
	tp = log(taue/taum)/((1/taum) - (1/taue))
	be = 1/(exp(-tp/taum) - exp(-tp/taue))
	ae = be*((taum/taue) - 1)

	
	tp = log(taui/taum)/((1/taum) - (1/taui))
	bi = 1/(exp(-tp/taum) - exp(-tp/taui))
	ai = bi*((taum/taui) - 1)
}

FUNCTION tf(bb) (ms) {
	if (bb > 1) { tf = taum*log((bb-m)/(bb-1)) }
	else { tf = 1e9 }
}

FUNCTION firetimebound() (ms) {
	LOCAL h, temp
	h = ae*e + ai*i
	if (b>1) {
		if (h>0) {
			firetimebound = tf(h+b)
		} else {
			
			temp = tf(b) 
			temp = tf(b + h*exp(-temp/taui)) 
			
			firetimebound = tf(b + h*exp(-temp/taui))
		}
	} else {
		

		
		if (h+b > 1) {
			firetimebound = tf(h+b)
		} else {
			firetimebound = 1e9
		}
	}

}

LOCAL total, nexttotal

INITIAL {
	factors()
	e = 0
	i = 0
	m = 0
	t0 = t
	on = 1
	net_send(firetimebound(), 1)

	nself=0
	nexcite=0
	ninhibit=0
	total = 0
	nexttotal = 1000
}

NET_RECEIVE (w) {
	newstates(t-t0)
	update()	

	if (flag == 1) { 
		nself = nself + 1
		if (m > 1-eps) { 
			if (m > 1+eps) { 
				printf("m>1 error in IF3 mechanism--m = %g\n", m)
				
			}
			
			net_event(t)
			on = 0
			m = 0
			net_send(refrac, 2) 
			total = total + 1
			if (total == nexttotal) {

				nexttotal = nexttotal + 1000
			}
		}else{ 
			net_send(firetimebound(), 1)
		}
	}else if (flag == 2) { 
		on = 1
		m = 0
		net_send(firetimebound(), 1)
	} else {
		if (w > 0) {
			nexcite = nexcite + 1
			e = e + w
		} else {
			ninhibit = ninhibit + 1
			i = i + w
		}

		if (on) {
			net_move(firetimebound() + t)
		}
	}
}