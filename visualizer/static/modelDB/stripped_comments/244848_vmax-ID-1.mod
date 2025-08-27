NEURON{
	SUFFIX vmax
	RANGE vm
}

PARAMETER{
	v(millivolt)
}

ASSIGNED{
	vm(millivolt)
}

INITIAL{
	vm=v
}

BREAKPOINT{
	if(v>vm){vm=v}
}