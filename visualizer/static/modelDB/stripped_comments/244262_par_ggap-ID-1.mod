NEURON {
	POINT_PROCESS gGapPar
	RANGE g, i, vgap
    NONSPECIFIC_CURRENT	i
}
PARAMETER { g = 1e-10 (1/megohm) }
ASSIGNED {
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
}
BREAKPOINT { i = (v-vgap)*g }