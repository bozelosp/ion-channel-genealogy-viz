NEURON { SUFFIX nothing }

VERBATIM
char* secname();
ENDVERBATIM

PROCEDURE scale_connection_coef(x, factor) {
VERBATIM {
	Section* sec;
	Node* nd;
#if defined(t)
	_NrnThread* _nt = nrn_threads;
#endif
	sec = chk_access();
	if (_lx <= 0. || _lx > 1.) {
		hoc_execerror("out of range, must be 0 < x <= 1", (char*)0);
	}
	
	
	if (_lx == 1.) {
		nd = sec->pnode[sec->nnode-1];
	}else{
		nd = sec->pnode[(int) (_lx*(double)(sec->nnode-1))];
	}
	
	NODEA(nd) *= _lfactor;
	NODEB(nd) *= _lfactor;
}
ENDVERBATIM
}