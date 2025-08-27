NEURON {
 	SUFFIX nothing
 }
 
 VERBATIM
 static subtree(sec, sym) Section* sec; Symbol* sym; {
 	Section* child;
 
 	nrn_pushsec(sec);	
 	hoc_run_stmt(sym);	
 	nrn_popsec();
 
 	for (child = sec->child; child; child = child->sibling) {
 		subtree(child, sym);
 	}
 }
 ENDVERBATIM
 
 PROCEDURE subtree_traverse_all() {
   VERBATIM
   {
 	Section* chk_access();
 	Symbol* hoc_parse_stmt();
 	Symlist* symlist = (Symlist*)0;
 	subtree(chk_access(), hoc_parse_stmt(gargstr(1), &symlist));
 	
 	hoc_free_list(&symlist);
   }
   ENDVERBATIM
 }