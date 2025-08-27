NEURON {
        SUFFIX nothing
}

VERBATIM
static subtree(sec, sym) Section* sec; Symbol* sym; {
        Section* child;

 

        for (child = sec->child; child; child = child->sibling) {
       nrn_pushsec(child);       
        hoc_run_stmt(sym);      
        nrn_popsec(); 

        }
}
ENDVERBATIM

PROCEDURE subtree_traverse() {
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