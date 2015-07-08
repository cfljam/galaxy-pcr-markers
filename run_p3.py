#!/usr/bin/python
##run primer3 by passing Python dictionary

#run_P3.py

import primer3

##call P3 with dict of args, returns dict, no exception handling
def run_P3(target_dict,global_dict):
    P3_dict=primer3.bindings.designPrimers(target_dict,global_dict)
    ##return iterable list
    primer_list=[dict(PRIMER_RIGHT_SEQUENCE=P3_dict.get('PRIMER_RIGHT_'+ str(X) + '_SEQUENCE'),\
    PRIMER_LEFT=P3_dict.get('PRIMER_LEFT_'+ str(X) ),\
    PRIMER_RIGHT=P3_dict.get('PRIMER_RIGHT_'+ str(X) ),\
    PRIMER_LEFT_SEQUENCE=P3_dict.get('PRIMER_LEFT_'+ str(X) + '_SEQUENCE')) for X in range(0,int(P3_dict.get('PRIMER_RIGHT_NUM_RETURNED'))-1)]
    return(primer_list)
