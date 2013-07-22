#!/usr/bin/pyrthon
##run primer3 by passing Python dictionary

#un_P3.py

import subprocess as sp
import copy


if "check_output" not in dir(sp): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = sp.Popen(stdout=sp.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise sp.CalledProcessError(retcode,cmd)
        return output
    sp.check_output = f




##call P3 with dict of args, returns dict, no exception handling 
def run_P3(target_dict):
    p3_str=''
    for key in target_dict:
        p3_str+=key
        p3_str+='='
        p3_str+=str(target_dict[key])
        p3_str+='\n'
    p3_str+='='
    input_str='echo -e \"' + p3_str + '\" | primer3_core '
###exception handling to be added here
    output = sp.check_output(input_str,shell=True)
    output_fields=output.split('\n')
    ##put output into a dict, omitting trailing =
    P3_dict=dict([X.split('=') for X in output_fields][:len(output_fields)-2])
    ##return iterable list
    primer_list=[dict(PRIMER_RIGHT_SEQUENCE=P3_dict.get('PRIMER_RIGHT_'+ str(X) + '_SEQUENCE'),PRIMER_LEFT=P3_dict.get('PRIMER_LEFT_'+ str(X) ),PRIMER_RIGHT=P3_dict.get('PRIMER_RIGHT_'+ str(X) ),PRIMER_LEFT_SEQUENCE=P3_dict.get('PRIMER_LEFT_'+ str(X) + '_SEQUENCE')) for X in range(0,int(P3_dict.get('PRIMER_RIGHT_NUM_RETURNED'))-1)]
    return(primer_list)




