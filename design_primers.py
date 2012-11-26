#!/usr/local/bin/python2.6
##design primers to features in multiple sequences
##usage: python  design_primers.py <fasta-file> <gff file> <file of target IDs> <prod_min_size> <prod_max_size>

##CAUTION will only reliably work with  Primer3 version 1.1.4 or earlier ##


#Copyright 2012 John McCallum & Leshi Chen
#New Zealand Institute for Plant and Food Research
#This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import StringIO
import re
import tempfile
import subprocess
import copy
import sys
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Emboss import Primer3


in_file = sys.argv[1]
gff_file = sys.argv[2] 
target_file =  sys.argv[3]
prod_min_size = int(sys.argv[4])
prod_max_size = int(sys.argv[5])

max_tm_diff=1                                        ##
opt_GC_percent=50                                    ##
maxpolyx=4                                           ##
optimum_length=20
##target is specified in start, end format 
productsizerange = str(prod_min_size) + "," + str(prod_max_size)
#open input files
in_seq_handle = open(in_file)
in_gff_handle = open(gff_file)
in_target_handle=open(target_file)
#read  target feature IDs into list
targets=in_target_handle.readlines()
in_target_handle.close()
##and create a hit list of sequences from this
target_seq_id_list = list(set([line.split(":")[0] for line in targets]))
##create iterator returning sequence records
for myrec in SeqIO.parse(in_seq_handle, "fasta"):
    #check if this sequence is included in the target list
    if myrec.id in target_seq_id_list:
        ##create sequence dictionary so we can add in gff annotations
        seq_dict = {myrec.id : myrec}
        ##just limit to gff annotations for this sequence
        limit_info = dict(gff_id = [ myrec.id ])
        ##rewind gff filehandle
        in_gff_handle.seek(0)
        ##read annotations into sequence dictionary for this sequence record only 
        annotations = [r for r in GFF.parse(in_gff_handle, base_dict=seq_dict,limit_info=limit_info)]
        ##if there are any annotations, then  proceed. 
        if annotations:
            rec=annotations[0]
            ##iterate over list of target IDs
            for t in targets:
                target_ID = t.strip('\n')
                target_annotations = [f for f in rec.features if f.id == target_ID ]
                if target_annotations:
                    mytarget = target_annotations[0]
                    #create temporary files
                    tempfastaFile = tempfile.mktemp()
                    tempproutfile = tempfile.mktemp()
                    #just consider slice of sequence in a window of +/- prod_max_size  bp
                    ##from feature UNLESS feature is close to end
                    ##Note that slice is zero-based
                    featLocation = mytarget.location.start.position 
                    if featLocation > prod_max_size:
                        slice_start = featLocation - prod_max_size
                        featPosition = prod_max_size  
                    else:
                        slice_start = 0
                        featPosition = featLocation
                    if (len(rec) - featLocation) < prod_max_size:
                        slice_end = len(rec)
                    else:
                        slice_end = featLocation + prod_max_size
                    ###grab slice of sequence fom this window.
                    targetRec = rec[slice_start:slice_end]
                    matching_feature = [f for f in targetRec.features if f.id == mytarget.id]
                    if matching_feature:
                        target_feat = matching_feature[0]
                        if target_feat.location.start.position == 0:
                            target_feat.location.start.position = 1
                        #we get the mask features by removing the target...all features are masked as just using snp and indels
                        ##a smarter filter could be added 
                        ##note use of list copy to avoid possible side-effects
                        exclude_feat = list(targetRec.features)
                        exclude_feat.remove(target_feat)
                        ##print'targetRec.features',  targetRec.features ##for debug
                        mask_str=map(lambda f: str(f.location.start.position+1) + "," + str(f.location.end.position + 1) ,exclude_feat)
                        #mask_str=map(lambda f: str(f.location).strip('[]'),exclude_feat)
                        p3_exclude_str = str(mask_str).replace('\', \'',':')
                        p3_target = str(target_feat.location.start.position +1)  + "," + str(target_feat.location.end.position +1)
                        #write sequence record into template file as  fasta
                        t_output_handle = open(tempfastaFile, "w")
                        SeqIO.write([targetRec], t_output_handle, "fasta")
                        t_output_handle.close()
                        #create Primer3Commandline() for emboss tool
                        primer_cl = Primer3Commandline()
                        #set the emboss tool to suppress  output as this will make Galaxy  think it is error message although it is a message to state run success
                        primer_cl.set_parameter("-auto",'1')
                        ###pass  sequence file to emboss. FAILS for primer3_core V2 ##
                        primer_cl.set_parameter("-sequence",tempfastaFile)
                        #add target location
                        primer_cl.set_parameter("-target", p3_target)
                        ##mask off other features...dumb masking of everything at present, beware
                        if (p3_exclude_str != ""):
                            primer_cl.set_parameter("-excludedregion", p3_exclude_str)
                        #add temporary output file to get the result
                        primer_cl.set_parameter("-outfile", tempproutfile)
                        #specify maximum different of tm
                        primer_cl.set_parameter("-maxdifftm",max_tm_diff )
                        #other useful parameters
                        primer_cl.set_parameter("-ogcpercent", opt_GC_percent)
                        primer_cl.set_parameter("-opolyxmax", maxpolyx)  
                        primer_cl.set_parameter("-osize", optimum_length)
                        #set product size range
                        primer_cl.set_parameter("-prange", productsizerange)
                        #using python subprocess method to run emboss command line programe with the parameters given
                        fnull = open(os.devnull, 'w')
                        result=subprocess.check_call(str(primer_cl),shell=True ,stdout = fnull, stderr = fnull)
                        #read temporary outputfile
                        handle = open(tempproutfile)
                        record = Primer3.read(handle)
                        ##just return first set, if there is one
                        if len(record.primers) > 0:
                            primer= record.primers[0]
                            outputstr=[mytarget.id, primer.forward_seq,primer.reverse_seq,primer.size]
                        else:
                            outputstr=[mytarget.id,"NONE","NONE","NONE"]
                        print('\t'.join(map(str,outputstr)))

                        
in_gff_handle.close()
in_seq_handle.close()
