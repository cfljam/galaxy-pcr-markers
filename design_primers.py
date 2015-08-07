
#!/usr/bin/env python
##design primers to features in multiple sequences, with option to predict melting
#usage: design_HRM_primers.py [-h] -i IN_FILE -g GFF_FILE -T TARGET_FILE [-u]
#                             [-n MAX_PRIMERS] [-p PROD_MIN_SIZE]
#                             [-P PROD_MAX_SIZE] [-l OPT_PRIMER_LENGTH]
#                             [-m MAX_TM_DIFF] [-t OPTIMUM_TM]
#                             [-G OPT_GC_PERCENT] [-x MAXPOLYX] [-c GC_CLAMP]

#Copyright 2013 John McCallum & Susan Thomson
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
import copy
import sys
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
import run_p3 as P3
#import umelt_service as umelts
import argparse

##Primer3 defaults or additional options defined as dictionary
def_dict={
'PRIMER_MIN_SIZE':18 ,
'PRIMER_MAX_SIZE':25,
'PRIMER_MAX_NS_ACCEPTED':1}

#parse arguments
parser = argparse.ArgumentParser(description='Primer set design and melt prediction parameters')
parser.add_argument('-i', type=argparse.FileType('r'), help="input sequence file, required", dest='in_file', required=True)
parser.add_argument('-g', type=argparse.FileType('r'), help="input gff file with SNP and indels, required", dest='gff_file', required=True)
parser.add_argument('-T', type=argparse.FileType('r'), help="input target SNP file, required", dest='target_file', required=True)
parser.add_argument('-u',  help="do uMelt prediction, optional", dest='run_uMelt',action='store_true', default=False )
parser.add_argument('-n', type=int, help="maximum number of primer pairs to return, default=5", dest='max_primers', default=5) ## PRIMER_NUM_RETURN
parser.add_argument('-p', type=int, help="minimum product size", dest='prod_min_size', default=100)                             ## PRIMER_PRODUCT_SIZE_RANGE min
parser.add_argument('-P', type=int, help="maximum product size", dest='prod_max_size', default=300)                            ## PRIMER_PRODUCT_SIZE_RANGE max
parser.add_argument('-l', type=int, help="optimum primer length", dest='opt_primer_length', default=20)                        ## PRIMER_OPT_SIZE
parser.add_argument('-m', type=int, help="maximum tm difference between primers", dest='max_tm_diff', default=1)               ## PRIMER_PAIR_MAX_DIFF_TM
parser.add_argument('-t', type=int, help="optimum Tm for primers, recommend range 59 to 61", dest='optimum_tm', default=59)    ## PRIMER_OPT_TM
parser.add_argument('-G', type=int, help="optimum GC percentage of primers", dest='opt_GC_percent', default=50)                ## PRIMER_OPT_GC_PERCENT
parser.add_argument('-x', type=int, help="maximum polyx, recommend less than 4", dest='maxpolyx', default=3)                   ## PRIMER_MAX_POLY_X
parser.add_argument('-c', type=int, help="number of C/Gs at end, recommend 2", dest='gc_clamp', default=1)                     ## PRIMER_GC_CLAMP

parser.add_argument('-e', type=int, help="maximum allowable 3'-anchored complementarity", dest='maxselfend', default=3)                     ## PRIMER_MAX_SELF_END
parser.add_argument('-a', type=int, help="maximum complementarity between left and right or self", dest='maxselfany', default=8)                     ## PRIMER_MAX_SELF_ANY
parser.add_argument('-maxgc', type=float, help="Maximum allowable percentage of Gs and Cs in any primer.", dest='maxgc', default=80.0)                     ## PRIMER_MAX_GC
parser.add_argument('-mingc', type=float, help="Minimum allowable percentage of Gs and Cs in any primer.", dest='mingc', default=20.0)                     ## PRIMER_MIN_GC


parser.add_argument('-d', type=str, help="variant indentifier delimiter, used to separate sequence ID from rest ", dest='target_delim', default=':')
try:
        my_args = parser.parse_args()
except SystemExit:
        print("\nOops, an argument is missing/invalid, exiting...\n")
        sys.exit(0)

##update from args. NEEDS TO BE FINISHED
productsizerange = [my_args.prod_min_size,my_args.prod_max_size]
def_dict['PRIMER_PRODUCT_SIZE_RANGE']=productsizerange
def_dict['PRIMER_NUM_RETURN']=my_args.max_primers +1
def_dict['PRIMER_OPT_SIZE']=my_args.opt_primer_length
def_dict['PRIMER_PAIR_MAX_DIFF_TM']=my_args.max_tm_diff
def_dict['PRIMER_OPT_TM']=my_args.optimum_tm
def_dict['PRIMER_OPT_GC_PERCENT']=my_args.opt_GC_percent
def_dict['PRIMER_MAX_POLY_X']= my_args.maxpolyx
def_dict['PRIMER_GC_CLAMP']=my_args.gc_clamp

def_dict['PRIMER_MAX_SELF_END']=my_args.maxselfend
def_dict['PRIMER_MAX_SELF_ANY']=my_args.maxselfany
def_dict['PRIMER_MAX_GC']=my_args.maxgc
def_dict['PRIMER_MIN_GC']=my_args.mingc



##conditional import of umelt
if my_args.run_uMelt:
    import umelt_service as umelts

#open input files

targets=[line.rstrip() for line in my_args.target_file.readlines()]
my_args.target_file.close()
##and create a hit list of sequences from this
target_seq_id_list = [re.split(my_args.target_delim,X)[0] for X in targets] ## target_delimiter defaults to ':'  e.g. ABC:SNP:SAMTOOL:1234

##print header
print "SNP_Target_ID", "Position","Ref_base","Variant_base" ,"Amplicon_bp","PRIMER_LEFT_SEQUENCE",'PRIMER_RIGHT_SEQUENCE', "ref_melt_Tm","var_melt_Tm","Tm_difference"
##create iterator returning sequence records
for myrec in SeqIO.parse(my_args.in_file, "fasta"):
    #check if this sequence is included in the target list
    if myrec.id in target_seq_id_list:
        ##create sequence dictionary so we can add in gff annotations
        seq_dict = {myrec.id : myrec}
        ##just limit to gff annotations for this sequence
        limit_info = dict(gff_id = [ myrec.id ])
        ##rewind gff filehandle
        my_args.gff_file.seek(0)
        ##read annotations into sequence dictionary for this sequence record only
        annotations = [r for r in GFF.parse(my_args.gff_file, base_dict=seq_dict,limit_info=limit_info)]
        ##if there are any annotations, then  proceed.
        if annotations:
            rec=annotations[0]
            ##iterate over list of target IDs
            for t in targets:
                target_ID = t.strip('\n')
                target_annotations = [f for f in rec.features if f.id == target_ID ]
                if target_annotations:
                    mytarget = target_annotations[0]
                    #just consider slice of sequence in a window of +/- prod_max_size  bp
                    ##from feature UNLESS feature is close to end
                    ##Note that slice is zero-based
                    featLocation = mytarget.location.start.position
                    if featLocation > my_args.prod_max_size:
                        slice_start = featLocation -  my_args.prod_max_size
                        featPosition =  my_args.prod_max_size
                    else:
                        slice_start = 0
                        featPosition = featLocation
                    if (len(rec) - featLocation) <  my_args.prod_max_size:
                        slice_end = len(rec)
                    else:
                        slice_end = featLocation +  my_args.prod_max_size
                    ###grab slice of sequence fom this window.
                    targetRec = rec[slice_start:slice_end]
                    matching_feature = [f for f in targetRec.features if f.id == mytarget.id]
                    if matching_feature:
                        target_feat = matching_feature[0]
                        my_target_dict={} # re-initialize target dictionary
                        if target_feat.location.start.position == 0:
                            target_feat.location.start.position = 1
                        #get the mask features by removing  target...all features are masked as just using snp and indels, a smarter filter could be added
			exclude_feat = list(targetRec.features) ##list copy to avoid possible side-effects
                        exclude_feat.remove(target_feat)
                        my_target_dict={'SEQUENCE_ID' : rec.name,\
                         'SEQUENCE_TEMPLATE': targetRec.seq.tostring().upper(),\
                         'SEQUENCE_TARGET': [target_feat.location.start.position,1],\
                         'SEQUENCE_EXCLUDED_REGION': [[x.location.start.position,x.location.end.position -x.location.start.position] for x in exclude_feat]}
 			result=P3.run_P3(target_dict=my_target_dict,global_dict=def_dict)
                        if my_args.run_uMelt:
                            amp_seq=targetRec.seq ##need to make this conditional on getting a result >0 and melt=True
                            mutamp_seq=targetRec.seq.tomutable()
                            mutamp_seq[target_feat.location.start:target_feat.location.end]=target_feat.qualifiers['Variant_seq'][0] #mutate to variant
                            other_SNP=[f for f in targetRec.features if f.id != mytarget.id]
                            if other_SNP:
                                for snp in other_SNP:
                                    mutamp_seq[snp.location.start:snp.location.end]=snp.qualifiers['Variant_seq'][0]
			for primerset in result:
				amp_start=int(primerset['PRIMER_LEFT'][0])
				amp_end=int(primerset['PRIMER_RIGHT'][0])
                                ref_melt_Tm=0
                                var_melt_Tm=0
                                diff_melt=0
                                if my_args.run_uMelt:
                                    try:
                                        ref_melt_Tm=umelts.getTm(umelts.getmelt(amp_seq.tostring()[amp_start:amp_end+1]))
					var_melt_Tm=umelts.getTm(umelts.getmelt(mutamp_seq.tostring()[amp_start:amp_end+1]))
                                        diff_melt=abs(ref_melt_Tm - var_melt_Tm)
                                    except:
					ref_melt_Tm="NA" ##preferably something more informative?
					var_melt_Tm="NA" ##exception handling to be added
                                        diff_melt="NA"
                                if target_feat.qualifiers.has_key('Reference_seq'):
                                    reference_seq=target_feat.qualifiers['Reference_seq'][0]
                                else:
                                    reference_seq="NA"
                                if target_feat.qualifiers.has_key('Variant_seq'):
                                    variant_seq=target_feat.qualifiers['Variant_seq'][0]
                                else:
                                    variant_seq="NA"
                                print mytarget.id, featLocation + 1 ,reference_seq, variant_seq,\
                                amp_end-amp_start,primerset['PRIMER_LEFT_SEQUENCE'],\
                                primerset['PRIMER_RIGHT_SEQUENCE'], ref_melt_Tm,var_melt_Tm,diff_melt#, amp_seq.tostring()[amp_start:amp_end+1], mutamp_seq.tostring()[amp_start:amp_end+1]

my_args.gff_file.close()
my_args.in_file.close()
