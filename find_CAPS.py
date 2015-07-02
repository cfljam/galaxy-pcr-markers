#!/usr/bin/python2.6
##find snps that condition CAPS
##usage  find_CAPS.py [-h] -i <reference file>  -g <gff file>


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



import sys
from Bio import SeqIO
from BCBio import GFF
from Bio.Restriction import *
import argparse

###This list is limited to economical enzymes performing well in PCR buffer
rest_batch = RestrictionBatch(
    [AluI, ApaI, BamHI, BbrPI, BfrI, ClaI, DdeI, DpnII, DraI, EcoRI,
     HaeIII, HindII, HinfI, HpaI, PvuII, RsaI, SacI, Sau3AI, SmaI, TaqI])

parser = argparse.ArgumentParser(description='Identify SNPs that condition restriction polymorphisms')
parser.add_argument('-i', type=argparse.FileType('r'), help="input sequence file, required", dest='in_file', required=True)
parser.add_argument('-g', type=argparse.FileType('r'), help="input gff file with SNP and indels, required", dest='gff_file', required=True)
try:
    my_args = parser.parse_args()
except SystemExit:
    print("argument is missing or invalid, exiting...\n")
    sys.exit(0)

in_seq_handle = my_args.in_file
in_gff_handle = my_args.gff_file

out_file=open("find_caps_output.txt",'w')

##use iterator
for myrec in SeqIO.parse(in_seq_handle, "fasta"):
    ##create single-entry dictionary to accept gff annotations from parser
    seq_dict = {myrec.id:myrec}

    ##note that this filters out only SNP features
    limit_info = dict(gff_id = [myrec.id] ,gff_type = ['SNP'])
    in_gff_handle.seek(0)

    ##parse annotations into 
    annotations = [r for r 
                   in GFF.parse(in_gff_handle, 
                                base_dict=seq_dict, 
                                limit_info=limit_info)]

    ##if there are any for this sequence, proceed
    if annotations:
        rec=annotations[0]
        for feat in rec.features:
            fstart=feat.location.start.position
            fend=feat.location.end.position

            if   20 < fstart < len(rec) - 20:
                #just work with +/- 20 bp, ignoring SNPS within this 
                #distance from ends
                fseq=rec.seq[fstart-20:fstart+20]
                ref_seq = rec.seq[fstart-20:fstart+20]
                variant_seq = ref_seq.tomutable()

                #mutate the variant
                variant_seq[20]= feat.qualifiers['Variant_seq'][0]
                variant_seq = variant_seq.toseq()

                #digest the sequences
                ref_cuts =  rest_batch.search(ref_seq)
                var_cuts =  rest_batch.search(variant_seq)

                #print 
                for enz in ref_cuts:
                    kr = set(ref_cuts[enz])
                    km = set(var_cuts[enz])
                    outputstr=[rec.id, fstart +1,fend+1,feat.id,enz]
                    if len(kr) > len(km):
                        outputstr.append("reference")
                        print('\t'.join(map(str,outputstr)))
			out_file.write('\t'.join(str(x) for x in outputstr)+'\n')
                    elif len(kr) < len(km):
                        outputstr.append("variant")
                        print('\t'.join(map(str,outputstr)))
			out_file.write('\t'.join(str(x) for x in outputstr)+'\n')

in_gff_handle.close()
in_seq_handle.close()                                                                              
out_file.close()


