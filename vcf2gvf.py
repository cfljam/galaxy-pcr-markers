#!/usr/bin/env python
"""
Usage vcf2gvf.py <vcf_file input> <gff_file output>
"""
import sys,re

##simplify the variant descriptor for GVF
def get_vartype(type_str):
    type_str=type_str.upper()
    for V in ['COMPLEX', 'MNP', 'DEL', 'INDEL','INS','SNP']:
        if V in type_str:
            return V
            break


                                                                              
in_vcf_file = open(sys.argv[1], 'r')
out_gff_file = open(sys.argv[2], 'w')
out_gff_file.write("#gff-version 3\n") ####gvf-version 1.06

##only deal with tidy source lines, samtools is special needs
source_pat = re.compile('##source=(\S+)', re.IGNORECASE)


while True:
    line=in_vcf_file.readline()
    if line[0]!='#':
        break
    elif line[2:10]=='samtools':
        source_type='SAMTOOLS'
        break
    else:
        m=source_pat.match(line)
        if m:
            source_type=m.group(1).upper()
            break


##now read the data
##samtools only distinguishes SNPs and indels
##This is much complicated in highly heterozygous systems
##see http://www.sequenceontology.org/miso/current_release/term/SO:0001059
##GVF column 3 requires a type from this ontology
###vcf type may be snp, mnp, ins, del, or complex.
##for current purposes map multiple types to most complex state
##ie SNP,SNP -> SNP; SNP,MNP-> MNP; SNP,complex 
try:
    for line in in_vcf_file:
        if line[0] != '#':
            var=line.split() ##see if var[7] has indel at start. If so add TYPE=INDEL for samtools case
            info_field=var[7].split(';')
            if info_field[0]=='INDEL':
                info_field[0]="TYPE=INDEL"
            info_dict=dict(X.split('=') for X in info_field)
            ##if no TYPE key, then add TYPE=SNP for samtools case
            if not info_dict.has_key('TYPE'):
                info_dict['TYPE']='SNP';
            var_type=get_vartype(info_dict['TYPE'])
            ID=":".join([var[0],source_type,var_type,var[1]])
            start=int(var[1])
            length=len(var[3]) ## using reference length in this case
            attributes=";".join(['ID='+ID,'Reference_seq='+var[3],'Variant_seq='+var[4]])
            output_line=[var[0], source_type, var_type,  str(start), str(start+length-1) ,'.','.','.',attributes,"\n"]
            out_gff_file.write("\t".join([str(X) for X in output_line]))
finally:
    in_vcf_file.close()
    out_gff_file.close()

