#!/usr/local/bin/python2.6
"""
Usage vcf_gff.py <vcf_file input> <gff_file output>

Input vcf file format:
   CHROM POS ID REF ALT QUAL FILTER INFO FORMAT

Note: Generating vcf from a single merged bam file, using multiple bam files results in multiple FORMAT columns!!!

Output gff format:
    SEQID SOURCE TYPE START END SCORE STRAND PHASE ATTRIBUTES

"""
#Copyright 2012 Susan Thomson
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
                                                                              
in_vcf_file = open(sys.argv[1], 'r')
out_gff_file = open(sys.argv[2], 'w')

def get_info(attribute_input):
    """
    Get the record_type by determining if INDEL is stated in column INFO; get
    raw read count
    """
    INFO = {}
    rec = attribute_input.split(";")
    if "INDEL" in rec:
        record_type = "INDEL"
        rec = rec[1:]
    else:
        record_type = "SNP"
    for entry in rec:
        detail = entry.split("=")
        if len(detail) < 2:
            continue
        INFO[detail[0]] = detail[1]
    if INFO.has_key("DP"):
        reads = INFO.get("DP")
    else:
        reads = "NA"
    data = (record_type, reads)
    return data

def get_gen(formatcols, ref):
    """
    Get info on heterozyosity or homozygosity (could be useful later),
    by estimating genotype(GT) calling ref allele=0, variant allele=1

    """
    formats = []
    sample_dict = {}
    genos = ""
    format_types = formatcols[0].split(":")
    samples = formatcols[1:]
    for entry in format_types:
        formats.append(entry)
    for sample in samples:
        var = ""
        data = sample.split(":")
        sample_dict = dict(zip(formats, data))
        if sample_dict.has_key("DP"):
            reads = sample_dict.get("DP")
        else:
            reads = "NA"
        if sample_dict.has_key("GT"):
            """
            mpileup output, recommend ALWAYS have GT, note this only good for scoring diploids too!
            """
            genotypes = sample_dict.get("GT")
            if genotypes == "1/1":
                gen = "HOM_mut"
            if genotypes == "0/1":
                gen = "HET"
            if genotypes == "0/0":
                gen = "HOM_ref"
        try: # set gen to 'NA' if still unset
            gen
        except NameError:
            gen = "NA"
        geno = ("%s:%s " % (reads, gen))
        genos += geno
        sample_dict = {}
    return genos
    
attributes = {}
"""
Get relevant info from vcf file and put to proper gff columns
"""

out_gff_file.write("#gff-version 3\n")
for line in in_vcf_file:
    if line.startswith("#") == False:
        info = line.split()
        seqid = info[0].strip()
        source = "SAMTOOLS"
        start = int(info[1])
        score = info[5]
        strand = "."
        phase = "."
        reference = info[3].strip()
        variant = info[4].strip()
        attr = get_info(info[7])
        record_type = attr[0]
        reads = attr[1]
        if record_type == "INDEL":
            end = start + len(reference)
        else:
            end = start 
        gen = get_gen(info[8:], reference)
        out_gff_file.write(
            ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tID=%s:%s:%s:%d;Variant" +
             "_seq=%s;Reference_seq=%s;Total_reads=%s;Zygosity=%s\n") %
            ( seqid, source,record_type, start, end, score, strand, phase,seqid, 
              source, record_type, start, variant, reference, reads, gen))

out_gff_file.close()    


