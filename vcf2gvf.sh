#!/bin/sh
##convert vcf to gvf
##NOTE This is a very simple basic parser for a complex format.
#It is intended for use with mpileup output where -g or -u flags are NOT used.

##usage vcf2gvf.sh <vcf file> <outputfile>

#Copyright 2012 John McCallum & Leshi Chen
#New Zealand Institute for Plant and Food Research

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



inputfile=$1
outputfile=$2

echo  "##gvf-version 1.05" > $outputfile

awk '
BEGIN {OFS="\t"}

##get feature type 
{if (index($8,"INDEL")== 1) {type="INDEL"} else {type="SNP"} }
##get feature length
{if (type=="SNP") 
    {feat_length=1}
    else {feat_length=length($4)} 
}
{end=($2+feat_length)}

!/^#/  { print $1 ,"SAMTOOLS",type,$2,end,$6,".",".","ID="$1":SAMTOOLS:"type":"$2";Variant_seq="$5";Reference_seq="$4";"$8}

END {print ""} 
' "$inputfile" > "$outputfile"