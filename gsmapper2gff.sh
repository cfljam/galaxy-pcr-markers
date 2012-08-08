#!/bin/sh
##convert gsMapper output into gff3/GVF format

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


infile=$1
outfile=$2

awk '
BEGIN {OFS="\t"}
/^>/  && sub(/%/,"",$7)  {
  ID=substr($1,2)
  if  (length($4) > 1 || match($4,"-") || length($5) > 1 || match($5,"-")) 
    type="indel"
    else
      type="SNP"
start=$2
end=$3
Col9_ID=ID ":gsmapper:" type ":"start

Reference_seq=$4
Variant_seq=$5
Total_reads=$6
Variant_reads=Total_reads * $7 /100 - (Total_reads * $7 % 100)/100 



 print ID,"gsmapper",type,start,end,".",".",".","ID="Col9_ID";Reference_seq="Reference_seq";Variant_seq="Variant_seq";Total_reads="Total_reads";Variant_reads="Variant_reads
}' "$infile" > "$outfile" 








