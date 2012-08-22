#!/usr/bin/perl
#parse_primersearch.pl
#reformat EMBOSS primersearch output into columnar Galaxy interval format 

#Copyright 2012 John McCallum 
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

open (IN, "<$ARGV[0]");
open (OUT, ">$ARGV[1]");

#print OUT  "primerset_id","\t","sequence_id","\t","hit_start","\","mismatches","\t","amplimer_size",\n";



while (<IN>) {
         /^Primer name (\S+)/  && ($name = $1);  # get primer set name
         # Modified to cope with unnamed sequence input 28/7/05
        /Sequence: (\S+)/ && print OUT  $name,"\t",$1;
        /Sequence:(\s{4,})/ && print OUT $name,"\t","unnamed_seq";
        /hits forward strand at (\d+) with (\d) mismatches/ && ($start = $1) &&  print OUT  "\t",$2,"\t",$start,;
        /Amplimer length: (\S+)/ && ($amp_length = $1) &&  print OUT "\t",$start + $amp_length,"\t",$1,"\n";
        }

close( IN );
close( OUT );
