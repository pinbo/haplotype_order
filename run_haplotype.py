#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  run_getkasp.py
#  
#  Copyright 2017 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# Just to run all the steps at once for designing KASP primers
# Need to use in MacPro

# change the reference location accordingly.

# example: run_getkasp.py for_polymarker.csv 3 200 1 1 1 65

#########################
from glob import glob


def main(args):
	chrom = args[1]
	from_bp =  args[2]
	to_bp = args[3]
	vcf_file = args[4] # vcf.gz
	meta_data = args[5] # metadata for lines to merge to SNP table

	script_path = os.path.dirname(os.path.realpath(__file__)) + "/bin/" # scripts folder
	vcf_file = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
	
	# step 1: extract SNPs
	cmd0 = "gzcat " + vcf_file + " | gawk '!/^#/{exit} /^#CHROM/'> test.txt"
	cmd1 = "gzcat " + vcf_file + ' | gawk \'$1 == "' + chrom + '" &&  $2 > ' + to_bp + ' {exit} $1 == "' + chrom +  '" &&  $2 >= ' + from_bp + '" && $2 <= ' + to_bp + "' >> test.txt"
	print "Step 1: Extract SNPs\n", cmd0, "\n", cmd1
	call(cmd1, shell=True)
	
	#step 2: convert from GT to SNP alleles
	cmd2 = script_path + "convert_vcf_calls_to_SNP.py test.txt out.txt"
	print "Step 2: Convert vcf to SNP table\n", cmd2
	call(cmd2, shell=True)
	
	# Step 3: sort lines by similarity for haplotype analysis
	cmd3 = script_path + "haplotype_analysis.R out.txt " + meta_data
	print "Step 3: Get the sorted SNP table for haplotype analysis:\n", cmd3
	call(cmd3, shell=True)
	
	print "\n\n\n Finish successfully! Output data is 'haplotype_ordered_genotable.txt' "
	return 0

if __name__ == '__main__':
	import sys, os
	from subprocess import call
	sys.exit(main(sys.argv))
