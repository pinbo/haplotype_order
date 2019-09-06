#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  Copyright 2018 Junli Zhang <zhjl86@gmail.com>
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

# example
# ./add_Blosum62_score.py input_file output_file

### Imported
import sys

### read files
in_file = sys.argv[1] # extracted vcf file with #CHROM header
out_file = sys.argv[2]
geno_starts = 10 # vcf v4.2 geno starts from column 10
out = open(out_file, "w")


def gt2snp(allele_list, gt):
	c = "" # snp calls
	if "." in gt: # missing "./." or "."
		c = "N"
	else:
		a, b = gt.split("/")
		if a == b: # homozygous
			c = allele_list[int(a)]
		else:
			c = "H" # heterozygous
	return c

with open(in_file) as file_one:
	#header_line = next(file_one)
	for line in file_one:
		if line.startswith("#"):
			if line.startswith("#CHROM"):
				 header_line = line.strip("#")
			         ll0 = header_line.split("\t")
 			         out.write("\t".join(ll0[0:5] + ll0[(geno_starts - 1):]))
			continue
		line = line.strip()
		if line: # and not line.startswith('#'):
			ll = line.split("\t")
			# convert GT to SNPs
			GTs = ll[(geno_starts - 1):] # all GT calls
			ref = ll[3]
			alt = ll[4].split(",") # there may be more than one alternative alleles
			alleles = [ref] + alt # this way, 0 will be ref, and 1, 2 ... will be alternative allele
			SNPs = [gt2snp(alleles, x) for x in GTs]
			out.write("\t".join(ll[0:5] + SNPs) + "\n")

out.close()
