# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

import commands
import sys

if len(sys.argv) != 5:
	print "Usage: python msa2block.py msa_fname perl_dir fasta_fname outfname"
	sys.exit(1)
# end if 

msa_fname   = sys.argv[1]
fasta_fname = sys.argv[2]
perl_dir    = sys.argv[3]
outfname    = sys.argv[4]



def read_fasta(fname):
	lines = open(fname).read().strip().split('\n')
	lines = lines[1:] # skip '>'
	seq = ''.join(lines)
	return seq
# end def


# check number of iteration, get the last's
fin = open(msa_fname)
for line in fin:
	if 'Results from round' in line:
		iter_n = line.split()[-1]
	# end if
# end for
fin.close()

#outfname = outfname+'.'+iter_n+'.blk'
# generate command to run perl script (align to block form)
cmd = ' '.join(['perl ', perl_dir.strip()+'/msa2block.pl', msa_fname, outfname, iter_n])
print cmd
commands.getoutput(cmd)

# read fasta file and write it on the top of block
seq = read_fasta(fasta_fname)
content = open(outfname).read()
content = seq + '\n' + content
fout = open(outfname, 'w')
fout.writelines(content)
fout.close()
