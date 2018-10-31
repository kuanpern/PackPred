# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# this script automate the process of getting entropy values for a pdb file
import commands
import sys
fname               = sys.argv[1]
final_entropy_fname = sys.argv[2]
blosum62            = sys.argv[3] #'../database/BLOSUM62'
bin_dir             = sys.argv[4] #'../'
blast_db            = sys.argv[5] #'../database/uniref90.db'
blast_exe           = sys.argv[6] #'../blast-2.2.22/bin/blastpgp'
blast_iter_n        = sys.argv[7] #'5'
e_thr               = sys.argv[8] #'0.0001'
max_seq_n           = sys.argv[9] #'1000'


# 1. get fasta sequence for every chain
inputs = ['python', bin_dir+'/pdb2fasta.py', fname, fname]
cmd = ' '.join(inputs)
print cmd
fastas = commands.getoutput(cmd).strip().split('\n')

# 2. blast for all chains
# copy blosum62 table here
inputs = ['cp', blosum62, '.']
cmd = ' '.join(inputs)
commands.getoutput(cmd)


out_lines = []
for fasta_fname in fastas:
	msa_fname = fasta_fname+'.msa'
	inputs = [blast_exe, '-i ', fasta_fname, '-b', max_seq_n, '-d',  blast_db, '-j', blast_iter_n, '-o', msa_fname+' -m 0', '-e', e_thr]
	cmd = ' '.join(inputs)
	print cmd
	commands.getoutput(cmd)

	# 3. msa2block
	block_fname = fasta_fname+'.blk'
	inputs = ['python', bin_dir+'/msa2block.py', msa_fname, fasta_fname, bin_dir, block_fname]
	cmd = ' '.join(inputs)
	print cmd
	commands.getoutput(cmd)

	# 4. msa_entropy.py
	entropy_fname = msa_fname+'.entropy'
	inputs = ['python', bin_dir+'/get_msa_entropy.py', block_fname, entropy_fname]
	cmd = ' '.join(inputs)
	print cmd
	commands.getoutput(cmd)

	# 5. combine
	lines = open(entropy_fname).read().strip().split('\n')
	header, content = [lines[0], lines[1:]]

	out_lines.extend(content)
# end for

# write output
fout = open(final_entropy_fname, 'w')
fout.writelines(header+'\n')
fout.writelines('\n'.join(out_lines)+'\n')
fout.close()


