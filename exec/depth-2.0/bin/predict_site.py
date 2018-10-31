# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# binding site prediction using depth and accessible area
import sys
import commands

pdb_name        = sys.argv[1] 
depth_exe       = sys.argv[2] #'../DEPTH'
ASA_exe         = sys.argv[3] #'../ASA'
hotspot_exe     = sys.argv[4] #'../predict-binding-site'
parameter_dir   = sys.argv[5] #'../'
threshold       = sys.argv[6] #'0.8'
cavity_size     = sys.argv[7] #'4.2'
cavity_pred_pdb = sys.argv[8] 
cavity_pred_out = sys.argv[9]
bin_dir         = sys.argv[10] #'../'
iteration_sol   = sys.argv[11] #'5'
dist_sol        = sys.argv[12] #'4.2'
iteration_depth = sys.argv[13] #'25'
survive_n       = sys.argv[14] #'4'
resolution      = sys.argv[15] #'92'
probe_radius    = sys.argv[16] #'1.4'
use_msa         = sys.argv[17] #'1'
blosum62        = sys.argv[18] #'../../database/BLOSUM62'
blast_db        = sys.argv[19] #'../../database/uniref90.db'
blast_exe       = sys.argv[20] #'../blast-2.2.22/bin/blastpgp'
blast_iter_n    = sys.argv[21] #'5'
e_thr           = sys.argv[22] #'0.0001'
max_seq_n       = sys.argv[23] #'1000'


# get depth
inputs = [depth_exe, '-i', pdb_name, '-o', pdb_name, '-n', iteration_depth, '-survive', survive_n, '-cavity', cavity_size]
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)

# get accessibility
inputs = [ASA_exe, pdb_name, resolution, probe_radius, pdb_name + '.asa']
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)

# hot-spot prediction
inputs = [hotspot_exe, '-d', pdb_name+'-residue.depth', '-a', pdb_name+'.asa', '-p ', parameter_dir +'/par/:'+survive_n, '-c', pdb_name, '-e',  pdb_name+'.depth-asa', '-o', cavity_pred_out, '-y', cavity_pred_pdb]
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)

if use_msa == '1':
	# backup
	cmd  = 'cp ' + cavity_pred_pdb + ' ' +cavity_pred_pdb+'.backup'

	# get entropy values
	entropy_fname = pdb_name+'.entropies'

	inputs = ['python', bin_dir+'/get_blast_entropy.py ', pdb_name, entropy_fname, blosum62, bin_dir, blast_db, blast_db, blast_exe, blast_iter_n, e_thr, max_seq_n]

	cmd = ' '.join(inputs)
	print cmd
	commands.getoutput(cmd)
	
	# combine entropy values
	new_cavity_pred_pdb = cavity_pred_pdb + '.entropy.pdb'
	inputs = ['python', bin_dir+'/combine_pred_entropy.py ', cavity_pred_pdb, entropy_fname, cavity_pred_pdb]
	cmd = ' '.join(inputs)

	print cmd
	commands.getoutput(cmd)
# end if

outfname = pdb_name+'.binary.pdb'
# classify binding cavity prediction
inputs = ['python', bin_dir+'/site_finder.py', cavity_pred_pdb, 'F', outfname, threshold, '7.0', cavity_size, '6.5', '5.6', bin_dir]
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)

