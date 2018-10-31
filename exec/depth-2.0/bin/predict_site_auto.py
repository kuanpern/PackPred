# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# automated script to run binding site prediction program
import sys
import commands

if len(sys.argv) != 4:
	print 'Usage: python predict_site_auto.py input.pdb output parameter_file'
	sys.exit(1)
# end if

pdb_name  = sys.argv[1]
output    = sys.argv[2]
par_fname = sys.argv[3]

# output filenames
cavity_pred_pdb = output+'.pred.pdb'
cavity_pred_out = output+'.pred.out'

# read input parameters
input_labels = {'DEPTH_EXECUTABLE':	0, 'PARAMETER_DIR':	1, 'DETECTION_THRESHOLD':	2, 'CAVITY_SIZE':	3, 'RESOLVATION_CYCLE':	4, 'SOLVENT_SHELL_SIZE':	5, 'DEPTH_CYCLE':	6, 'MIN_SOLVENT_NEIGHBOURS':	7, 'ASA_EXECUTABLE':	8, 'ASA_RESOLUTION':	9, 'ASA_PROBE_RADIUS':	10, 'USE_MSA':	11, 'HOTSPOT_EXECUTABLE':	12, 'BLOSUM62':	13, 'BLAST_DATABASE':	14, 'BLASTPGP_EXECUTABLE':	15, 'BLASTPGP_CYCLE':	16, 'BLASTPGP_E_THRESHOLD':	17, 'BLASTPGP_MAX_SEQ':	18, 'SCRIPTS_DIR':	19}
I = [None]*20
fin = open(par_fname)
for line in fin:
	key, value = line.split()
	try:
		I[input_labels[key]] = value
	except KeyError:
		print 'input label "'+key+'" not recognized. Ignored.'
	# end try
# end for
fin.close()

if None in I:
	print 'Not all input parameters have been specified. Please check input parameter file.'
	print 'Note: all parameters have to be provided:', ' '.join(input_labels.keys())
	sys.exit(1)
# end if

# arrange parameter and run
depth_exe, parameter_dir, threshold, cavity_size, iteration_sol, dist_sol, iteration_depth, survive_n, ASA_exe, resolution, probe_radius, use_msa, hotspot_exe, blosum62, blast_db, blast_exe, blast_iter_n, e_thr, max_seq_n, bin_dir = I
inputs = [bin_dir+'/'+'predict_site.py'] + [pdb_name, depth_exe, ASA_exe, hotspot_exe, parameter_dir, threshold, cavity_size, cavity_pred_pdb, cavity_pred_out, bin_dir, iteration_sol, dist_sol, iteration_depth, survive_n, resolution, probe_radius, use_msa, blosum62, blast_db, blast_exe, blast_iter_n, e_thr, max_seq_n]
cmd = ' '.join(['python'] + inputs)
print cmd
commands.getoutput(cmd)
