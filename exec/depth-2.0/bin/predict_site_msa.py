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

pdb_name        = sys.argv[1] #'1MVT.pdb'
depth_exe       = sys.argv[2] #'DEPTH'
ASA_exe         = sys.argv[3] #'ASA'
hotspot_exe     = sys.argv[4] #'/home/kuanpern/Downloads/server_backup/cgi-bin/predict-binding-site'
parameter_dir   = sys.argv[5] #'/home/kuanpern/Downloads/server_backup/cgi-bin/'
distmatrix_exe  = sys.argv[6] #'/home/kuanpern/Downloads/server_backup/cgi-bin/distmatrix_self'
threshold       = sys.argv[7] #'0.5'
cavity_size     = sys.argv[8] #'4.2'
cavity_pred_pdb = sys.argv[9] #'/home/kuanpern/Downloads/server_backup/cgi-bin/DEPTH'
cavity_pred_out = sys.argv[10]
bin_dir         = sys.argv[11] #'/home/kuanpern/Downloads/server_backup/cgi-bin/'
iteration_sol   = sys.argv[12] #'5'
dist_sol        = sys.argv[13] #'4.2'
iteration_depth = sys.argv[14] #'25'
survive_n       = sys.argv[15] #'4'
resolution      = sys.argv[16] #'92'
probe_radius    = sys.argv[17] #'1.4'

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
# resolvation and get neighbours
outflist = pdb_name+'.sol.distlist'
inputs = ['python', bin_dir+'/get_distlist.py', pdb_name, iteration_sol, cavity_size, depth_exe, distmatrix_exe, dist_sol, outflist]
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)
# classify binding cavity prediction
inputs = ['python', bin_dir+'/cavity_finder.py', outflist, dist_sol, cavity_pred_pdb, threshold, cavity_size, pdb_name+'.binary']
cmd = ' '.join(inputs)
print cmd
commands.getoutput(cmd)
