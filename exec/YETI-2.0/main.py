# - the script is modified to better incorporated into automation
# main program of clique building
import sys
import commands
install_dir = "/home/tankp/githubs/packpred_server/exec/YETI-2.0"

def line_count(fname):
	fin = open(fname)
	n = 0
	for line in fin:
		n = n + 1
	# end for
	return n
# end def 

# naming

if len(sys.argv) != 7:
	print 'python main.py input.pdb order(max) cutoff distance_selector linkage_builder linkage_par_file'
	sys.exit(1)
# end if

# read input
pdb_in            = sys.argv[1]
order             = sys.argv[2]
cutoff            = sys.argv[3]
distance_selector = sys.argv[4]
linkage_builder   = sys.argv[5]
linkage_par_file  = sys.argv[6]

# parse input and generate internal parameter names
order = int(order)
cutoff = float(cutoff)

root = pdb_in.split('/')[-1].split('.pdb')[0]
pdb_select = root+'_select.pdb'
pdb_xyz = root+'.xyz'
pdb_distmap_pre = root+'.pre.distmap'
pdb_distmap = root+'.distmap'
pdb_distmap_sorted = root+'.sorted.distmap'
pdb_link = root+'.link'
clique_root = root+'.clique'


# 1. select specific atom from PDB
atom_selector = install_dir+'/exec/PDB_selector.py'
cmd = 'python '+atom_selector+' '+pdb_in+' '+pdb_select
print cmd
commands.getoutput(cmd)
# 2. compute distances

# 2.1 convert PDB coordinate to xyz
xyz_generator = install_dir+'/exec/xyz_generator.py'
cmd = 'python '+xyz_generator+' '+pdb_select+' '+pdb_xyz
print cmd
commands.getoutput(cmd)
# 2.2 invoke cell-list algorithm
cell_list_exe = install_dir+'/exec/g_distmatrix'
cmd = cell_list_exe+' '+pdb_xyz+' '+pdb_xyz+' '+str(cutoff)+' '+pdb_distmap_pre
print cmd
commands.getoutput(cmd)

# 3. build primary linkage
# 3.1 select distances
distance_selector = install_dir+'/exec/distance_selector.py'
cmd = 'python '+distance_selector+' '+pdb_distmap_pre+' '+pdb_distmap
print cmd
commands.getoutput(cmd)
# 3.2 sort distances according to distance magnitude
distmap_sorter = install_dir+'/exec/sort_distmatrix'
cmd = distmap_sorter+' '+pdb_distmap+' '+pdb_distmap_sorted
print cmd
commands.getoutput(cmd)
# 3.3 build linkage
linkage_builder = install_dir+'/exec/linkage'
cmd = linkage_builder+' '+pdb_distmap_sorted+' '+linkage_par_file+' '+pdb_link
print cmd
commands.getoutput(cmd)
# 4. build clique
# 4.1 specify highest order to build, determine filenames
clique_fnames = dict([[i, clique_root+str(i)] for i in range(2, order+1)])
# 4.2 second order clique = linkage
cmd = 'cp '+pdb_link+' '+clique_root+'2'
print cmd
commands.getoutput(cmd)
# 4.2 iteratively build clique
clique_builder = install_dir+'/exec/build_clique'
for i in range(2, order):
	cmd = clique_builder+' '+pdb_link+' '+clique_fnames[i]+' '+clique_root
	print cmd
	commands.getoutput(cmd)
	if line_count(clique_fnames[i]) == 0: # no more clique
		break
	# end if
# end for
