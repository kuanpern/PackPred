# append information to a clique, save in pickle format
import cPickle as pickle
#from shorthand import *
from modeller import *
log.none()
env = environ()
from numpy import *

import sys

def isint(s):
	if s - int(s) == 0:
		return True
	else:
		return False
	# end if
# end def

def read_table(fname, FS = None, label = [], check_num = True, check_int = True, header = True):
	output = []
	fin = open(fname)

	if header == True:
		headers = fin.next().split()
	# end if


	for line in fin:
		output.append(line.split(FS))
	# end for
	fin.close()

	if check_num == False:
		return output
	# end if

	if type(label) == int:
		label = [label]
	# end for
	for i in range(len(output)):
		for j in range(len(output[i])):
			if j not in label:
				try:
					output[i][j] = float(output[i][j])
					if isint(output[i][j]):
						output[i][j] = int(output[i][j])
					# end if
				except ValueError:
					pass
				# end try
			# end if
		# end for
	# end for

	if header == True:
		out_dict = {}
		for i in range(len(headers)):
			content = [output[t][i] for t in range(len(output))]
			out_dict.update({headers[i]:content})
		# end for
		output = out_dict
	# end if

	return output
# end def


if len(sys.argv) != 5:
	print 'python', sys.argv[0], 'input.clique input.pdb input.residue-depth Cliques.clique.json'
	sys.exit(0)
# end if
clique_fname = sys.argv[1]
pdb_fname    = sys.argv[2]
depth_fname  = sys.argv[3]
out_fname    = sys.argv[4]
   
# read residue amino acid type
mdl = model(env, file = pdb_fname)
aa_dict = dict([[r.chain.name + r.num, r.code] for r in mdl.residues])

# read residue depth
depths = read_table(depth_fname, header=False)
depths = depths[1:] # skip header
depths = [[dat[0].replace(':', ''), [dat[2], dat[3]]] for dat in depths]
depths = dict(depths)

# memory for Residues 
Residues = {}
for resname in depths.keys():
	Residues[resname] = {
		'depth': depths[resname],
		'clique_memberships': [],
		'num': resname,
		'aatype': aa_dict[resname],
	} # end Residue
# end for

# memory for Cliques

Cliques = []
data = read_table(clique_fname, check_num = False)
for i in range(len(data)):
	dat = data[i]
	nums, csize = [dat[:-1], float(dat[-1])]
	order = len(nums)

	# update child residues' membership
	for resname in nums:
		Residues[resname]['clique_memberships'].append(i)
	# end for

	residues, aatypes, depth_data = [[], [], []]
	for res in nums:
		aatype = aa_dict[res]
		aatypes.append(aatype)

		depth  = depths[res]
		depth_data.append(depth)
	# end for

	# clique type
	clique_type = ''.join(sorted(aatypes))

	# clique depth
	U, S = zip(*depth_data)
	u_depth, s_depth = [mean(U), sqrt(mean(array(S)**2))]
	u_depth, s_depth = [round(t, 2) for t in [u_depth, s_depth]]

	_clique = {
		'size'    : csize,
		'order'   : order, 
		'depth'   : [u_depth, s_depth],
		'residues': nums,
		'clique type': clique_type,
	} # end _clique
	Cliques.append(_clique)
# end for

Data = {
	'Residues': Residues,
	'Cliques' : Cliques,
}

# save to pickle format output
pickle.dump(Data, open(out_fname, 'w'), pickle.HIGHEST_PROTOCOL)

