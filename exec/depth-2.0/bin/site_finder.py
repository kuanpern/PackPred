# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

from PDB import *
import commands
import sys
from shorthand import *

def separate_molsol(solvated_fname):
	root = solvated_fname.split('.pdb')[0]+'.pdb'
	molfname = root + '.mol.pdb'
	solfname = root + '.sol.pdb'
	fout_mol = open(molfname, 'w')
	fout_sol = open(solfname, 'w')
	fin = open(solvated_fname)
	for line in fin:
		if 'SOL' in line:
			fout_sol.writelines(line)
		else:
			fout_mol.writelines(line)
		# end if
	# end for
	fin.close()
	fout_mol.close()
	fout_sol.close()
	return molfname, solfname
# end def

def assign_T(mdl, D):
	for i in range(len(mdl)):
		res = str(mdl.resSeq(i)) + ':' + mdl.chainID(i)
		try:
			mdl.T()[i] = D[res]
		except KeyError:
			pass
		# end try
	# end for
	return mdl
# end def

fname          = sys.argv[1]
write_water    = sys.argv[2] # dummy
outfname       = sys.argv[3]
prob_cutoff    = float(sys.argv[4]) # 0.75
dist_thr       = float(sys.argv[5]) # 7.0
d_thr1         = float(sys.argv[6]) # 4.2
d_thr2         = float(sys.argv[7]) # 6.5
d_thr3         = float(sys.argv[8]) # 5.6
exedir         = sys.argv[9]

# solvate model
cmd = exedir+'DEPTH -i ' + fname + ' -o ' + fname +' -keep '+ fname+'.sol -n 1 -survive 1'
print cmd
commands.getoutput(cmd)
solvated_fname = fname+'.sol-0.pdb'


molsol_distmap_name = fname.split('.pdb')[0]+'.pdb.molsol.distmap'
solsol_distmap_name = fname.split('.pdb')[0]+'.pdb.solsol.distmap'

files = commands.getoutput('ls *distmap').split('\n')
if not(molsol_distmap_name in files) or not(solsol_distmap_name in files): # compute if not exists
	print 'computing distmaps:', molsol_distmap_name, solsol_distmap_name
	# get distance matrix
	molfname, solfname = separate_molsol(solvated_fname)
	cmd = 'python ' + exedir + 'pdb2xyz.py '+molfname
	print cmd
	commands.getoutput(cmd)
	cmd = 'python ' + exedir + 'pdb2xyz.py '+solfname
	print cmd
	commands.getoutput(cmd)


	cmd = exedir + 'g_distmatrix '+molfname+'.xyz'+' '+solfname+'.xyz '+str(dist_thr)+' '+molsol_distmap_name
	print cmd
	commands.getoutput(cmd)
	cmd = exedir + 'g_distmatrix '+solfname+'.xyz'+' '+solfname+'.xyz '+str(dist_thr)+' '+solsol_distmap_name
	print cmd
	commands.getoutput(cmd)
else:
	print 'distmaps already exist'
# end if

# assign neighbours 
z = read_table(molsol_distmap_name, header=False)
print len(z)
D = {}
for entry in z:
	try:
		mol_name, sol_name, d = entry
	except:
		print entry
		sys.exit(1)
	# end try
	mol = ':'.join(mol_name.split(':')[1:])
	sol = ':'.join(sol_name.split(':')[1:])
	try:
		D[(mol, sol)].append(d)
	except KeyError:
		D.update({(mol, sol):[d]})
	# end try
# end for

mol_neighbours, sol_neighbours = [{}, {}]
for key in D.keys():
	mol, sol = key
	d = min(D[key]) # (use average distance)

	if d <= d_thr2:
		try:
			mol_neighbours[sol].append(mol)
		except KeyError:
			mol_neighbours.update({sol:[mol]})
		# end try
	# end if
	if d <= d_thr1:
		try:
			sol_neighbours[mol].append(sol)
		except KeyError:
			sol_neighbours.update({mol:[sol]})
		# end try
	# end if
# end for

for key in mol_neighbours:
	mol_neighbours[key] = list(set(mol_neighbours[key]))
# end for
	
for key in sol_neighbours:
	sol_neighbours[key] = list(set(sol_neighbours[key]))
# end for

# Assign solsol neighbours
# assign neighbours 
z = read_table(solsol_distmap_name, header=False)
D = {}
for entry in z:
	sol1_name, sol2_name, d = entry
	sol1 = ':'.join(sol1_name.split(':')[1:])
	sol2 = ':'.join(sol2_name.split(':')[1:])
	try:
		D[(sol1, sol2)].append(d)
	except KeyError:
		D.update({(sol1, sol2):[d]})
	# end try
# end for

solsol_neighbours = {}
for key in D.keys():
	sol1, sol2 = key
	d = min(D[key]) # (use min distance)

	if d <= d_thr3:
		try:
			solsol_neighbours[sol1].append(sol2)
		except KeyError:
			solsol_neighbours.update({sol1:[sol2]})
		# end try
		try:
			solsol_neighbours[sol2].append(sol1)
		except KeyError:
			solsol_neighbours.update({sol2:[sol1]})
		# end try
	# end if
# end for

for key in solsol_neighbours:
	solsol_neighbours[key] = list(set(solsol_neighbours[key]))
# end for
	
solvents = mol_neighbours.keys()
residues = sol_neighbours.keys()

mdl = PDB(fname)
Pred = dict([[str(mdl.resSeq(i)) + ':' + mdl.chainID(i), max(mdl.T()[i], 0)] for i in range(len(mdl))])

def separate_by_chains(S):
	keys = list(set([res.split(':')[-1] for res in S]))
	chains = dict([ [key, []] for key in keys ])
	for res in S:
		key = res.split(':')[-1]
		chains[key].append(res)
	# end for
	return chains
# end def

# Assign primary probability measure to solvent
Prob_sol = {}
for solvent in solvents:
	neighbour_residues = mol_neighbours[solvent]
	# separate by chains
	chains = separate_by_chains(neighbour_residues)
	P_chains = dict([[key, 0] for key in chains])

	for chain in chains:
		probs = [Pred[res] for res in chains[chain]]
		P = 1 - product([1 - p for p in probs])

		if len(probs) <= 1: # binding is cooperative event, needs more than 1 residues
			P = 0
		# end if
		P_chains[chain] = P
	# end for

	# product of all probability
	P = product([P_chains[key] for key in P_chains.keys()])
	Prob_sol.update({solvent:P})
# end for


# dynamic cut-off implemented here
Probs = [Prob_sol[t] for t in Prob_sol.keys()]
prob_sorted = sorted(Probs)
min_p = mean(prob_sorted[:5])
max_p = mean(prob_sorted[-5:])
data_range = max_p - min_p
for key in Prob_sol.keys():
	Prob_sol[key] = (Prob_sol[key] - min_p) / data_range
	Prob_sol[key] = max(min(Prob_sol[key], 1), 0)
# end for

# average probability among solvent neighbours
hot_solvents = []
for solvent in solvents:
	try:
		probs = [Prob_sol[sol_nearby] for sol_nearby in solsol_neighbours[solvent]]
	except KeyError: # no neighbour at all
		probs = [0] * 2
	# end try

	P = mean(sorted(probs)[-2:]) # arithmetric mean of largest two
	if P >= prob_cutoff:
		hot_solvents.append(solvent)
	# end if
# end for

# find nearby residues
hot_residues = []
[hot_residues.extend(mol_neighbours[sol]) for sol in hot_solvents]
hot_residues = list(set(hot_residues))

# write value to PDB
mdl2 = PDB(solvated_fname, keywords = ['ATOM', 'HETATM'], remove_non_std = False, remove_H = False)
for i in range(len(mdl2)):
	if mdl2.resName(i) == 'SOL':
		mdl2.T()[i] = 0
	# end if
# end for

mdl2 = assign_T(mdl2, Prob_sol)
for i in range(len(mdl2)):
	if mdl2.resName(i) == 'SOL':
		continue
	# end if

	key = str(mdl2.resSeq(i)) + ':' + mdl2.chainID(i)
	if key in hot_residues:
		mdl2.T()[i] = 1
	else:
		mdl2.T()[i] = 0
	# end if
# end for

print 'writing output', outfname
print 'current directory', commands.getoutput('pwd')

if write_water[0].upper() == 'T':
	mdl2.write(outfname, write_water = True)
else:
	mdl2.write(outfname, write_water = False)
# end if
