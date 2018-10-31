#!/usr/bin/python
# mutate cliques
import sys
from copy import deepcopy
import json

if len(sys.argv) != 3:
	print 'python', sys.argv[0], 'input.json mutation_list.json'
	sys.exit(0)
# end if

# read and parse input
infname  = sys.argv[1]
mutfname = sys.argv[2]

mutation_list = json.loads(open(mutfname).read())
mutations = mutation_list['mutations']


# read original clique
Cliques_template = json.load(open(infname))

for mutation in mutations:
	# parse mutation
	chain, mut = mutation.split(':')
	from_res = mut[0]
	to_res   = mut[-1]
	resnum = chain + str(int(mut[1:-1]))

	# actually mutate
	Cliques = deepcopy(Cliques_template)
	for clique in Cliques:
		residues = clique['residues']
		for residue in residues:
			if residue['num'] == resnum:
				assert residue['aatype'] == from_res
				residue['aatype'] = to_res
			# end if
		# end for
		# print 'mutate from', clique['clique type'],
	
		new_clique_type = ''.join(sorted([residue['aatype'] for residue in residues]))
		# print 'to', new_clique_type

		clique['clique type'] = new_clique_type
	# end for


	# write output 
	outfname = 'mutant_'+mutation.replace(':', '_')+'.'+infname
	open(outfname, 'w').writelines(json.dumps(Cliques))
# end for
