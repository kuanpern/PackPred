import sys
import json
from modeller import *
log.none() # asks modeller to shut up, in case we were discovered
env = environ()

# read input
if len(sys.argv) != 3:
	print 'python', sys.argv[0], 'input.json mutation.list'
	sys.exit(0)
# end if

# read input
input_fname  = sys.argv[1]
output_fname = sys.argv[2]

# parse input
inputs = json.loads(open(input_fname).read())
pdb_name      = inputs['pdb name']
mutation_type = inputs['mutation type']
positions     = inputs['positions']
mutations     = inputs['mutations']

# read pdb using modeller
mdl = model(env, file = str(pdb_name))

# assert mutation type
assert mutation_type in ['Manual Input', 'Alanine Scanning Mutagenesis', 'Saturated Scanning Mutagenesis'], 'Error: do not understand the mutation type %s' % mutation_type
output = {'mutations': [], 'errors': []}

# for manual input
matcheds, mismatcheds = [[], []]
if mutation_type == 'Manual Input':
	if len(mutations) == 0:
		print 'Error: no mutation is specified. Please check input.'
		open(output_fname, 'w').writelines('')
		sys.exit(1)
	# end if
	
	# check if the mutations and PDB reconciles
	for mutation in mutations:
		chain_name, mut = mutation.split(':')
		chain = mdl.chains[str(chain_name)]
		_from, _resnum = mut[0], mut[1:-1]
		if _from == chain.residues[str(_resnum)].code:
			matcheds.append(mutation)
		else:
			mismatcheds.append(mutation)
		# end if
	# end for

	if len(matcheds) == 0:
		print 'Error: No valid mutation is specified, exiting program'
		open(output_fname, 'w').writelines('')
		sys.exit(1)
	else:
		output['mutations'] =    matcheds
		output['errors']    = mismatcheds
		open(output_fname, 'w').writelines(json.dumps(output))
		sys.exit(0)
	# end if
# end if

# for 'all' positions, explicitly list all positions
if positions == 'all':
	swaps = []
	for residue in mdl.residues:
		position = residue.chain.name + ':' + residue.num
		swaps.append(position)
	# end for
	positions = swaps
# end if

if mutation_type == 'Alanine Scanning Mutagenesis':
	for residue_name in positions:
		resname = ':'.join(list(reversed(residue_name.split(':')))) # modeller format
		try:
			residue = mdl.residues[str(resname)]
		except KeyError:
			output['errors'].append(residue_name)
			continue
		# end try

#		if residue.code == 'A':
#			continue
#		# end if

		chain_name = residue.chain.name
		mutation = chain_name + ':' + residue.code + residue.num + 'A'
		output['mutations'].append(mutation)
	# end for
	open(output_fname, 'w').writelines(json.dumps(output))
	sys.exit(0)
# end if

if mutation_type == 'Saturated Scanning Mutagenesis':
	aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

	for residue_name in positions:
		resname = ':'.join(list(reversed(residue_name.split(':')))) # modeller format
		try:
			residue = mdl.residues[str(resname)]
		except KeyError:
			output['errors'].append(residue_name)
			continue
		# end try

		chain_name = residue.chain.name
		for _to in aa_list:
#			if residue.code == _to:
#				continue
#			# end if

			mutation = chain_name + ':' + residue.code + residue.num + _to
			output['mutations'].append(mutation)
		# end for
	# end for
	open(output_fname, 'w').writelines(json.dumps(output))
	sys.exit(0)
# end if

