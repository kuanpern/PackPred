from bioparser import PDB
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

# select atom from pdb file (here select side chain)
mdl = PDB(infile)
remove_list = []
for i in range(len(mdl)):
	if mdl.name()[i] in ('N', 'C', 'O'):
		remove_list.append(i)
	elif mdl.name()[i] == 'CA':
		if mdl.resName()[i] != 'GLY':
			remove_list.append(i)
		# end if
	# end if
# end for
mdl.write(outfile, remove_list)
