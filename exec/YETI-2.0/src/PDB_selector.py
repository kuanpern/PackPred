from bioparser import *
import sys

fin = open(sys.argv[1])
fout = open(sys.argv[2], 'w')

for line in fin:
	name = line[12:16].strip()
	resname = line[17:20].strip()

	if name in ['N', 'C', 'O']:
		continue
	elif name == 'CA':
		if resname != 'GLY':
			continue
		# end if
	# end if

	fout.writelines(line)
# end for
fin.close()
fout.close()
