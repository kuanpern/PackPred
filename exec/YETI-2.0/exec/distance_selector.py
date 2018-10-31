import sys

fname = sys.argv[1]
outfname = sys.argv[2]
fin = open(fname)
fout = open(outfname, 'w')
for line in fin:
	# A9:ALA:CB	A7:VAL:CB	8.45181
	bufferline = line.split()
	res1 = bufferline[0].split(':')[0]
	res2 = bufferline[1].split(':')[0]
	if res1 == res2: # meaningless to define linkage with self
		continue
	# end if
	else:
		fout.writelines(line)
	# end if
# end for
fin.close()
fout.close()


