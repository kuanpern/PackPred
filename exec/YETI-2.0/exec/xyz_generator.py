from bioparser import *
import sys

fname = sys.argv[1]
outfname = sys.argv[2]
mdl = PDB(fname)

fout = open(outfname, 'w')
for i in range(len(mdl)):
	string = mdl.chainID(i)+str(mdl.resSeq(i))+str(mdl.iCode(i))+':'+mdl.resName(i)+':'+mdl.name(i)+'\t'+str(mdl.x(i))+'\t'+str(mdl.y(i))+'\t'+str(mdl.z(i))
	fout.writelines(string+'\n')
# end for
fout.close()
