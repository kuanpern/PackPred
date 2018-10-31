from numpy import *
from modeller import *
log.none()
env = environ()
import sys

score_path = './lib/' # default parameter


class clique_id:
	def __init__(self):
		self.residues = []
		self.ctype = None
		self.order = 0
		self.csize = None
		self.score = None
	# end def
# end class
		

cliquefname = sys.argv[1]
pdbfname = sys.argv[2]
if len(sys.argv) == 4:
	score_path = sys.argv[3]
# end if

# read model
mdl = model(env, file = pdbfname)
restype_dict = dict([[r.chain.name+r.num, r.code] for r in mdl.residues])

# annotate cliques
C = []
clique_lines = open(cliquefname).read().strip().split('\n')
for line in clique_lines:
	bufferline = line.split()
	cres = bufferline[:-1]
	csize = float(bufferline[-1])

	tmp = sorted([[restype_dict[r], r] for r in cres])
	ctype = ''.join([t[0] for t in tmp])
	cres = [t[1] for t in tmp]

	C.append(clique_id())
	C[-1].residues = cres
	C[-1].ctype = ctype
	C[-1].order = len(cres)
	C[-1].csize = csize
# end for

# score cliques
latin_n = {2:'BI', 3:'TRI', 4:'TETRA', 5:'PENTA'}
for c in C:
	fname = score_path+latin_n[c.order]+'CLIQUES/scores/'+c.ctype+'.score'
	try:
		fin = open(fname)
		for line in fin:
			bufferline = line.split()
			c.score = float(bufferline[2])
			if float(bufferline[0]) > c.csize:
				break
			# end if
		# end for
	except:
		c.score = None
	# end try
# end def

# map cliques to residues (somewhat memory intensive here)
residues = dict([[r, []] for r in sorted(restype_dict.keys())])
for c in C:
	if c.score == None:
		print 'warning:', c.residues, c.ctype, 'do not have a valid clique score'
		continue
	# end if

	for r in c.residues:
		residues[r].append(c)
	# end for
# end for

# score for each residues
score_res = {}
for r in residues.keys():
	s = [c.score for c in residues[r] if str(c.score) != 'nan']
#	s = [c.score for c in residues[r]] # no check version
	if s != []: # cases where no cliques is available for a residue
		score_res.update({r: [mean(s), std(s)]})
	# end if
# end for


# map to structure
outlines = ['']
for atom in mdl.atoms:
	r = atom.residue
	r_name = r.chain.name+r.num
	try:
		atom.biso = score_res[r_name][0]
		out = r_name + ' '+ ' '.join([str(round(t, 3)) for t in score_res[r_name]])
		if out != outlines[-1]:
			outlines.append(out)
		# end if
	except KeyError:
		atom.biso = 0
	# end try
# end for

# print output
for line in outlines:
	print line
# end for

# map color to pdb
outfname = '.pdb'.join(pdbfname.split('.pdb')[:-1])+'.clique.'+str(C[-1].order)+'.pdb'
mdl.write(outfname)
