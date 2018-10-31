from numpy import mean, sqrt, pi
from bioparser import *
from shorthand import *
import commands
from modeller import *
from modeller.scripts import complete_pdb
from modeller.automodel import *
log.none()
env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# radius of gyration
def gyration(x, y, z):
	x_mean, y_mean, z_mean = [mean(x), mean(y), mean(z)]
	gyr = 0.0
	for i in range(len(x)):
		gyr = gyr + (dist(x[i],y[i],z[i],x_mean,y_mean,z_mean))**2
	# end for
	gyr = sqrt(gyr/(len(x)+0.0))
	return gyr
# end def

def dist(x1,y1,z1,x2,y2,z2):
	'''distance between two Cartesian coordinate'''
	return sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
# END def

# STATISTICS DEFINITIONS
def histogram_fixed_size(x, bin_size):	# build a histogram by specifying bin-size
	max_x = max(x) + bin_size / 2.0
	min_x = min(x) - bin_size / 2.0
	bin_n = int((max_x - min_x) / bin_size) + 1
	bins = linspace(min_x, max_x, bin_n)
	raw = histogram(x, bins)
	return [[p, q] for p,q in zip(raw[1], raw[0])]
# end def


def sphericity(sasa, vol):
	return pi**(1/3.0) * (6*vol)**(2/3.0) / sasa
# end def

# MODELLER DEPENDENT
# from bioparser import *

def reindex_pdb(pdbfname, outfname):
	mdl = PDB(pdbfname)
	t = [mdl.resSeq(0)]
	for i in range(1, len(mdl)):
		if mdl.resSeq(i) != t[-1]:
			t.append(mdl.resSeq(i))
		# end if
	# end for
	t = dict([ [a[1], a[0]] for a in enumerate(t)])
	for i in range(len(mdl)):
		mdl.__serial[i] = i + 1
		mdl.__resSeq[i] = t[mdl.resSeq(i)]
	# end for
	mdl.write(outfname)
# end def

def atom_dist(A, B):
	'''distance between two atom (modeller object)'''
	p = array([A.x, A.y, A.z])
	q = array([B.x, B.y, B.z])
	return linalg.norm(p - q)
# end def

def atom_angle(A, B, C):
	p = array([A.x, A.y, A.z])
	q = array([B.x, B.y, B.z])
	r = array([C.x, C.y, C.z])
	M, N = [p - q, r - q]
	v = dot(M, N) / (linalg.norm(M) * linalg.norm(N))
	return arccos(v)
# end def

def is_donor(a):
	return (a.residue.name, a.name) in donor_lib
# end def

def is_acceptor(a):
	return (a.residue.name, a.name) in acceptor_lib
# end def

def find_antecedent(a):
	name = a_acceptor_lib[(a.residue.name, a.name)]
	return a.residue.atoms[name]
# end def

def RMSD(mdl1, mdl2, atom_types = ['N', 'CA', 'C', 'O']):
	tmp1 = unique_name()+'_1.pdb'
	tmp2 = unique_name()+'_2.pdb'
	tmp3 = unique_name()+'.ali'
	mdl1.write(tmp1)
	mdl2.write(tmp2)

	env = environ()
	env.io.atom_files_directory = ['../atom_files']
	mdlA = model(env, file = tmp1)
	mdlB = model(env, file = tmp2)
	aln = alignment(env)
	aln.append_model(mdlA, tmp1, atom_files = tmp1)
	aln.append_model(mdlB, tmp2, atom_files = tmp2)
	aln.align()
	atmsel = selection(mdlA).only_atom_types(' '.join(atom_types))
	r = atmsel.superpose(mdlB, aln)
	commands.getoutput('rm '+tmp1+' '+tmp2+' '+tmp3)

	return r.rms, r.drms, output
# end def

def RMSD_per_atom(pdb_1, pdb_2, outfname, atom_list = ['*'], ignore_loop = False):
	'''return per-atom (main-chain) RMSD for input pdb files'''
	mdl1 = model(env, file = pdb_1)
	mdl2 = model(env, file = pdb_2)
	aln = alignment(env)
	aln.append_model(mdl1, pdb_1, atom_files=pdb_1)
	aln.append_model(mdl2, pdb_2, atom_files=pdb_2)
	aln.align() # auto-alignment

	# align only non-loop parts of structures
	s = selection()
	if ignore_loop == False:
		s.add(mdl1)
	else:
		A, B = [SSM(pdb_1), SSM(pdb_2)]
		non_loop1, non_loop2 = [[], []]
		for i in range(len(A)):
			if A.name(i).strip() != 'CA':
				continue
			name = ':'.join([A.name(i), str(A.resSeq(i)), A.chainID(i)])
			if A.T(i) != 0:
				non_loop1.append(aln[pdb_1].atoms[name].residue.num)
		for i in range(len(B)):
			if B.name(i).strip() != 'CA':
				continue
			name = ':'.join([B.name(i), str(B.resSeq(i)), B.chainID(i)])
			if B.T(i) != 0:
				non_loop2.append(aln[pdb_2].atoms[name].residue.get_aligned_residue(aln[pdb_1]).num)

		non_loop1, non_loop2 = [sorted(list(set(non_loop1))), sorted(list(set(non_loop2)))]
		non_loop = list(set(non_loop1) & set(non_loop2))
 		for a in non_loop:
			s.add(mdl1.residues[a])
	# end if

	atomtypes = ' '.join(atom_list)
	if atomtypes == '*':
		atmsel = s
	else:
		atmsel = s.only_atom_types(atomtypes)
	# end if
	r = atmsel.superpose(mdl2, aln)

	# compute rmsd for each atom pair
	aln = alignment(env)
	aln.append_model(mdl1, pdb_1, atom_files=pdb_1)
	aln.append_model(mdl2, pdb_2, atom_files=pdb_2)

	for a in aln[pdb_2].atoms:
		a.biso = 0

	for t in range(len(aln[pdb_2])):
		template_res = aln[pdb_2].residues[t]
		target_res   = template_res.get_aligned_residue(aln[pdb_1])
		for c in ['N', 'CA', 'C', 'O']:
			d = atom_dist(template_res.atoms[c], target_res.atoms[c])
			aln[pdb_2].residues[t].atoms[c].biso = d
		# end for
	# end for
	aln[pdb_2].write(outfname)
# end def


def get_restraints(pdbfname, outrsrname):
	'''get all restraints (stereochemical and homology derivied) from a PDB file'''
	# step 1. generate self alignment
	root = pdbfname.split('.pdb')[0]
	mdl = model(env, file = pdbfname)
	aln = alignment(env)
	aln.append_model(mdl, root, pdbfname)
	aln.append_model(mdl, root+'.target')
	aln.align()
	alnfname = unique_name()+'.self.aln'
	aln.write(alnfname)
	# step 2. build model using self alignment
	a = automodel(env, alnfile = alnfname, knowns = root, sequence = root+'.target')
	a.starting_model = 1
	a.ending_model   = 1
	a.make()
	# step 3. write out restraints from the automodel module
#	commands.getoutput('rm '+alnfname+' '+root+'.target.*')
	a.restraints.condense()
	a.restraints.write(outrsrname)

	commands.getoutput('rm '+alnfname)
	return a.restraints, root+'.target.B99990001.pdb'
# end def

def mdl2xyz(mdl, xyz_fname):
	fout = open(xyz_fname, 'w')
	for atom in mdl.atoms:
		name = atom.name+':'+atom.residue.num+':'+atom.residue.chain.name
		out = [name, str(round(atom.x, 3)), str(round(atom.y, 3)), str(round(atom.z, 3))]
		fout.writelines('\t'.join(out)+'\n')
	# end for
	fout.close()
# end def

def build_xyz_neighbourhood(xyz_fnames, cutoff, exe_dir=""):
	'''note: make sure the index from xyz_fnames do not repeat, or it leads to confusion'''

	# step 1. read files
	if type(xyz_fnames) == str:
		xyz_fnames = [xyz_fnames, xyz_fnames]
	elif type(xyz_fnames) in (list, tuple):
		if len(xyz_fnames) > 2:
			print 'not support more than 2 xyz files'
			return None
		# end if
	# end if
	output = {}

	tmp_distmap = unique_name()+'.distmap'
	cmd = exe_dir+'g_distmatrix '+xyz_fnames[0]+' '+xyz_fnames[1]+' '+str(cutoff)+' '+tmp_distmap
	print cmd
	commands.getoutput(cmd)

	# build neighbourhood library
	fin = open(tmp_distmap)
	for line in fin:
		bufferline = line.split()
		p, q, d = bufferline
		d = float(d)
		if p not in output.keys():
			output.update({p:[]})
		# end if
		if q not in output.keys():
			output.update({q:[]})
		# end if
		if d <= cutoff:
			output[p].append((d,q))
			output[q].append((d,p))
		# end if
	# end for
	fin.close()
#	commands.getoutput('rm '+tmp_distmap)

	# sort unique
	for k in output.keys():
		output[k] = [[t[1], t[0]] for t in sorted(list(set(output[k])))]
	# end for
	output = dict(output)
	return output
# end def

def build_mdl_neighbourhood(mdls, cutoff, exe_dir=""):
	'''note: make sure the index from pdb_fnames do not repeat, or it leads to confusion'''

	# step 1. read files
	if str(type(mdls)) == "<type 'instance'>":
		mdls = [mdls, mdls]
	elif type(mdls) in (list, tuple):
		if len(mdls) > 2:
			print 'not support more than 2 xyz files'
			return None
		# end if
	# end if
	output = {}

	# step 2: convert to xyz
	xyz_fnames = [unique_name()+'.1.xyz', unique_name()+'.2.xyz']
	mdl2xyz(mdls[0], xyz_fnames[0])
	mdl2xyz(mdls[1], xyz_fnames[1])

	# step 3: compute neighbourhood
	z = build_xyz_neighbourhood(xyz_fnames, cutoff, exe_dir)
	# serial number is integer
	output = {}
	for i in z.keys():
		output.update({int(i):[[int(x[0]), float(x[1])] for x in z[i]]})
	# end for
	return output
# end def


def build_pdb_neighbourhood(pdb_fnames, cutoff, exe_dir=""):
	'''note: make sure the index from pdb_fnames do not repeat, or it leads to confusion'''

	# step 1. read files
	if type(pdb_fnames) == str:
		pdb_fnames = [pdb_fnames, pdb_fnames]
	elif type(pdb_fnames) in (list, tuple):
		if len(pdb_fnames) > 2:
			print 'not support more than 2 xyz files'
			return None
		# end if
	# end if
	output = {}

	# step 2: convert to mdl
	mdls = [PDB(pdb_fnames[0], remove_H = False), PDB(pdb_fnames[1], remove_H = False)]

	# call build_mdl_neighbourhood
	return build_mdl_neighbourhood(mdls, cutoff, exe_dir="")
# end def

# constants
class amino_acid:	# this is an abstract class
	def __init__(self):
		self.__code = {"A":"ALA", "R":"ARG", "D":"ASP", "N":"ASN", "C":"CYS", "E":"GLU", "Q":"GLN", "G":"GLY", "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS", "M":"MET", "F":"PHE", "P":"PRO", "S":"SER", "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL","ALA":"A", "ARG":"R", "ASP":"D", "ASN":"N", "CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}
		self.__KD_table = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
		self.__mass = {"A":71.09, "R":156.19, "D":115.09, "N":114.11, "C":103.15, "E":129.12, "Q":128.14, "G":57.05, "H":137.14, "I":113.16, "L":113.16, "K":128.17, "M":131.19, "F":147.18, "P":97.12, "S":87.08, "T":101.11, "W":186.12, "Y":163.18, "V":99.14}

		self.__surface = {'C': 135.8, 'D': 141.3, 'S': 118.2, 'Q': 180.6, 'K': 204.1, 'I': 177.1, 'P': 140.9, 'T': 142.1, 'F': 199.8, 'A': 109.1, 'G': 81.2, 'H': 183.9, 'E': 174.5, 'L': 180.30273159366294, 'R': 243.0, 'W': 248.1, 'V': 153.3, 'N': 144.9, 'Y': 212.5, 'M': 197.1}
		self.__ASA = {'A': 85.4931, 'R': 201.192, 'N': 123.064, 'D': 121.811, 'C': 109.152, 'Q': 153.715, 'E': 151.585, 'G': 58.5531, 'H': 148.807, 'I': 141.436, 'L': 144.628, 'K': 172.384, 'M': 162.441, 'F': 163.043, 'P': 104.763, 'S': 96.081, 'T': 115.097, 'W': 206.168, 'Y': 174.237, 'V': 115.681}
		self.__ASA_polar_side = {'A': 0.0, 'C': 0.0, 'E': 77.719999999999999, 'D': 60.130000000000003, 'G': 0.0, 'F': 0.0, 'I': 0.0, 'H': 33.18, 'K': 56.880000000000003, 'M': 0.0, 'L': 0.0, 'N': 61.5, 'Q': 83.969999999999999, 'P': 0.0, 'S': 33.560000000000002, 'R': 103.97, 'T': 29.32, 'W': 9.8900000000000006, 'V': 0.0, 'Y': 34.509999999999998}

		self.__pKa = {'A': None, 'C': 8.18, 'D': 3.9, 'E': 4.07, 'F': None, 'G': None, 'H': 6.04, 'I': None, 'K': 10.54, 'L': None, 'M': None, 'N': None,  'P': None, 'Q': None, 'R': 12.48, 'S': None, 'T': None, 'U': 5.73, 'V': None, 'W': None, 'Y': 10.46}
		self.__charge = {'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0, 'G': 0, 'H': 1, 'I': 0, 'K': 1, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1, 'S': 0, 'T': 0, 'U': 0, 'V': 0, 'W': 0, 'Y': 0}

		self.__surface_side = {'C': 98.4, 'D': 103.3, 'S': 80.0, 'Q': 142.8, 'K': 166.3, 'I': 139.8, 'P': 120.3, 'T': 104.2, 'F': 164.3, 'A': 70.9, 'G': 33.0, 'H': 148.2, 'E': 136.7, 'L': 142.4, 'R': 205.2, 'W': 209.5, 'V': 116.0, 'N': 106.9, 'Y': 177.0, 'M': 159.6}

		self.__surface_main = {'C': 37.8, 'D': 37.8, 'S': 37.8, 'Q': 37.8, 'K': 37.8, 'I': 37.8, 'P': 20.8, 'T': 37.8, 'F': 37.8, 'A': 37.8, 'G': 48.0, 'H': 37.8, 'E': 37.8, 'L': 37.8, 'R': 37.8, 'W': 37.8, 'V': 37.8, 'N': 37.8, 'Y': 37.8, 'M': 37.8}

		self.__volume = {"A":88.6, "R":173.4, "D":111.1, "N":114.1, "C":108.5, "E":138.4, "Q":143.8, "G":60.1, "H":153.2, "I":166.7, "L":166.7, "K":168.6, "M":162.9, "F":189.9, "P":112.7, "S":89.0, "T":116.1, "W":227.8, "Y":193.6, "V":140.0}
	# end def

	def code(self,s):
		try:
			return self.__code[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def KD(self,s):
		'''Kite-Dolittle hydrophobicity table'''
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__KD_table[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def mass(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__mass[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def surface(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__surface[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def pKa(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__pKa[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def charge(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__charge[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def


	def accessibility(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return [self.__ASA[s.upper()],self.__ASA_polar_side[s.upper()]]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def surface_side(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__surface[s.upper()] - self.__surface['G']
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def

	def volume(self,s):
		if len(s.strip()) == 3:
			s = self.code(s)
		# end if
		try:
			return self.__volume[s.upper()]
		except KeyError:
			print 'unmatched code name'
			return None
		# end try
	# end def
# end class

def RMSD(pdb_1, pdb_2, atom_list = ['*']):
	'''return RMSD and DRMS for input pdb files'''
	mdl1 = model(env, file = pdb_1)
	mdl2 = model(env, file = pdb_2)
	aln = alignment(env)
	aln.append_model(mdl1, pdb_1, atom_files=pdb_1)
	aln.append_model(mdl2, pdb_2, atom_files=pdb_2)
	aln.align() # auto-alignment
	atomtypes = ' '.join(atom_list)
	if atomtypes == '*':
		atmsel = selection(mdl1)
	else:
		atmsel = selection(mdl1).only_atom_types(atomtypes)
	# end if
	r = atmsel.superpose(mdl2, aln)
	return r.rms, r.drms
# end def

def get_atom_name(A):
	return ':'.join([A.name, A.residue.num, A.residue.chain.name])
# end def


def RMSD_HBOND(pdb_1, pdb_2, outfname):
	'''R.M.S.D. of hydrogen bonded atoms (as specified in the first pdb'''

	mdl1 = model(env, file = pdb_1)
	mdl2 = model(env, file = pdb_2)
	aln = alignment(env)
	aln.append_model(mdl1, pdb_1, atom_files=pdb_1)
	aln.append_model(mdl2, pdb_2, atom_files=pdb_2)
	aln.align() # auto-alignment


	# build hydrogen bond for first pdb
	mdl1.write_data('HBONDS', pdb_1)
	A = H_bonds(pdb_1+'.hbnds')
	select_atoms = [A.bond[i].donor_name for i in range(len(A))] + [A.bond[i].acceptor_name for i in range(len(A))]
	select_atoms = sorted(list(set(select_atoms)))
	# select only H-bonded atoms
	atmsel = selection()
	for atom_name in select_atoms:
		atmsel.add(mdl1.atoms[atom_name])
	# end for

	# compute rmsd for each atom pair
	r = atmsel.superpose(mdl2, aln)
	aln = alignment(env)
	aln.append_model(mdl1, pdb_1, atom_files=pdb_1)
	aln.append_model(mdl2, pdb_2, atom_files=pdb_2)
	aln.align()

	# output colored pdb file
	for a in aln[pdb_2].atoms:
		a.biso = 0
	for t in range(len(aln[pdb_2])):
		template_res = aln[pdb_2].residues[t]
		target_res   = template_res.get_aligned_residue(aln[pdb_1])
		atoms = template_res.atoms
		for atom in atoms:	
			if get_atom_name(atom) in select_atoms: # if in H-bonded atom list
				try:
					d = atom_dist(template_res.atoms[atom.name], target_res.atoms[atom.name])
					template_res.atoms[atom.name].biso = d
				except AttributeError: # misalignment
					pass
				# end try
			# end if
		# end for
	# end for
	aln[pdb_2].write(outfname)

	# output R.M.S.D.
	return r.rms
# end def


donor_lib = [('ALA', 'N'), ('ARG', 'N'), ('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'N'), ('ASN', 'ND2'), ('ASP', 'N'), ('ASX', 'N'), ('ASX', 'ND2'), ('CSS', 'N'), ('CYS', 'N'), ('CYS', 'SG'), ('GLN', 'N'), ('GLN', 'NE2'), ('GLU', 'N'), ('GLX', 'N'), ('GLX', 'NE2'), ('GLY', 'N'), ('HIS', 'N'), ('HIS', 'ND1'), ('HSE', 'N'), ('HSE', 'NE2'), ('HSP', 'N'), ('HSP', 'ND1'), ('HSP', 'NE2'), ('ILE', 'N'), ('LEU', 'N'), ('LYS', 'N'), ('LYS', 'NZ'), ('MET', 'N'), ('PHE', 'N'), ('SER', 'N'), ('SER', 'OG'), ('THR', 'N'), ('THR', 'OG1'), ('TRP', 'N'), ('TRP', 'NE1'), ('TYR', 'N'), ('TYR', 'OH'), ('UNK', 'N'), ('VAL', 'N')]
acceptor_lib = [('ALA', 'O'), ('ARG', 'O'), ('ASN', 'O'), ('ASN', 'OD1'), ('ASP', 'O'), ('ASP', 'OD1'), ('ASP', 'OD2'), ('ASX', 'O'), ('ASX', 'OD1'), ('CSS', 'O'), ('CYS', 'O'), ('GLN', 'O'), ('GLN', 'OE1'), ('GLU', 'O'), ('GLU', 'OE1'), ('GLU', 'OE2'), ('GLX', 'O'), ('GLX', 'OE1'), ('GLY', 'O'), ('HIS', 'NE2'), ('HIS', 'O'), ('HSE', 'ND1'), ('HSE', 'O'), ('HSP', 'O'), ('ILE', 'O'), ('LEU', 'O'), ('LYS', 'O'), ('MET', 'O'), ('PHE', 'O'), ('PRO', 'O'), ('SER', 'O'), ('SER', 'OG'), ('THR', 'O'), ('THR', 'OG1'), ('TRP', 'O'), ('TYR', 'O'), ('TYR', 'OH'), ('UNK', 'O'), ('VAL', 'O')]
a_acceptor_lib = {('ALA', 'O'): 'C', ('ARG', 'O'): 'C', ('ASN', 'OD1'): 'CG', ('ASN', 'O'): 'C', ('ASX', 'OD1'): 'CG', ('ASX', 'O'): 'C', ('ASP', 'OD1'): 'CG', ('ASP', 'OD2'): 'CG', ('ASP', 'O'): 'C', ('CYS', 'O'): 'C', ('CSS', 'O'): 'C', ('GLN', 'OE1'): 'CD', ('GLN', 'O'): 'C', ('GLX', 'OE1'): 'CD', ('GLX', 'O'): 'C', ('GLU', 'OE1'): 'CD', ('GLU', 'OE2'): 'CD', ('GLU', 'O'): 'C', ('GLY', 'O'): 'C', ('UNK', 'O'): 'C', ('HIS', 'NE2'): 'CD2', ('HIS', 'O'): 'C', ('HSE', 'ND1'): 'CG', ('HSE', 'O'): 'C', ('HSP', 'O'): 'C', ('ILE', 'O'): 'C', ('LEU', 'O'): 'C', ('LYS', 'O'): 'C', ('MET', 'O'): 'C', ('PHE', 'O'): 'C', ('PRO', 'O'): 'C', ('SER', 'OG'): 'CB', ('SER', 'O'): 'C', ('THR', 'OG1'): 'CB', ('THR', 'O'): 'C', ('TRP', 'O'): 'C', ('TYR', 'OH'): 'CZ', ('TYR', 'O'): 'C', ('VAL', 'O'): 'C'}
acceptor_sc_lib = dict([('ASN', ['OD1']), ('ASP', ['OD1', 'OD2']), ('GLN', ['OE1']), ('GLU', ['OE1', 'OE2']), ('HIS',[ 'NE2']), ('SER', ['OG']), ('THR', ['OG1']), ('TYR', ['OH'])])
donor_sc_lib = dict([('ARG', ['NE', 'NH1']), ('ARG', ['NH2']), ('ASN', ['ND2']), ('CYS', ['SG']), ('GLN', ['NE2']), ('HIS', ['ND1']), ('LYS', ['NZ']), ('SER', ['OG']), ('THR', ['OG1']), ('TRP', ['NE1']), ('TYR', ['OH'])])

vdw = {'H':1.20, 'C':1.70, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.80, 'S':1.80, 'CL':1.75, 'CU':1.4}
chem_group = {('E', 'OE1'): ('CO', 'O'), ('F', 'CE2'): ('RS', 'C'), ('D', 'OD1'): ('CO', 'O'), ('F', 'CG'): ('RS', 'C'), ('I', 'CG2'): ('HO', 'C'), ('M', 'CG'): ('TE', 'C'), ('H', 'CE1'): ('HS', 'C'), ('K', 'CD'): ('H2', 'C'), ('T', 'CB'): ('OL', 'C'), ('Y', 'CE1'): ('RS', 'C'), ('H', 'ND1'): ('HS', 'N'), ('M', 'CB'): ('H1', 'C'), ('N', 'CG'): ('AD', 'C'), ('L', 'CG'): ('HO', 'C'), ('I', 'CG1'): ('HO', 'C'), ('L', 'CB'): ('HO', 'C'), ('Y', 'CZ'): ('RS', 'C'), ('W', 'CZ2'): ('RS', 'C'), ('Q', 'CD'): ('AD', 'C'), ('R', 'NE'): ('GS', 'N'), ('W', 'NE1'): ('RS', 'N'), ('W', 'CE2'): ('RS', 'C'), ('L', 'CD2'): ('HO', 'C'), ('M', 'CE'): ('TE', 'C'), ('I', 'CD1'): ('HO', 'C'), ('H', 'NE2'): ('HS', 'N'), ('R', 'NH1'): ('GS', 'N'), ('Q', 'CG'): ('AD', 'C'), ('E', 'OE2'): ('CO', 'O'), ('Q', 'OE1'): ('AD', 'O'), ('W', 'CG'): ('RS', 'C'), ('E', 'CD'): ('CO', 'C'), ('H', 'CB'): ('HS', 'C'), ('V', 'CB'): ('HO', 'C'), ('K', 'CE'): ('AS', 'C'), ('P', 'CG'): ('H2', 'C'), ('Y', 'CG'): ('RS', 'C'), ('A', 'CB'): ('HO', 'C'), ('V', 'CG2'): ('HO', 'C'), ('Y', 'CE2'): ('RS', 'C'), ('N', 'CB'): ('AD', 'C'), ('W', 'CD1'): ('RS', 'C'), ('D', 'CB'): ('CO', 'C'), ('Q', 'CB'): ('H1', 'C'), ('K', 'NZ'): ('AS', 'N'), ('I', 'CB'): ('HO', 'C'), ('R', 'CD'): ('GS', 'C'), ('W', 'CB'): ('H1', 'C'), ('E', 'CG'): ('CO', 'C'), ('F', 'CD2'): ('RS', 'C'), ('P', 'CD'): ('H2', 'C'), ('Y', 'CB'): ('H1', 'C'), ('N', 'ND2'): ('AD', 'N'), ('E', 'CB'): ('H1', 'C'), ('C', 'CB'): ('TO', 'C'), ('W', 'CD2'): ('RS', 'C'), ('D', 'CG'): ('CO', 'C'), ('K', 'CG'): ('H2', 'C'), ('M', 'SD'): ('TE', 'S'), ('S', 'CB'): ('OL', 'C'), ('K', 'CB'): ('H2', 'C'), ('P', 'CB'): ('H2', 'C'), ('F', 'CD1'): ('RS', 'C'), ('C', 'SG'): ('TO', 'S'), ('Y', 'OH'): ('PH', 'O'), ('W', 'CE3'): ('RS', 'C'), ('R', 'CZ'): ('GS', 'C'), ('H', 'CD2'): ('HS', 'C'), ('W', 'CZ3'): ('RS', 'C'), ('W', 'CH2'): ('RS', 'C'), ('S', 'OG'): ('OL', 'O'), ('Y', 'CD2'): ('RS', 'C'), ('F', 'CE1'): ('RS', 'C'), ('V', 'CG1'): ('HO', 'C'), ('L', 'CD1'): ('HO', 'C'), ('F', 'CB'): ('H1', 'C'), ('R', 'CB'): ('H2', 'C'), ('D', 'OD2'): ('CO', 'O'), ('F', 'CZ'): ('RS', 'C'), ('T', 'CG2'): ('OH', 'C'), ('Y', 'CD1'): ('RS', 'C'), ('Q', 'NE2'): ('AD', 'N'), ('R', 'CG'): ('H2', 'C'), ('H', 'CG'): ('HS', 'C'), ('N', 'OD1'): ('AD', 'O'), ('R', 'NH2'): ('GS', 'N'), ('T', 'OG1'): ('OL', 'O')}
