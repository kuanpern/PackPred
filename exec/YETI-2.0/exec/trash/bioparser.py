from modeller import *
log.none()
env = environ()
import commands
from shorthand import *
from biophysics import *

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


std_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

class PDB:	# class to read pdb. Functions are a subset to modeller model class
	def __init__(self, filename, keywords = ["ATOM"], remove_H = True, remove_alt = True, remove_non_std = True):
		self.__atom = []
		self.__serial = []
		self.__name = []
		self.__altLoc = []
		self.__resName = []
		self.__chainID = []
		self.__resSeq  = []
		self.__iCode = []
		self.__x = []
		self.__y = []
		self.__z = []
		self.__occupancy = []
		self.__T = []
		self.__element = []
		self.__charge = []

		pdb = open(filename)
		for line in pdb:
			line = line + ' '*80
			det = line[0:6].strip()
			if det not in keywords:
				continue
			# end if
			self.__atom.append(det)
			self.__serial.append(int(float(line[6:11])))
			self.__name.append(line[12:16].strip())
			self.__altLoc.append(line[16:17].strip())
			self.__resName.append(line[17:20].strip())
			self.__chainID.append(line[21:22].strip())
			self.__resSeq.append(int(float(line[22:26])))
			self.__iCode.append(line[26:27].strip())
			self.__x.append(float(line[30:38]))
			self.__y.append(float(line[38:46]))
			self.__z.append(float(line[46:54]))
			self.__occupancy.append(float(line[54:60]))
			self.__T.append(float(line[60:66]))
			self.__element.append(line[76:78].strip())
			self.__charge.append(line[78:80].strip())
		# end for
		pdb.close()

		# filter entries
		remove_indices = []
		if remove_H == True:
			for i in range(len(self.__element)):
				if self.__element[i] == 'H':
					remove_indices.append(i)
				# end if
			# end for
		# end if

		if remove_alt == True:
			for i in range(len(self.__element)):
				if self.__iCode[i] != "":
					if self.__altLoc[i] not in ("", "A", "1"):
						remove_indices.append(i)
					# end if
				# end if
			# end for
		# end for

		if remove_non_std == True:
			for i in range(len(self.__resName)):
				if self.__resName[i] not in std_aa:
					remove_indices.append(i)
				# end if
			# end for
		# end if

		self.remove(remove_indices)

		self.__aa_length = len(set(self.__resSeq))
		self.__iter_length = len(self.__x)
	# end def __init()

	# get functions


	def aa_length(self):
		return self.__aa_length
	# end def

	def serial(self, i = None):
		if i != None:
			return self.__serial[i]
		else:
			return self.__serial
		# end if
	# end def

	def name(self, i = None):
		if i != None:
			return self.__name[i]
		else:
			return self.__name
		# end if
	# end def

	def altLoc(self, i = None):
		if i != None:
			return self.__altLoc[i]
		else:
			return self.__altLoc
		# end if
	# end def

	def resName(self, i = None):
		if i != None:
			return self.__resName[i]
		else:
			return self.__resName
		# end if
	# end def


	def chainID(self, i = None):
		if i != None:
			return self.__chainID[i]
		else:
			return self.__chainID
		# end if
	# end def

	def resSeq(self, i = None):
		if i != None:
			return self.__resSeq[i]
		else:
			return self.__resSeq
		# end if
	# end def


	def iCode(self, i = None):
		if i != None:
			return self.__iCode[i]
		else:
			return self.__iCode
		# end if
	# end def

	def x(self, i = None):
		if i != None:
			return self.__x[i]
		else:
			return self.__x
		# end if
	# end def

	def y(self, i = None):
		if i != None:
			return self.__y[i]
		else:
			return self.__y
		# end if
	# end def

	def z(self, i = None):
		if i != None:
			return self.__z[i]
		else:
			return self.__z
		# end if
	# end def

	def occupancy(self, i = None):
		if i != None:
			return self.__occupancy[i]
		else:
			return self.__occupancy
		# end if
	# end def

	def T(self, i = None):
		if i != None:
			return self.__T[i]
		else:
			return self.__T
		# end if
	# end def

	def element(self, i = None):
		if i != None:
			return self.__element[i]
		else:
			return self.__element
		# end if
	# end def

	def charge(self, i = None):
		if i != None:
			return self.__charge[i]
		else:
			return self.__charge
		# end if
	# end def

	def sequence(self):
		aa = amino_acid()
		seq = ""
		for i in range(self.__iter_length):
			if self.__name[i] == 'CA':
				seq = seq + aa.code(self.__resName[i])
			# end if
		# end for
		return seq
	# end def

	def centerize(self):
		x_mean, y_mean, z_mean = [mean(p) for p in (self.__x, self.__y, self.__z)]
		self.__x = [p - x_mean for p in self.__x]
		self.__y = [p - y_mean for p in self.__y]
		self.__z = [p - z_mean for p in self.__z]
	# end def

	def remove(self, remove_indices): # remove certain element from memory
		self.__serial   = remove_list(self.__serial   , remove_indices)
		self.__name     = remove_list(self.__name     , remove_indices)
		self.__altLoc   = remove_list(self.__altLoc   , remove_indices)
		self.__resName  = remove_list(self.__resName  , remove_indices)
		self.__chainID  = remove_list(self.__chainID  , remove_indices)
		self.__resSeq   = remove_list(self.__resSeq   , remove_indices)
		self.__iCode    = remove_list(self.__iCode    , remove_indices)
		self.__x        = remove_list(self.__x        , remove_indices)
		self.__y        = remove_list(self.__y        , remove_indices)
		self.__z        = remove_list(self.__z        , remove_indices)
		self.__occupancy  = remove_list(self.__occupancy  , remove_indices)
		self.__T        = remove_list(self.__T        , remove_indices)
		self.__element  = remove_list(self.__element  , remove_indices)
		self.__charge   = remove_list(self.__charge   , remove_indices)
	# end def		

	def write(self, outfilename, remove_indices = []):
		fout = open(outfilename,'w')
		for i in range(len(self.__x)):
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			resSeq  = str(self.__resSeq[i])
			x = str(round(self.__x[i],3))
			y = str(round(self.__y[i],3))
			z = str(round(self.__z[i],3))
			occupancy = str(round(self.__occupancy[i],3))
			T = str(round(self.__T[i],3))
			element = self.__element[i]
			charge = str(self.__charge[i])
			altLoc = self.__altLoc[i]
			chainID = self.__chainID[i]
			iCode = self.__iCode[i]
			string = pdb_line(serial, name, resName, resSeq, x, y, z, occupancy, T, element, charge, altLoc, chainID, iCode)
			fout.writelines(string)
		# end for
		fout.close()
	# end def

	def __repr__(self):
		return 'structure contatins '+str(len(self.sequence()))+' residues, '+str(self.__iter_length)+' atoms'
	# end def

	def __len__(self):
		return self.__iter_length
	# end def
# end class

def pdb_line(serial = "", name = "", resName = "", resSeq  = "", x = "", y = "", z = "", occupancy = "", T = "", element = "", charge = "", altLoc = "", chainID = "", iCode = ""):
	'''serial, name, resName, resSeq , x, y, z, occupancy, T, element, charge, altLoc, chainID, iCode'''

	serial = str(serial)
	resSeq = str(resSeq)
	x = str(x)
	y = str(y)
	z = str(z)
	occupancy = str(occupancy)
	T = str(T)
	charge = str(charge)

	atom = 'ATOM  '
	serial =	serial.		rjust(len(range( 6,11)))
	name = 		name.		ljust(len(range(12,15)))
	altLoc = 	altLoc.		ljust(len(range(16,17)))
	resName = 	resName.	ljust(len(range(17,20)))
	chainID = 	chainID.	ljust(len(range(21,22)))
	resSeq = 	resSeq .	rjust(len(range(22,26)))
	iCode = 	iCode.		ljust(len(range(26,27)))
	x = 		x.		rjust(len(range(30,38)))
	y = 		y.		rjust(len(range(38,46)))
	z = 		z.		rjust(len(range(46,54)))
	occupancy = 	occupancy.	rjust(len(range(54,60)))
	T = 		T.		rjust(len(range(60,66)))
	element = 	element.	ljust(len(range(76,78)))
	charge = 	charge.		rjust(len(range(78,80)))

	line = atom + serial + '  ' + name + altLoc + resName + ' ' + chainID + resSeq + iCode + ' '*3 + x + y + z + occupancy + T + ' '*11 + element + charge
	line = line[:80].strip() + '\n'

	return line
# end def


def pdb2xyz(pdbfname, xyzfname):
	'''write xyz coordinate of pdbfname to xyzfname, label with atom index'''
	mdl_pdb = PDB(pdbfname)
	ftmp_out = open(xyzfname, 'w')
	for i in range(len(mdl_pdb)):
		ftmp_out.writelines(str(mdl_pdb.serial(i))+'\t'+str(mdl_pdb.x(i))+'\t'+str(mdl_pdb.y(i))+'\t'+str(mdl_pdb.z(i))+'\n')
	# end for
	ftmp_out.close()
# end def





# MODELLER hydrogen bond file parser
class h_bond:
	def __init__(self):
		self.donor = None
		self.acceptor = None
		self.donor_residue = None
		self.acceptor_residue = None
		self.type = None
		self.donor_atom = None
		self.acceptor_atom = None
		self.distance = None
		self.angle = None
		self.donor_name = None
		self.acceptor_name = None
		self.line = '???'
	# end def

	def __repr__(self):
		return self.line
	# end def

# end class

class H_bonds:	# class to read hydrogen bond file from MODELLER
	def __init__(self, filename):
		self.bond = []
		self.filename = filename
		fin = open(filename)
		for line in fin:
			if not(line.startswith('#')):
				new_hbond = h_bond()
				new_hbond.donor = line[:4].strip()
				new_hbond.donor_chain = line[6]
				new_hbond.acceptor = line[21:25].strip()
				new_hbond.acceptor_chain = line[27]
				new_hbond.donor_residue = line[8:11].strip()
				new_hbond.acceptor_residue = line[29:32].strip()
				new_hbond.type = line[40:48].strip()
				new_hbond.donor_atom = line[11:17].strip()
				new_hbond.acceptor_atom = line[34:40].strip()
				new_hbond.distance = float(line[52:59].strip())
				new_hbond.angle = float(line[61:].strip())
				new_hbond.donor_name = ':'.join([new_hbond.donor_atom, new_hbond.donor, new_hbond.donor_chain]).strip()
				new_hbond.acceptor_name = ':'.join([new_hbond.acceptor_atom, new_hbond.acceptor, new_hbond.acceptor_chain]).strip()
				new_hbond.line = line[:-1] # skip \n symbol

				self.bond.append(new_hbond)
			# end if
		# end for
		fin.close()
		self.length = len(self.bond)
	# end def

	def __len__(self):
		return self.length
	# end def

	def __repr__(self):
		return 'Hydrogen bond record of PDB file "'+self.filename+'", total of '+str(self.length)+' bonds'
	# end def

# end class



# genethreading parser (for html file)
def GenParser(filename):
	class GenAlign:
		def __init__(self,index):
			self.target = ""
			self.template = ""
			self.target_seq = ""
			self.ss_target = ""
			self.template_seq = ""
			self.ss_template = ""
			self.id = 0
		# end def
	# end class


	fin = open(filename)
	off_set = 9
	index = 0
	aligns = [GenAlign(index)]
	start_read = False
	for line in fin:
		if start_read == False:
			if len(line.strip()) == 0:
				continue
			else:
				start_read = True
			# end if
		# end if
		if line.startswith('Percentage'):
			aligns[-1].id = 0.01*float(line.split()[-1].split("%")[0])
			index = index + 1
			aligns.append(GenAlign(index))
			try:
				for i in range(3): # skip 3 lines
					line = fin.next()
				# end for
			except:
				pass
			# end try
			continue
		# end if

		aligns[-1].ss_template = aligns[-1].ss_template + line[off_set:].strip()
	#	print 'sstemp',line
		line = fin.next()
		aligns[-1].template = line[:off_set].strip()
		aligns[-1].template_seq = aligns[-1].template_seq + line[off_set:].strip()
	#	print 'seqtemp',line
		line = fin.next()
		aligns[-1].target = line[:off_set].strip()
		aligns[-1].target_seq = aligns[-1].target_seq + line[off_set:].strip()
	#	print 'seqtarget',line
		line = fin.next()
		aligns[-1].ss_target = aligns[-1].ss_target + line[off_set:].strip()
	#	print 'sstarget',line
		try:
			for i in range(3): # skip 3 lines
				line = fin.next()
			# end for
		except:
			pass
		# end try
	# end for
	fin.close()
	return aligns[:-1]
# end def


def GOR(seq, exedir = "/home/kuanpern/Downloads/GOR"):
   if exedir.endswith('/'):
       exedir = exedir[:-1]
   # end if
   name = unique_name()

   # preparing input file
   tmp_seq_name = name+'.seq'
   fin = open(tmp_seq_name, 'w')
   header = '!\t'+name
   fin.writelines(header+'\n')
   fin.writelines(seq+'\n')
   fin.close()

   outpredname = unique_name()+'.prediction'
   # run GOR
   commands.getoutput(exedir+'/SOURCE/gorIV -prd '+tmp_seq_name+' -seq '+exedir+'/DATABASE/New_KS.267.seq -obs '+exedir+'/DATABASE/New_KS.267.obs -pro '+outpredname)
   fin = open(outpredname)
   read = False
   ss = []
   prediction = []
   
   # parse prediciton
   for line in fin:
       if line.strip() != "":
           # detect first line of data
           if read == False:
               bufferline = line.split()
               if bufferline[0] == 'SEQ':
                   read = True
                   i = 0
               # end if
           else: # start reading
               i = i + 1
               bufferline = line.split()
               ss.append(bufferline[1])
               prediction.append({'H':float(bufferline[2]),'E':float(bufferline[3]),'C':float(bufferline[4])})
           # end if
       # end if
   # end for
   commands.getoutput('rm '+outpredname+' '+tmp_seq_name)

   return ss, prediction
# end def

def DEPTH(pdbfname, survive_n = 2):
   depth_exe = "/home/kuanpern/Downloads/DEPTH-CLONE-2.8.7/DEPTH" # !!! need to be configured
   cmd = depth_exe+" -i "+pdbfname+" -o "+pdbfname+" -n 5 -survive "+str(survive_n)
   commands.getoutput(cmd)
   depth_pdb_name = pdbfname+"-atomic_depth.pdb" # !!! could be DEPTH version dependent
   mdl = PDB(depth_pdb_name)
   return mdl
# end def

def PSA(pdbfname):
   psa_name = unique_name()
   mdl = model(env, file = pdbfname)
   mdl.write_data('PSA', psa_name) # can be modified here
   record = {}
   fin = open(psa_name+'.sol')
   for line in fin:
       if not(line.startswith('ATOM')):
           continue
       # end if
       sol = float(line[64:])
       index = int(line[6:11].strip())
       record.update({index:sol})
   # end for
   fin.close()
   commands.getoutput('rm '+psa_name+'.sol '+psa_name+'.psa')
   mdl = PDB(pdbfname)
   for i in range(len(mdl)):
       mdl.T()[i] = record[mdl.serial(i)]
   # end for
   return mdl
# end def

def SSM(pdbfname):
   psa_name = unique_name()
   mdl = model(env, file = pdbfname)
   mdl.write_data('SSM', psa_name) # can be modified here
   record = {}
   fin = open(psa_name+'.ssm')
   read = False
   for line in fin:
       if read == False:
           if line.startswith('---'):
               read = True
               continue
           else:
               continue
           # end if
       # end if

       resSeq = int(line[4:9].strip())
       ss = int(line.strip()[-1])

       record.update({resSeq:ss})
   # end for
   fin.close()

   mdl = PDB(pdbfname)
   for i in range(len(mdl)):
       mdl.T()[i] = record[mdl.resSeq(i)]
   # end for
   return mdl
# end def
