import sys

def remove_list(s, bag):
        s_tmp = []
        for i in s:
                if i not in bag:
                        s_tmp.append(i)
                # end if
        # end for
        return s_tmp
# end def

std_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
class PDB:	# class to read pdb. Functions are a subset to modeller model class
	def __init__(self, filename, keywords = ["ATOM"], remove_H = True, remove_alt = True, remove_non_std = True, remove_water = True):
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
			self.__resSeq.append(int(line[22:26]))
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

		if remove_water == True:
			for i in range(len(self.__resName)):
				if self.__resName[i] in ['HOH', 'DOD']:
					remove_indices.append(i)
				# end if
			# end for
		# end if

		remove_indices = sorted(list(set(remove_indices)))
		self.remove(remove_indices)

		self.__aa_length = len(set(self.__resSeq))
		self.__iter_length = len(self.__x)
	# end def __init()

	# get functions


	def aa_length(self):
		return self.__aa_length
	# end def

	def atom(self, i = None):
		if i != None:
			return self.__atom[i]
		else:
			return self.__atom
		# end if
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
		self.__atom     = remove_list(self.__atom    , remove_indices)
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

fname = sys.argv[1]
outfname = sys.argv[2]
mdl = PDB(fname)

fout = open(outfname, 'w')
for i in range(len(mdl)):
	string = mdl.chainID(i)+str(mdl.resSeq(i))+':'+mdl.resName(i)+':'+mdl.name(i)+'\t'+str(mdl.x(i))+'\t'+str(mdl.y(i))+'\t'+str(mdl.z(i))
	fout.writelines(string+'\n')
# end for
fout.close()
