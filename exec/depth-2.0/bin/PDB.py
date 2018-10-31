# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

def remove_list(s,remove_indices):
	output = []
	for i in range(len(s)):
		if i in remove_indices:
			continue
		# end if
		output.append(s[i])
	# end for
	return output
# end def

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
		self.__serial	= remove_list(self.__serial	, remove_indices)
		self.__name	 = remove_list(self.__name	 , remove_indices)
		self.__altLoc	= remove_list(self.__altLoc	, remove_indices)
		self.__resName  = remove_list(self.__resName  , remove_indices)
		self.__chainID  = remove_list(self.__chainID  , remove_indices)
		self.__resSeq	= remove_list(self.__resSeq	, remove_indices)
		self.__iCode	= remove_list(self.__iCode	, remove_indices)
		self.__x		= remove_list(self.__x		, remove_indices)
		self.__y		= remove_list(self.__y		, remove_indices)
		self.__z		= remove_list(self.__z		, remove_indices)
		self.__occupancy  = remove_list(self.__occupancy  , remove_indices)
		self.__T		= remove_list(self.__T		, remove_indices)
		self.__element  = remove_list(self.__element  , remove_indices)
		self.__charge	= remove_list(self.__charge	, remove_indices)
	# end def		

	def write(self, outfilename, remove_indices = [], write_water = True):
		fout = open(outfilename,'w')
		for i in range(len(self.__x)):
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			if write_water == False:
				if resName in ['SOL', 'HOH', 'DOD']:
					continue
				# end if
			# end if
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


class Atom:
	def __init__(self):
		self.x, self.y, self.z = [None, None, None]
		self.residue = None
		self.name = None
	# end def
# end class

class residue:
	def __init__(self):
		self.chain = None
		self.name = None
		self.atoms = {}
		self.chain = chain()
	# end def
# end class

class chain:
	def __self__(self):
		self.name = ''
	# end def
# end class

class model:
	def __init__(self, pdb_name, keywords = ["ATOM"], remove_H = True, remove_alt = True, remove_non_std = True):
		self.mdl = PDB(pdb_name, keywords, remove_H, remove_alt, remove_non_std)
		self.atoms, self.residues = [[], {}]

		for i in range(len(self.mdl)):
			residue_id = str(self.mdl.resSeq(i))+':'+self.mdl.chainID(i)
			if residue_id not in self.residues.keys():
				new_residue = residue()
				new_residue.name = self.mdl.resName(i)
				new_residue.num  = str(self.mdl.resSeq(i))
				new_residue.chain.name = self.mdl.chainID(i)
				self.residues.update({residue_id:new_residue})
			# end if

			new_atom = Atom()
			new_atom.name = self.mdl.name(i)
			new_atom.x, new_atom.y, new_atom.z = [self.mdl.x(i), self.mdl.y(i), self.mdl.z(i)] 
			self.residues[residue_id].atoms.update({new_atom.name:new_atom})
			new_atom.residue = self.residues[residue_id]

			self.atoms.append(new_atom)
		# end for
	# end def
# end class
			






