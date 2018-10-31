# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.


from numpy import *

q = [
[0.0215, 0.0016, 0.0022, 0.0030, 0.0016, 0.0058, 0.0011, 0.0032, 0.0033, 0.0044, 0.0013, 0.0019, 0.0022, 0.0019, 0.0023, 0.0063, 0.0037, 0.0051, 0.0004, 0.0013],
[0.0016, 0.0119, 0.0004, 0.0004, 0.0005, 0.0008, 0.0002, 0.0011, 0.0005, 0.0016, 0.0004, 0.0004, 0.0004, 0.0003, 0.0004, 0.0010, 0.0009, 0.0014, 0.0001, 0.0003],
[0.0022, 0.0004, 0.0213, 0.0049, 0.0008, 0.0025, 0.0010, 0.0012, 0.0024, 0.0015, 0.0005, 0.0037, 0.0012, 0.0016, 0.0016, 0.0028, 0.0019, 0.0013, 0.0002, 0.0006],
[0.0030, 0.0004, 0.0049, 0.0161, 0.0009, 0.0019, 0.0014, 0.0012, 0.0041, 0.0020, 0.0007, 0.0022, 0.0014, 0.0035, 0.0027, 0.0030, 0.0020, 0.0017, 0.0003, 0.0009],
[0.0016, 0.0005, 0.0008, 0.0009, 0.0183, 0.0012, 0.0008, 0.0030, 0.0009, 0.0054, 0.0012, 0.0008, 0.0005, 0.0005, 0.0009, 0.0012, 0.0012, 0.0026, 0.0008, 0.0042],
[0.0058, 0.0008, 0.0025, 0.0019, 0.0012, 0.0378, 0.0010, 0.0014, 0.0025, 0.0021, 0.0007, 0.0029, 0.0014, 0.0014, 0.0017, 0.0038, 0.0022, 0.0018, 0.0004, 0.0008],
[0.0011, 0.0002, 0.0010, 0.0014, 0.0008, 0.0010, 0.0093, 0.0006, 0.0012, 0.0010, 0.0004, 0.0014, 0.0005, 0.0010, 0.0012, 0.0011, 0.0007, 0.0006, 0.0002, 0.0015],
[0.0032, 0.0011, 0.0012, 0.0012, 0.0030, 0.0014, 0.0006, 0.0184, 0.0016, 0.0114, 0.0025, 0.0010, 0.0010, 0.0009, 0.0012, 0.0017, 0.0027, 0.0120, 0.0004, 0.0014],
[0.0033, 0.0005, 0.0024, 0.0041, 0.0009, 0.0025, 0.0012, 0.0016, 0.0161, 0.0025, 0.0009, 0.0024, 0.0016, 0.0031, 0.0062, 0.0031, 0.0023, 0.0019, 0.0003, 0.0010],
[0.0044, 0.0016, 0.0015, 0.0020, 0.0054, 0.0021, 0.0010, 0.0114, 0.0025, 0.0371, 0.0049, 0.0014, 0.0014, 0.0016, 0.0024, 0.0024, 0.0033, 0.0095, 0.0007, 0.0022],
[0.0013, 0.0004, 0.0005, 0.0007, 0.0012, 0.0007, 0.0004, 0.0025, 0.0009, 0.0049, 0.0040, 0.0005, 0.0004, 0.0007, 0.0008, 0.0009, 0.0010, 0.0023, 0.0002, 0.0006],
[0.0019, 0.0004, 0.0037, 0.0022, 0.0008, 0.0029, 0.0014, 0.0010, 0.0024, 0.0014, 0.0005, 0.0141, 0.0009, 0.0015, 0.0020, 0.0031, 0.0022, 0.0012, 0.0002, 0.0007],
[0.0022, 0.0004, 0.0012, 0.0014, 0.0005, 0.0014, 0.0005, 0.0010, 0.0016, 0.0014, 0.0004, 0.0009, 0.0191, 0.0008, 0.0010, 0.0017, 0.0014, 0.0012, 0.0001, 0.0005],
[0.0019, 0.0003, 0.0016, 0.0035, 0.0005, 0.0014, 0.0010, 0.0009, 0.0031, 0.0016, 0.0007, 0.0015, 0.0008, 0.0073, 0.0025, 0.0019, 0.0014, 0.0012, 0.0002, 0.0007],
[0.0023, 0.0004, 0.0016, 0.0027, 0.0009, 0.0017, 0.0012, 0.0012, 0.0062, 0.0024, 0.0008, 0.0020, 0.0010, 0.0025, 0.0178, 0.0023, 0.0018, 0.0016, 0.0003, 0.0009],
[0.0063, 0.0010, 0.0028, 0.0030, 0.0012, 0.0038, 0.0011, 0.0017, 0.0031, 0.0024, 0.0009, 0.0031, 0.0017, 0.0019, 0.0023, 0.0126, 0.0047, 0.0024, 0.0003, 0.0010],
[0.0037, 0.0009, 0.0019, 0.0020, 0.0012, 0.0022, 0.0007, 0.0027, 0.0023, 0.0033, 0.0010, 0.0022, 0.0014, 0.0014, 0.0018, 0.0047, 0.0125, 0.0036, 0.0003, 0.0009],
[0.0051, 0.0014, 0.0013, 0.0017, 0.0026, 0.0018, 0.0006, 0.0120, 0.0019, 0.0095, 0.0023, 0.0012, 0.0012, 0.0012, 0.0016, 0.0024, 0.0036, 0.0196, 0.0004, 0.0015],
[0.0004, 0.0001, 0.0002, 0.0003, 0.0008, 0.0004, 0.0002, 0.0004, 0.0003, 0.0007, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003, 0.0003, 0.0003, 0.0004, 0.0065, 0.0009],
[0.0013, 0.0003, 0.0006, 0.0009, 0.0042, 0.0008, 0.0015, 0.0014, 0.0010, 0.0022, 0.0006, 0.0007, 0.0005, 0.0007, 0.0009, 0.0010, 0.0009, 0.0015, 0.0009, 0.0102]]

aa_range = range(20)
aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Q the probability of amino acid i occuring in the background distribution
Q = [None for i in aa_range]
for i in aa_range:
	Q[i] = sum([q[i][a] for a in aa_range])
# end for

bg = [0.0813, 0.0142, 0.0540, 0.0673, 0.0388, 0.0704, 0.0228, 0.0592, 0.0588, 0.0967, 0.0241, 0.0405, 0.0477, 0.0395, 0.0550, 0.0667, 0.0535, 0.0682, 0.0109, 0.0293]



def MSA_entropy(M, k=1, m=5):
	M = [list(aln) for aln in M]

	'''column-wise entropy of a MSA with consideration of pseudo-count'''
	alignment_n = len(M)
	A = transpose(M).tolist() # transpose alignment for column wise operations
	aa_length   = len(A)

	# get memory
	KL1 = [None for i in range(len(A))] # entropy at position c
	KL2 = [None for i in range(len(A))] # entropy at position c
	H = [None for i in range(len(A))] # entropy at position c
	N = [None for i in range(len(A))] # number of   real count at position c
	B = [None for i in range(len(A))] # number of pseudo count at position c
	R = [None for i in range(len(A))] # number of different residue types at position c

	n = zeros([aa_length, len(aa_range)]) # real   count of each residue type
	b = zeros([aa_length, len(aa_range)]) # pseudo count of each residue type
	p = zeros([aa_length, len(aa_range)]) # probability that amino acid a is in position c


	# for each position
	for c in range(aa_length):
		R[c] = 0.0 + len(set(A[c]) - set(['-'])) # get number of different residue types at position c
		B[c] = 0.0 + m * R[c] # number of pseudo-counts
		N[c] = 0.0 + len(A[c]) - A[c].count('-') # number of real count (all res minus gap)
	# end for

	# real count of amino acid a in column c
	for c in range(aa_length):
		for a in aa_range:
			n[c][a] = 0.0 + A[c].count(aa_list[a])
		# end for
	# end for

	# pseudo count of amino acid a in column c
	for c in range(aa_length):
		for a in aa_range:
			b[c][a] = B[c]*sum([(n[c][a2] / N[c]) * (q[a2][a] / Q[a2]) for a2 in aa_range])
		# end for
	# end for

	# probability that amino acid a is in position c
	for c in range(aa_length):
		for a in aa_range:
			p[c][a] = (N[c] / (N[c] + B[c])) * (n[c][a]/N[c]) + (B[c] / (N[c]+B[c])) * (b[c][a] / B[c])
		# end for
	# end for

	# calculate JS divergence 
	for c in range(aa_length):
		r = [ 0.5*(p[c][a] + bg[a]) for a in aa_range ]
		KL1[c] = k * sum([p[c][a]*log(p[c][a]/r[a]) for a in aa_range])
		KL2[c] = k * sum([  bg[a]*log(  bg[a]/r[a]) for a in aa_range])		
#		KL1[c] = k * sum([p[c][a]*log( p[c][a]/bg[a]  ) for a in aa_range])
#		KL2[c] = k * sum([  bg[a]*log(bg[a]   /p[c][a]) for a in aa_range])
		H[c] = 0.5*(KL1[c] + KL2[c])
	# end for

	H = array(H) / log(2)

	return H, p
# end def


