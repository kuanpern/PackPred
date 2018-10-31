# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

import time
from numpy import *
from random import random

# MATH #

def normalize(data):
	s = sum(data) + 0.0
	return [t/s for t in data]
# end def

def normalize_dict(s):
	s = dict(s)
	total = sum([s[k] for k in s.keys()]) + 0.0
	return dict([[k, s[k] / total] for k in s.keys()])
# end def

def isint(s):
	if s - int(s) == 0:
		return True
	else:
		return False
	# end if
# end def

def crossing(price_in):
	'''return [index, unit-gradient]'''
	indexed_price = list(enumerate(price_in))
	det = [] # zero is not considered
	for i in range(1, len(indexed_price)):
		if indexed_price[i] != 0:
			det.append(indexed_price[i])
		# end if
	# end for

	output = []
	for i in range(1, len(det)):
		if det[i-1][1] < 0 and det[i][1] > 0:
			output.append([det[i][0], 1])
 		elif det[i-1][1] > 0 and det[i][1] < 0:
			output.append([det[i][0],-1])
		# end if
	# end for
	return output
# end def

def connect_points(pts):
	pts.sort()
	x, y = [[],[]]
	for t in range(len(pts) - 1):
		n = pts[t+1][0] - pts[t][0] + 1
		u = linspace(pts[t][0], pts[t+1][0], n)
		v = linspace(pts[t][1], pts[t+1][1], n)
		x.extend(u[:-1]) # exclude last point, which will be first from next secion
		y.extend(v[:-1])
	# end for
	output = [[p, q] for p,q in zip(x,y)]
	output.append([pts[-1][0], pts[-1][1]]) # include last point
	return output
# end def

# FILE IO #

def read_table(fname, FS = None, label = [], check_num = True, check_int = True, header = True):
	output = []
	fin = open(fname)

	if header == True:
		headers = fin.next().split()
	# end if


	for line in fin:
		output.append(line.split(FS))
	# end for
	fin.close()

	if check_num == False:
		return output
	# end if

	if type(label) == int:
		label = [label]
	# end for
	for i in range(len(output)):
		for j in range(len(output[i])):
			if j not in label:
				try:
					output[i][j] = float(output[i][j])
					if isint(output[i][j]):
						output[i][j] = int(output[i][j])
					# end if
				except ValueError:
					pass
				# end try
			# end if
		# end for
	# end for

	if header == True:
		out_dict = {}
		for i in range(len(headers)):
			content = [output[t][i] for t in range(len(output))]
			out_dict.update({headers[i]:content})
		# end for
		output = out_dict
	# end if

	return output
# end def


def read_dict(fname, FS = None, label = []):
	s = read_table(fname, FS, label)
	return dict([[t[0], t[1:]] for t in s])
# end def

def read_list(fname):
	f = open(fname)
	output = f.readlines()
	for i in range(len(output)):
		output[i] = output[i].strip()
	# end for
	f.close()
	return output
# end def

def print_table(s, fname):
	fout = open(fname, 'w')
	for i in range(len(s)):
		stmp = list(s[i])
		stmp = [str(t) for t in stmp]
		fout.writelines('\t'.join(stmp)+'\n')
	# end for
	fout.close()
# end def

def print_list(s, fname):
	fout = open(fname, 'w')
	for i in range(len(s)):
		fout.writelines(str(s[i])+'\n')
	# end for
	fout.close()
# end def

def print_string(s, fname):
	fout = open(fname, 'w')
	fout.writelines(s+'\n')
	fout.close()
# end def

def unique_name():
	return str(time.time()).replace('.','')
# end def

def shuffle(s):
	return [x[1] for x in sorted([[random(), s[i]] for i in range(len(s))])]

def dict2list(s):
	return [[x,s[x]] for x in sorted(s.keys())]
# end def

def remove_list(s, bag):
	s_tmp = []
	for i in s:
		if i not in bag:
			s_tmp.append(i)
		# end if
	# end for
	return s_tmp
# end def

def collapse_list(s, type=str):
	'''works only on list with numerical list'''
	s = str(s)
	s = s.replace('[', '')
	s = s.replace(']', '')
	s = s.replace(',', '')
	s = s.split()
	if type == int:
		s = [int(t) for t in s]
	elif type == float:
		s = [float(t) for t in s]
	# end if
	return s
# end def

def cluster_pairs(parts):

	D = [parts[0]]
	n_old = 0
	while True:
		for i in range(len(parts)):
			this_part = parts[i]
			used = False
			for i in range(len(D)):
				if intersect(D[i], this_part):
					D[i].extend(this_part)
					D[i] = list(set(D[i]))
					used = True
					break
				# end if
			# end for

			if used == False:
				D.append(this_part)
			# end if
		# end for

		n = len(D)
		if n == n_old or n == 1:
			break
		else:
			parts = D
			D = []
			n_old = n
		# end if
	# end while

	return D
# end def

def intersect(S1, S2):
	return len(set(S1) & set(S2)) > 0
# end def

