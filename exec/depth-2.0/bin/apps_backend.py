# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *

def read2floatarray(s):
	w = []
	for i in range(len(s)):
		try:
			t = float(s[i])
		except ValueError:
			t = None
		# end try
		w.append(t)
	# end for
	return array(w, dtype=float)
# end def

			


def isint(s):
	if s - int(s) == 0:
		return True
	else:
		return False
	# end if
# end def

def read_table(fname, FS = None, label = [], check_num = True, check_int = True):
	output = []
	fin = open(fname)
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
	return output
# end def

class bookmarks:
	def __init__(self):
		self.prefix = '<li> <a href='
		self.suffix = '</a> </li>'
		self.record = []
	# end def

	def add(self, name, title):
		line =  self.prefix+'"#'+name+'"> '+title+self.suffix
		self.record.append(line)
	# end def

	def generate(self):
		return self.record
	# end def

# end class

def matplotlib_1st_nan_bugfix(S):
	# first number cannot be 'NaN', matplotlib bug. Bug fix here
	start = 0
	n = len(S[0])
	for i in range(n):
		if list(set([isnan(s[i]) for s in S])) == [False]:
			start = i
			break
		# end if
	# end for
	W = [s[start:] for s in S]
	return W
# end def

def plot_type_I(y_mean, y_std, xnames, title, xlabel, ylabel, figname, xlim='auto', ylim='auto'):
	plt.cla()
	y_upper = y_mean + y_std
	y_lower = y_mean - y_std

	x = range(len(y_mean))

	x, y_mean, y_upper, y_lower = matplotlib_1st_nan_bugfix([x, y_mean, y_upper, y_lower])


	plt_mean = plt.plot(x, y_mean , 'r-', linewidth=2)
	plt_std  = plt.plot(x, y_upper, 'g:')
	plt_std  = plt.plot(x, y_lower, 'g:')

	plt.legend([plt_mean, plt_std], ["mean", "stdev"])

	x_nums = [int(t) for t in linspace(0, len(y_mean)-1, 10)]  # 10 intervals + last point
	x_tics = [xnames[i] for i in x_nums]
	plt.xticks(x_nums, x_tics, rotation=90)

	if xlim == 'auto':
		plt.xlim(0, len(y_mean)+1)
	else:
		plt.xlim(*zip(xlim))
	# end if

	if ylim == 'auto':
		pass
	elif ylim == 'reverse':
		max_y = max(max(y_upper), max(y_lower))
		min_y = min(min(y_upper), min(y_lower))
		plt.ylim(max_y, min_y)
	# end if

	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.savefig(figname)
# end def

def plot_type_II(y, xnames, hlines, title, xlabel, ylabel, figname, xlim='auto', ylim='auto'):
	plt.cla()

	x = range(len(y))

	x, y = matplotlib_1st_nan_bugfix([x, y])

	plt_y = plt.plot(x, y , 'r-', linewidth=2)

	if hlines != None:
		plt_hlines = []
		hvalues, hlegends = hlines
		for h in hvalues:
			plt_hlines.append(plt.axhline(h))
		# end for

		plt.legend(plt_hlines, hlegends)
	# end if

	x_nums = [int(t) for t in linspace(0, len(y)-1, 10)]  # 10 intervals + last point
	x_tics = [xnames[i] for i in x_nums]
	plt.xticks(x_nums, x_tics, rotation=90)

	if xlim == 'auto':
		plt.xlim(0, len(y)+1)
	else:
		plt.xlim(*zip(xlim))
	# end if

	if ylim == 'auto':
		stdev = std(y)
		plt.ylim(min(y)-stdev/3, max(y)+stdev/3)
	elif ylim == 'reverse':
		plt.ylim(max(y), min(y))
	# end if

	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.savefig(figname)
# end def

class Summary:
	def __init__(self):
		self.title = None
		self.out = []
	# end def

	def add_title(self, title):
		self.title = title
	# end def

	def add(self, pair):
		name, par = pair
		self.out.append('<tr><td>'+name+'</td><td>'+par+'</td><tr>')
	# end def

	def generate(self):
		self.lines = []
		self.lines.append('<table>')
		self.lines.append('<tbody><tr><td><b>'+self.title+'</td></tr>')
		for i in range(len(self.out)):
			self.lines.append('\t'+self.out[i])
		# end for
		self.lines.append('</table>')

		output = '\n'.join(self.lines)
		return output
	# end def
# end class

