from scipy.stats import norm
from numpy import *

def isint(s, warn_non_float = False):
	if s in [inf, -inf, nan]:
		if warn_non_float:
			print ' - warning: non-float representation'
		# end if
		return False
	# end if
	if s - int(s) == 0:
		return True
	else:
		return False
	# end if
# end def

def read_table(fname, FS = None, label = [], check_num = True, check_int = True, header = False):
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

def read_Scores_lib(score_files_dir, d_cuts, order):
	Scores_lib = {}
	depth_levels = []
	for d_cut in d_cuts:
		# read residue clique score at different cut-off distances
		score_fname = score_files_dir + '/clique_score.'+str(order)+'-'+str(d_cut)+'.dat'
		scores_lib = read_table(score_fname, header = True)
		_depth_levels = scores_lib['DEPTH']
		scores_lib = dict(zip(zip(scores_lib['CLIQUE'], scores_lib['DEPTH']), scores_lib['LOG_ODD']))
	
		Scores_lib[d_cut] = scores_lib
		if depth_levels == []:
			depth_levels = set(_depth_levels)
		else:
			if depth_levels != set(_depth_levels):
				print 'error: different score files have different depth levels. Exiting'
				sys.exit(0)
			# end if
		# end if
	# end for
	Scores_lib['depth levels'] = list(depth_levels)

	return Scores_lib
# end def

def score_residue(res, Cliques, d_cuts):
	clique_indexes = res['clique_memberships']
	output = {}
	for d_cut in d_cuts:
		scores = [Cliques[index]['scores'][d_cut] for index in clique_indexes]
		scores = [s for s in scores if s != None]
#		scores = [s for s in scores if not(isnan(s))]
		output[d_cut] = mean(scores)
	# end for

	res['score'] = output
# end def

def score_structure(R, d_cuts):
	output = {}

	# mean score of all residues
	for d_cut in d_cuts:
		S = [res['score'][d_cut] for res in R.values()]
		S = [s for s in S if not(isnan(s))] # remove unscoreables
		output[d_cut] = mean(S)
	# end for

	return output
# end def 

def score_clique(clique, d_cuts, Scores_lib):
	assert type(clique) == dict
	if 'scores' not in clique.keys():
		clique['scores'] = dict([[d_cut, None] for d_cut in d_cuts])
	# end if

	c_size = clique['size']
	c_type = clique['clique type']
	depth_u, depth_s = clique['depth']

	# depth level parameters
	depth_levels = list(Scores_lib['depth levels'])
	depth_cuts = []
	for d in depth_levels:
		cut = [float(t) for t in d.split(',')]
		depth_cuts.append(cut)
	# end for

	# generate pdf for depth
	rv = norm(loc=depth_u, scale=depth_s)

	# score at different level of cut-off distances
	for d_cut in d_cuts:
		if c_size > d_cut: continue

		values = [Scores_lib[d_cut][(c_type, depth_cut)] for depth_cut in depth_levels]
		probs  = empty(len(values))
		for i in range(len(depth_cuts)):
			start_pt, end_pt = depth_cuts[i]
			prob = rv.cdf(end_pt) - rv.cdf(start_pt)
			probs[i] = prob
		# end for
		probs = probs / sum(probs)
		assert abs(1 - sum(probs)) <= 0.05
		score = dot(probs, values)

		# update to appropriate data-holder
		clique['scores'][d_cut] = score
	# end for 
	return clique
# end def
