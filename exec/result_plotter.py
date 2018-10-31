# plot mutation results into a matrix form
import sys
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from modeller import *
log.none() # tell modeller import shut up
env = environ()

# read input
if len(sys.argv) != 5:
	print 'python', sys.argv[0], 'input.tsv input.pdb outfigname.png outerror.txt'
	sys.exit(0)
# end if

scorefname = sys.argv[1]
pdbfname   = sys.argv[2]
outfigname = sys.argv[3]
errorfname = sys.argv[4]

# read model and construct residue list (y-axis)
mdl = model(env, file = pdbfname)
reslist = []
for chain in mdl.chains:
	for res in chain.residues:
		reslist.append(chain.name + ':' + res.num)
	# end for
# end for

output_table = [[nan for i in range(20)] for j in range(len(reslist))]

aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# read score data
score_data = open(scorefname).read().splitlines()[1:] # skip header

errors = []
for dat in score_data:
	mutant, score = dat.split()[:2]
	score = float(score)
	chain, mut = mutant.split(':')
	_from, _num, _to = [mut[0], mut[1:-1], mut[-1]]
	try:
		_num = str(int(_num))
	except:
		pass
	# end try
	res = chain + ':' + _num
	try:
		i = reslist.index(res)
		j = aa_list.index(_to)
		output_table[i][j] = score
	except:
		errors.append(dat)
	# end try
# end for
		
# actually plotting the table
# - setting the figure dimensions
x_width = 20
y_bit_width = 0.90

plt.figure(figsize = (x_width, y_bit_width*len(reslist)))

masked_table = ma.array(output_table, mask = isnan(output_table))
plt.matshow(masked_table)

# setting the x- and y- ticks
plt.xticks(arange(20), aa_list)
plt.yticks(arange(len(reslist)), reslist)
ax = plt.gca()
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

plt.savefig(outfigname, dpi=150, bbox_inches='tight')


# output to error file
if len(errors) != 0:
	ferr = open(errorfname, 'a')
	ferr.writelines('The following mutations cannot be mapped to the heatmap:\n')
	ferr.writelines(str(errors)+'\n')
	ferr.close()
# end if
